# Project notes and decisions

## Step 1: Background concepts

**1. What does one number in the count table represent?**
It is the amount of reads that matched to one gene in a single sample. It measures the amount of RNA from the specific gene in the tissue. A high number correlates to higher gene activation. 

**2. If Gene A has 500 counts and Gene B has 50 in the same sample, is Gene A “10x more active”?**
Realistically, gene activity does not rely on the length of the genes. Instead, 

**3. Why is it fine to compare the same gene across two samples?**
Because the gene’s length is identical in both samples, so it cancels out. The only thing you need to correct for is sequencing depth, which is what DESeq2’s normalization handles. This is exactly what your project does: compare each gene in Flight vs. Ground.

**4. Sample 1 has 20M total reads and Sample 2 has 40M. A gene has 100 counts in Sample 1 and 200 in Sample 2. Is it more active in Sample 2?**
No. In both samples, the gene makes up 5 reads per million total reads. Sample 2 was simply sequenced twice as deeply. After normalization, they’re equal.

**5. Why does DESeq2 need raw counts, not already-normalized values?**
DESeq2’s statistical model is built specifically for count data. It uses the actual size of the counts to judge how reliable each measurement is, since a gene with 5 reads is far noisier than one with 5,000. Normalized values lose that information. DESeq2 also calculates its own normalization internally, so pre-normalizing would distort its results.

**6. Why doesn’t a regular t-test work well for RNA-seq with 5 samples per group?**
* RNA-seq counts aren’t normally distributed, which a t-test assumes.
* The variance grows with the mean, and there’s extra variation between biological replicates (overdispersion).
* With only 5 samples, the variance estimate for each gene is very unreliable, leading to false positives and missed genes.

DESeq2 fixes this with a negative binomial model and by borrowing information across all genes to stabilize each gene’s variance estimate.

**7. A gene has a log2 fold change of −2. What does that mean?**
The gene is one-quarter as active in Flight mice as in Ground mice (4x lower), since 2² = 4.

**8. WIf you test 20,000 genes at p < 0.05 and nothing is truly different, how many look significant by chance?**
About 1,000 (20,000 × 0.05). That’s why raw p-values alone are misleading here.

**9. A gene has padj = 0.03. What does that mean?**
If you call every gene with padj ≤ 0.03 significant, you’d expect about 3% of that list to be false positives. Technically, it’s the lowest false discovery rate at which this gene would make your list. It does not mean “there’s a 3% chance this specific gene is a false positive,” a common misinterpretation worth avoiding in interviews.

**10. On your PCA plot, samples split into “processed in March” vs. “April” rather than Flight vs. Ground. What’s the problem?**
That’s a batch effect: a technical difference (processing date) is causing more variation than the biology you care about. It can hide real flight effects or create fake ones. It’s especially bad if batch overlaps with condition, like all Flight samples processed in March, because then you can’t separate the two at all. If batch is balanced across groups, you can account for it by adding it to the DESeq2 design (~ batch + condition).

**11. What’s one advantage of GSEA over only analyzing your list of significant genes?**
GSEA uses all genes ranked, not just ones past an arbitrary cutoff. That means it can detect pathways where many genes each change a little in the same direction, even if few are individually significant. It also avoids results depending on whether you picked padj < 0.05 or < 0.1.

**12. Name three biological effects of spaceflight you might expect in your pathway results.**
* Mitochondrial dysfunction and oxidative stress
* Changes in lipid (fat) metabolism, especially in the liver
* Immune system and inflammation changes

**13. If a stress-related pathway is up in Flight mice, why can’t you conclude microgravity caused it?**
Flight mice experience more than microgravity: space radiation, launch and landing stress, higher CO₂ levels on the ISS, different handling, and sometimes different timing or methods of tissue collection. Any of these could cause stress responses. The experiment shows what changed with spaceflight as a whole, not microgravity specifically. Stating this clearly in your README’s limitations section shows good scientific judgment.



## Step 3: Data Exploration


**1. How many genes and samples?**
The counts table has [first number] genes and 27 samples (9 Flight, 9 Ground Control, 9 Vivarium). Each row is a gene, identified by an Ensembl ID (e.g., ENSMUSG00000000031), and each column is one mouse’s muscle sample.

**2. What do the sample names look like?**
You can distinguish between space flight, ground, and vivarium by the acronyms; FLT, GC, VIV.

**3. Whole numbers or decimals?**
The counts are stored as decimals (all 27 columns are float64). Most values are whole numbers displayed as .0, but a check across the whole table found that 5,785 genes contain at least one real decimal value (71,683 values in total). The first 5 rows happened to be whole numbers, which is why the preview looked misleading.

These decimals come from RSEM, which produces “expected counts”: when a read matches more than one gene, RSEM splits it between them based on probability (e.g., 0.7 to one gene and 0.3 to another).

The decimals are small relative to the counts (e.g., 87.51 vs. 115.0), so rounding has almost no effect. I’ll round all values to integers before running DESeq2, since it requires whole-number counts.

**4: How many genes have zero counts in every sample?**
20,976 genes (out of [total from Q1]) have zero counts in all 27 samples. This is expected: the table includes every annotated mouse gene, and many are only expressed in other tissues or are rarely expressed non-coding genes. These genes carry no information for comparing groups, so I’ll filter them out (along with very low-count genes) during QC. This also reduces the multiple-testing burden.

**5: Smallest and largest total reads per sample?**
Across all 27 samples, total reads range from 68.1 million (RR23_R-TA_GC_G6) to 116.4 million (RR23_R-TA_VIV_V7), a 1.71x difference. For the two groups I’m analyzing (Flight and Ground Control), the range is 68.1 million to 102.5 million (RR23_R-TA_FLT_F5), about 1.5x.

No sample stands out as unusually low, and the groups are well mixed in the sorted list. Flight and Ground Control have similar median depths (about 86.6 and 84.7 million), so neither group was sequenced systematically deeper. These modest depth differences are what normalization corrects for before comparing samples.

## Step 4: Quality control

- **ERCC removal:** Removed 92 ERCC spike-in controls before filtering. These are synthetic
  RNA controls, not mouse genes (see Decisions log).
- **Gene filtering:** Kept genes with ≥10 reads in ≥9 samples (the size of one group, so
  genes active in only one group survive). Kept 16,054 of 56,832 genes.
- **Library sizes:** 68–103 million reads per sample, similar between groups. No sample flagged.
- **PCA:** Flight and Ground separate completely along PC1 (41.8% of variance), indicating a
  strong group effect. This is down slightly from 43.4% with spike-ins included, so the
  separation is driven by real gene expression, not the spike-ins. Flight samples are more
  spread out on PC2 (13.9%), suggesting more variable responses among flight mice.
- **Correlation heatmap:** Hierarchical clustering splits samples perfectly into Flight and
  Ground. F9 is least similar to Ground samples; G1 and G8 are slightly less similar to other
  Ground samples. Results were unchanged after removing spike-ins.


## Step 5: Differential expression (DESeq2)

- **Design:** ~ condition, with Ground as the reference (positive log2FC = higher in Flight).
- **Results at adjusted p < 0.05:** 4,234 significant genes out of 16,054 tested (26%):
  2,179 up and 2,055 down in Flight.
- **ERCC check:** Before removing spike-ins, 4,293 genes were significant; nearly all 61
  spike-ins that passed filtering were among them. Top real genes were unchanged after removal.
- **Outliers:** 0 genes flagged by Cook's distance, supporting the decision to keep all samples.
- **Largest effects in the top 10:** ENSMUSG00000038216 (log2FC +4.0, ~16x up),
  ENSMUSG00000004939 (−3.9, ~15x down), ENSMUSG00000026822 (+3.8, ~14x up).
- **Comparison:** A 2026 reanalysis of this dataset reported ~4,400 significant genes; my
  result is consistent (difference likely due to my stricter filter and spike-in removal).
- **Interpretation note:** Significance doesn't imply a large change; many genes may have
  small fold changes. Some changes may also reflect landing and dissection timing.

## Step 6: Visualization and top genes

- **Effect sizes:** Of 4,234 significant genes, 306 changed at least 2-fold (160 up, 146 down
  in Flight). Most changes are modest; e.g., Clpx is the 2nd most significant gene but
  changed less than 2x.
- **Heatmap:** The top 30 genes cluster samples perfectly into Flight and Ground, with
  consistent patterns across mice. F8 and F9 show the strongest changes, matching the PCA.

### Top genes
- **Pnmt** (up ~16x): enzyme that makes adrenaline; increased by stress hormones.
- **Clpx** (up ~1.9x): mitochondrial protein-quality-control enzyme.
- **Nmrk2** (down ~15x): muscle-specific gene for making NAD+, needed for energy production.
- **Lcn2** (up ~14x): inflammation and stress-response protein.
- **Fkbp5** (up ~2.5x): classic marker of stress hormone (glucocorticoid) signaling.

### Emerging hypothesis
Several top genes (Fkbp5, Pnmt, Lcn2, Cdkn1a, Sesn1) respond to stress hormones, suggesting
a strong stress signature in Flight mice. Because mice were dissected the day after
splashdown and transport, part of this signal may reflect the return to Earth rather than
microgravity alone. Dbp (circadian clock gene, down ~4x) may reflect circadian disruption
or differences in dissection timing. To test in Step 7 with pathway analysis.




## Decisions log

### Study choice
**Chose OSD-576** (Rodent Research-23, right tibialis anterior muscle, ~38 days on the ISS).
- Reasons: strong differential expression signal in prior reanalysis, 9 mice per group,
  and muscle results are more consistent across missions than liver or retina.
- Ruled out: OSD-242 (ground controls from a different cohort), OSD-173 (13-day flight,
  collected hours after landing), OSD-137 (weaker, mixed liver signal), and
  OSD-326 / OSD-255 (couldn't verify design details).
- Limitations: mice were returned to Earth alive before dissection, and groups were
  dissected on different days (Flight Jan 14, Ground Control Jan 17, Vivarium Jan 20, 2021),
  so landing and timing effects may be mixed with spaceflight effects.

### Groups
**Comparing 9 Space Flight vs 9 Ground Control (habitat).**
- Excluding the 9 Vivarium controls, because they lived in standard cages rather than
  the same habitat hardware as the flight mice.

### Sex
**All male** (confirmed in the OSDR sample table).
- A 2026 PLOS ONE paper described these mice as female, but the OSDR metadata says male.
  Using OSDR as the source of truth.

### Counts file
**GLDS-570_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv**
- Using RSEM counts with ribosomal RNA genes removed, because leftover rRNA varies
  between samples and can distort normalization.
- Note: the file uses the ID GLDS-570 while the study is OSD-576. OSDR uses two ID
  systems that don't always match; confirmed the file is from the OSD-576 study page.
- Did not use normalized counts, because DESeq2 requires raw counts.

### Data handling
- **Rounding:** 5,785 genes contain decimal values (RSEM expected counts). Will round to
  integers before DESeq2, which requires whole-number counts.
- **Zero-count genes:** 20,976 genes are zero in all samples. Will filter these out
  (along with very low-count genes) during QC.
- **Sequencing depth:** Flight and Ground Control have similar median depths
  (~86.6M and ~84.7M reads), and no sample is unusually low. No samples excluded.

### Sample sheet
**Built data/samples.csv** from the OSDR sample table: 9 Flight and 9 Ground samples.
- Verified all sample names exist in the counts table and that group labels match
  the codes in the sample names (FLT, GC).

  ### Outliers
**Kept all 18 samples.** No sample crosses into the other group on PCA, and none has an
unusual library size. F9 is the most extreme Flight sample, but its difference goes in the
same direction as the spaceflight effect, suggesting biology rather than a technical problem.
Removing samples without a technical reason risks making results look cleaner than they are.

### Caveat on group separation
The PC1 separation reflects all differences between groups, not only spaceflight: Flight mice
returned to Earth before dissection and were dissected on a different day (Jan 14) than
Ground Controls (Jan 17).