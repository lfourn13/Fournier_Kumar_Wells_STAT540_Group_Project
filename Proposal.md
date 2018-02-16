Proposal
================
Louis-Alexandre Fournier
2/13/2018

**Research Proposal STAT 540**
==============================

Team: "Splice and dice" Members: Arun Kumar @arunk95, James Wells @jameswells96 and Louis-Alexandre Fournier @lfourn13

**Motivation and background work**
----------------------------------

Our final project will involve RNA-Seq data from the Peter Stirling lab (Focus on genomic instability and cancer). Yeast wild-type data will be compared to two different mutants in the yeast homolog of SF3B1 (splicing factor 3b subunit 1) called HsH155. SF3B1 is a component of several splicesomes and is important translational modification of proteins to ensure proper functionality. SF3B1 mutations have been detected in several cancers including breast, uveal melanoma and myeloid leukemias. Two loss of function mutations commonly seen in the HEAT repeat domains of cancer patients: K335N and P369E. These mutations can be modelled in a yeast system in HSH155 to understand the consequence of loss of splicesomal function. Specifically, we want to know how these splice mutants are implicated in gene expression and intron retention rate and how this may induce a cancerous phenotype. SF3B1 is the most mutated gene in myeloid maligancies and is indicative of better overal survival compared to wildtype counterparts. However, SF3B1 is a poor prognostic indicator in progesterone receptor-negative breast cancer. This comparison highlights the gap in our understanding of SF3B1 mutation.

**Division of labor**
---------------------

Yeast wild-type data will be compared to two different splice mutants, all previously collected in triplicate. Our project aims to uncover specific gene expression profile variants between the control (WT) and splice mutants that can be associated with the diseased phenotypes. In order to accomplish this, our project will be divided into three parts, each associated with one member of our team:

1.  *Differential gene expression* = person responsible to write the code that will identify the top genes from the RNA-seq data that are differentially expressed between the control and experimental groups. This will allow to correlate different gene expression patterns to observed phenotypes.

2.  *Quality Control* = person responsible for inspecting the dataset and analyzing the quality control data generated from the RNA-seq experiment. This person will also oversee the the progression of the project with the other team members.

3.  *Intron retention* = person responsible to identify genes from the RNA-seq data that show significantly high or low retention rates as a means to correlate intron retention to observed phenotypes. Investigating intron retention is important since splice variants can have different functions leading to distinct phenotypes.

``` r
member <- matrix(c("Microbiology and Immunology", "Biotechnology","Biochemistry", "MSc IOP", "MSc Medical Genetics", "MSc IOP","Stirling Lab", "Stirling Lab", "Stirling Lab", "Differential gene expression", "Quality control", "Intron retention"), ncol=3, byrow=TRUE)
colnames(member) <- c("Louis-Alexandre Fournier", "Arun Kumar", "James Wells")
row.names(member) <- c("background", "degree", "affiliation", "job assignment")

member <- as.table(member)
member
```

    ##                Louis-Alexandre Fournier     Arun Kumar          
    ## background     Microbiology and Immunology  Biotechnology       
    ## degree         MSc IOP                      MSc Medical Genetics
    ## affiliation    Stirling Lab                 Stirling Lab        
    ## job assignment Differential gene expression Quality control     
    ##                James Wells     
    ## background     Biochemistry    
    ## degree         MSc IOP         
    ## affiliation    Stirling Lab    
    ## job assignment Intron retention

**Dataset**
-----------

RNA was extracted from Saccharomyces cerevisiae using the RiboPure RNA purification kit for Yeast. Steps were performed as per the manufacturer's instructions along with an additional DNase step for maximum purity. 9 samples in total - one wild type, one point mutant K335N, and another point mutant P369E (in triplicates) were sent to Genome Quebec Innovation Centre for sequencing. The library was prepared using Illumina TruSeq rRNA depleted stranded library prep. Based on a literature survey, rRNA depletion was chosen over PolyA since rRNA depletion protocols appear to measure immature transcripts (pre-mRNAs) and therefore provide more information on splicing patterns and possible splice junctions. Furthermore, since the main area of interest were rare splice variants, the strongest library prep method was chosen. Sequencing was performed by HiSeq 4000 Paired End 100bp sequencing lane yielding &gt; 300M reads, ~25M reads per library.

**Aims and methodology**
------------------------

intron retention: using the SpliceR package <http://www.bioconductor.org/packages/2.13/bioc/html/spliceR.html> IRFinder <https://omictools.com/irfinder-tool>

Please let me know if you require anymore information about our proposed project. Group members: Me, @arunk95, and @lfourn13