Proposal
================
Louis-Alexandre Fournier
2/13/2018

**Research Proposal STAT 540**
==============================

Team: "Splice and dice" Members: Arun Kumar @arunk95, James Wells @jameswells96 and Louis-Alexandre Fournier @lfourn13

**Motivation and background work**
----------------------------------

Our final project will involve RNA-Seq data from the Peter Stirling lab (Focus on genomic instability and cancer). Yeast wild-type data will be compared to two different mutants in the yeast homolog of SF3B1 (splicing factor 3b subunit 1) called HSH155. SF3B1 is a component of several splicesomes and is important for translational modification of proteins that ensure proper functionality. SF3B1 mutations have been detected in several cancers including breast, uveal melanoma and myeloid leukemias. Two SF3B1 loss of function mutations that have been commonly seen in cancer patients are K666T/R and K700E/R. These mutations have been modelled in a yeast system in HSH155 (K335N and P369E respectively) to understand the consequence of loss of splicesomal function. Specifically, we want to know how these splice mutants are implicated in gene expression and intron retention rate and how this may induce a cancerous phenotype. We expect loss of HSH155 function in yeast will result in higher intron retention and/or altered gene expression in genes associated with growth, proliferation or DNA repair, all of which are pathways commonly altered in cancer. SF3B1 is the most mutated gene in myeloid malignancies and is indicative of better overall survival compared to wildtype counterparts. However, SF3B1 is a poor prognostic indicator in progesterone receptor-negative breast cancer. This comparison highlights the gap in our understanding of SF3B1 mutation. Ultimately the differences in gene expression and intron retention identified from the RNAseq data will provide insight on the impact of abberrations in splicing factors, such as HSH155, in the context of cancer. Due to the homologous nature of HSH155 and SF3B1, these findings may have the potential to be translated to human systems.

**Division of labor**
---------------------

Yeast wild-type data will be compared to two different splice mutants, all previously collected in triplicate. Our project aims to uncover specific gene expression profile variants between the control (WT) and splice mutants that can be associated with the diseased phenotypes. In order to accomplish this, our project will be divided into three parts, each associated with one member of our team:

1.  *Quality Control* = person responsible for inspecting the dataset and analyzing the quality control data generated from the RNA-seq experiment. This person will also import the data, create the design matrix and oversee the progression of the project with the other team members.

2.  *Differential gene expression* = person responsible for writing the code that will identify the top genes from the RNA-seq data that are differentially expressed between the control and experimental groups. This will allow us to correlate different gene expression patterns to observed phenotypes.

3.  *Intron retention* = person responsible for identifing genes from the RNA-seq data that show significantly high or low retention rates as a means to correlate intron retention to observed phenotypes. Investigating intron retention is important since splice variants can have different functions leading to distinct phenotypes.

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

Questions/Aims:

1.  What genes are differentially expressed between the control and experimental groups. Knowing that the human homolog of yeast HSH155, SF3B1, binds DNA repair proteins, we suspect DNA repair coding genes to be effected in the mutant groups. Initially we will convert the RNA-seq data in BAM file format to Fastq files and format it in R. The data will be filtered to remove lowly expressed genes and visualize the data using hierarchical clustering with heatmaps. The data will be normalized to eliminate composition bias. Next we will use [limma with voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29) on the normalized data to assess differential gene expression. To achieve this, a design matrix will be formed which will then be subjected to voom transformation. Depending on the nuimber of differentially expressed gene hits, we will perform gene set testing using the goana function in limma to identify specific pathways effected by the point mutations ine experimental samples.

2.  Is there a difference in splicing efficiency between the control and experimental groups. This can be determined by comparing intron retention rates. Previous groups have had success using the [SpliceR package](http://www.bioconductor.org/packages/2.13/bioc/html/spliceR.html) a bioconducter R package that can determine alternative splicing from RNA-seq data. Additionally, [IRFinder](https://omictools.com/irfinder-tool) is another option which has been shown to have identify intron retention which [higher accuracy](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4) than other comparable software tools such as MISO and DEXseq.

Please let us know if you require anymore information about our proposed project. Group members: @Jameswells96, @arunk95, and @lfourn13
