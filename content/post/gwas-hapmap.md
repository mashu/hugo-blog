
# Table of Contents

1.  [Background](#org4b2cdd0)
2.  [Data](#orgf1a3d6e)

&#x2014;
title: "Genome wide association studies"
date: 2018-05-10T19:08:36+02:00
draft: false
&#x2014;


<a id="org4b2cdd0"></a>

# Background

Genome Wide Association Studies (GWAS) make use of genetic data that are a collection of genetic variants such as single nucleotide polymorphism (SNP) across large population. Such data can be used to find which genetic variants are associated with certain phenotype traits.

Analysis requires two sources of information, genetic variant data for instance SNP for multiple individuals and data about their phenotype. The goal is to find which regions in the genome measured by means of SNPs co-vary with certain trait.


<a id="orgf1a3d6e"></a>

# Data

I will be using HapMap data <sup id="2a4905183ee50ac8fc1b1dbb3f5f38a2"><a href="#HapMap2010" title="Consortium HapMap3, Integrating common and rare genetic variation in  diverse human populations, {Nature}, v(7311), 52&#8211;58 (2010).">HapMap2010</a></sup> which can be downloaded from [link](ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/).
The **map** file contains SNP names, and **ped** file is for genotype information. 

    hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.map
    hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.ped

I will use R package GenABEL <sup id="e2cc37ccc1dfabf7f6a1dff52a849675"><a href="#Aulchenko2007" title="Aulchenko, Ripke, Isaacs, \&amp; van Duijn, GenABEL: an R library for genome-wide association  analysis, {Bioinformatics}, v(10), 1294&#8211;1296 (2007).">Aulchenko2007</a></sup>.
First convert data to GenABEL format

    library(GenABEL)
    convert.snp.ped(
        pedfile = '/home/mateusz/Downloads/GWAS-Phaselll-01-2009/hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.ped',
        mapfile = '/home/mateusz/Downloads/GWAS-Phaselll-01-2009/hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.map',
        outfile = '/home/mateusz/Downloads/GWAS-Phaselll-01-2009/hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.out')

Load phenotype data, here sex and pedigree 

    library(GenABEL)
    # Load tab-delimited file with phenotype data and fix factor levels
    pheno = read.table('/home/mateusz/Downloads/GWAS-Phaselll-01-2009/relationships_w_pops_121708.txt',
                       header = TRUE,stringsAsFactors = FALSE)
    pheno$sex = pheno$sex -1
    write.table(pheno,
                '/home/mateusz/Downloads/GWAS-Phaselll-01-2009/pheno.txt')
    
    # Load combined data
    data = load.gwaa.data(
        phenofile= "/home/mateusz/Downloads/GWAS-Phaselll-01-2009/pheno.txt",
        genofile = "/home/mateusz/Downloads/GWAS-Phaselll-01-2009/hapmap3_pop/hapmap3_r2_b36_fwd.CEU.qc.poly.out",
        id='IID')
    
    # Quality check
    qc = check.marker(data)
    
    # Select samples and SNPs that are ok by default filtering
    fdata = data[qc$idok,qc$snpok]
    mdata = phdata(fdata)
    
    nids(fdata)  # Number of individuals
    nsnps(fdata) # Number of SNPs
    
    # Fit logistic model, makes no sense, but just test the method
    result = scan.glm("sex~CRSNP",
                   data=fdata,
                   family=binomial(),
                   snps=(1:1000)) # Otherwise too long for this test
    
    result$map[result$P1df < 0.05] # P-values testing association to select SNPs


# Bibliography
<a id="HapMap2010"></a>[HapMap2010] Consortium HapMap3, Integrating common and rare genetic variation in  diverse human populations, <i>{Nature}</i>, <b>467(7311)</b>, 52–58 (2010). <a href="http://dx.doi.org/10.1038/nature09298">link</a>. <a href="http://dx.doi.org/10.1038/nature09298">doi</a>. [↩](#2a4905183ee50ac8fc1b1dbb3f5f38a2)

<a id="Aulchenko2007"></a>[Aulchenko2007] Aulchenko, Ripke, Isaacs, \& van Duijn, GenABEL: an R library for genome-wide association  analysis, <i>{Bioinformatics}</i>, <b>23(10)</b>, 1294–1296 (2007). <a href="http://dx.doi.org/10.1093/bioinformatics/btm108">link</a>. <a href="http://dx.doi.org/10.1093/bioinformatics/btm108">doi</a>. [↩](#e2cc37ccc1dfabf7f6a1dff52a849675)

