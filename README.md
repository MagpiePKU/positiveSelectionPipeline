positiveSelectionPipeline
=========================

Random perl and shell scripts plus notes for whole genome CDS positive selection detection

Zhang Yi, 20140825, synapse@pku.edu.cn

Please simply cite this github if you use the scripts. 

Note I dont know why the github readme looks bad on web. In UltraEdit the lines and spaces seems fine...

1. Prerequistiques:
  The programs run under CentOS 6.5 with perl 5.10. Requires the following softwares to be installed to run properly:
  a. orthoMCL
  b. SatsumaSynteny
  c. Kraken
  d. UCSC faSomeRecords
  e. MUSCLE (I use muscle3.8.31_i86linux64)
  e. PAL2NAL (I use v12)
  f. PhyML (I use the 20120412 version)
  g. PAML (I use 4.8)
  h. Cufflinks with the gffread program


2. Working pipeline:
  a. protein grouped by orthoMCL (times undefined. For what we know the orthoMCL system on CentOS is very problematic. Good luck.)
  b. genomes mapped by SatsumaSynteny (~1 day for 1.2G vs 1.2G genome)
  c. Cross-check synteny with Kraken using the genome pairwise maps produced in b. (~8 hrs for 1.2G genome)
      if you have only GFF3 annotation file, change it with gffread into GTF. 
  (the scripts doing a,b,c are currently not out yet)

  d. Annotate the groups with annotate.genelist.pl and Biomart-downloaded gene list. 
  e. use ./runAll.sh [modify ./runAll.sh, and treelabel.pl before using them!]:
    i.    splits orthoMCL-and-synteny passed gene groups from the output file
    ii.   creates mrna and protein fasta bundle file of your chosen species
    iii.  fetches mrna and protein fasta files for each individual gene group 
    iv.   do MUSCLE alignment on the protein fasta files
    v.    use PAL2NAL to change protein alignment into codon-matched alignment
    vi.   PhyML builds tree with 100 bootstrap loop
    vii.  label tree branch using the option you choose (HERE always change treelabel.pl before use)
    viii. run PAML on model H0/H1, both are model=2, NSsites=2, runMode=0. The difference is H0 is fixed omega=1 and H1 is free omega  
    ix.   use shell script and chi2test.pl to get LRT results for each gene group. Significant hits are recorded under Chi2Significant* outputs. Narrowed down LRT significant hit list to those with positively selected sites under Bayes Empirical Bayes are also processed. 
    
3. Outputs:
  a. Chi2RawOutput.txt  //Chi2 result for every gene group
  b. Chi2SignificantGroup.annotated.txt  //Significant Chi2 (P<0.05) groups, with annotation. 
  c. Chi2SignificantGroups.PassedBEB.txt  //Significant Chi2 (P<0.05) groups, with significant BEB hit
  d. Chi2SignificantGroup.txt  //Significant Chi2 (P<0.05) groups, group name only list
  e. Chi2SignificantResultDetails.txt //Significant Chi2 (P<0.05) groups, detailed codeml output
  f. Group_xxxx.H0.result and Group_xxxx.H1.result  //codeml outputs 
  g. Group_xxxx.H1.result.hit //codeml significant BEB hit only output of Significant Chi2 (P<0.05) groups
