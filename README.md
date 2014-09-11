positiveSelectionPipeline
=========================
updated 20140911


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
  f. PAL2NAL (I use v12)
  g. PhyML (I use the 20120412 version)
  h. PAML (I use 4.8)
  i. Cufflinks with the gffread program

Special notes on orthoMCL: we now recommend a "correct blast" and a "correct installation of MySQL". 

By "correct blast", we mean blasting all proteins from all species at once. Splitting the blast to different machines, such as blasting proteins from a single species against all protein database in one single node, might have cause an unknown problem. A large all-vs-all blast is very expensive, and we have developed a tool for you that you could select a few "seed" species, subdivide your species to groups that contain any of these seeds, perform positive selection test on these subgroups, and then combine and converge those "selected gene families" from species groups. 

By "correct installation of MySQL", we now found that MySQL might be quite problematic. So titrate the MySQL config file before running large sets of data. As a standard, we performed 8-mammalian-species orthomcl on a CentOS machine (E5-4680) with single CPU within 1.5 hours. Anything longer than this is abnormal and should be stopped and check. As a superstitution, we used:

myisam_sort_buffer_size=700G
myisam_max_sort_file_size=4000G

and nothing else. This worked. Anything other than this, like changing the key buffer, join buffer, sort buffer etc, does not work in our hand. We manually set MySQL to use MySIAM but not InnoDB. 

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
    viii. run PAML on model H0/H1, both are model=2, NSsites=2, runMode=0. The difference is H0 is fixed omega=1 and H1 is free omega  (NOTE now we provide a script that runs PAML in parallel)
    ix.   use shell script and chi2test.pl to get LRT results for each gene group. Significant hits are recorded under Chi2Significant* outputs. Narrowed down LRT significant hit list to those with positively selected sites under Bayes Empirical Bayes are also processed. 
  
  f.  run a perl script which goes through your 'seed' species, combine the positive selection results, and converge them to find any converging evolution signature. 
    
3. Outputs:
  a. Chi2RawOutput.txt  //Chi2 result for every gene group
  b. Chi2SignificantGroup.annotated.txt  //Significant Chi2 (P<0.05) groups, with annotation. 
  c. Chi2SignificantGroups.PassedBEB.txt  //Significant Chi2 (P<0.05) groups, with significant BEB hit
  d. Chi2SignificantGroup.txt  //Significant Chi2 (P<0.05) groups, group name only list
  e. Chi2SignificantResultDetails.txt //Significant Chi2 (P<0.05) groups, detailed codeml output
  f. Group_xxxx.H0.result and Group_xxxx.H1.result  //codeml outputs 
  g. Group_xxxx.H1.result.hit //codeml significant BEB hit only output of Significant Chi2 (P<0.05) groups
  h. Outputs from the CompareSeedFromTrees.pl (labelled in that script). 

4. References:
  a. PAML: http://abacus.gene.ucl.ac.uk/software/paml.html
  b. Satsuma: https://www.broadinstitute.org/science/programs/genome-biology/spines
  c. kraken: https://github.com/nedaz/kraken
  d. orthoMCL: http://orthomcl.org/
  e. UCSC utilities: http://hgdownload.cse.ucsc.edu/admin/exe/
  f. MUSCLE: http://www.drive5.com/muscle/
  g. PAL2NAL: http://www.bork.embl.de/pal2nal/
  h. PhyML: http://code.google.com/p/phyml/
  i. cufflinks: http://cufflinks.cbcb.umd.edu/
  j. model H0/H1: Clark, A. G. et al. 2003. Inferring non-neutral evolution from human-chimp-mouse orthologous gene trios. Science 302: 1960-1961.
