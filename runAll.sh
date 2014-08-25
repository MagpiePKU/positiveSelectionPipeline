### Zhang Yi
### 20140823
###

##. 1 Annotate MCL groups, prepare group and sequence input files

echo "start";
date;
time;
echo "MCL groups, magpie, groundtit, chicken, 3 branches";
echo "branch add on PICA";

cat ../mcl_birds_groups|grep PICA|grep PHUM|grep GGAL|sed -e 's/PICA|//g'|sed -e 's/PHUM|//g'|sed -e 's/GGAL|//g' >workingGroups
cp workingGroups Magpie.PHUM.GGAL.mcl.groups;

##sequence: chicken, crow, pica, groundjay. chicken: G, crow: C, pica: P, groundjay: H

rm mrna.fa;
touch mrna.fa;
cat /work/zy/genomicsequences/CDS/Gallus_gallus.Galgal4.75.cds.all.fa  >>mrna.fa
#cat /work/zy/genomicsequences/CDS/Corvus_brachyrhynchos.gene.cds.cleared.fa  >>mrna.fa
cat /work/zy/genomicsequences/Magpie/evm_nucl.curated.20140602version.fa  >>mrna.fa
cat /work/zy/genomicsequences/CDS/Ground_Tit.cds  >>mrna.fa

rm protein.fa;
touch protein.fa;
cat /work/zy/genomicsequences/pep/Gallus_gallus.Galgal4.75.pep.all.fa >>protein.fa
#cat /work/zy/genomicsequences/pep/Corvus_brachyrhynchos.gene.pep >>protein.fa
cat /work/zy/genomicsequences/Magpie/step1_prot.fa |sed -e 's/EVM_prediction_/evm.model./' >>protein.fa
cat /work/zy/genomicsequences/pep/Ground_Tit.pep >>protein.fa


##. 2 split groups and prepare AA, NA files
rm Group* -f;
perl splitgroups.pl workingGroups;
rm workingGroups -f;
ls Group_???? > groups.txt;
ls Group_????? >> groups.txt;
for a in `cat groups.txt`; do cat $a | awk '{FS=" "; OFS="\n"; for (i=1;i<=NF;i++) print $i}'>$a.column;done ;
for a in `cat groups.txt`; do faSomeRecords protein.fa $a.column $a.fasta; done;
for a in `cat groups.txt`;do cat $a.fasta |grep ">"|sed -e 's/>.\{0,1000\}transcript://'|sed -e 's/ .\{0,1000\}$//'|sed -e 's/\t.\{0,1000\}$//'|sed -e 's/>//' > $a.mrna.list; done; ## do the chicken protein-transript mapping. other species may not need this step
for a in `cat groups.txt`; do faSomeRecords mrna.fa $a.mrna.list $a.mrna.fasta; done;

##. 3 get MUSCLE alignment
rm *.muscle -f;
perl testMuscle.pl groups.txt;

##. 4 PAL2NAL
rm *.aa.alignment -f;
for a in `cat groups.txt`; do pal2nal.pl $a.muscle $a.mrna.fasta -nogap -codontable 1 -output paml >$a.aa.alignment; done;

##. 5 PhyML
ls *.aa.alignment > alignment.txt
rm *phyml -f;
perl runPhyml.pl;

##. 6 modify tree to add branch information
perl treelabel.pl; ### here add on pica branch only for now. 20140825

##. 7 PAML, process Ln0/Ln1, get LRT and Chi2, process to output final results;
perl testPAML2model.pl;
./Ln0Ln1.sh;
perl chi2test.pl;

echo "finished";
date;
time;
