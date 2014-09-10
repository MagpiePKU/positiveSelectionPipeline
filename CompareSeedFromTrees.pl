###################################################################################################################
# The program compares two "seed" species between different runs of positive selection scans
# What you need is a config file, ncbi blast+, orthomcl, any gene groups of your species (like generated from a 
# previous orthomcl run with A/B/C species, and D/E/F species, and now you want to combine A/B/C/D/E/F), peptide
# sequences of these species. Here we use "species n chi2" because these groups are generated from 
# positive selection scan and are chi2-test-positive gene groups (under positive selection). 
# We use the program for comparing two chi2-test-positive set of groups, seeing whether one family of gene undergoes 
# positive selection in both trees or not. 
#
# run program as:
# perl CompareSeedFromTrees.pl your-config-file your-sql-password & 
# (I suppose you use root for sql! Otherwise you may change the -u in )
#
# provide config file as:
# ==========================
# species 1     chi2 XXX
# species 2     chi2 YYY
# species 1 fasta       XXX.fa
# species 2 fasta       YYY.fa
# species 1 REGX        XXX-regx
# species 2 REGX        YYY-regx
# ==========================
# species 1 is the query, species 2 is the subject!
# The program will run:
# 1.Fetch FA files for the hit list you provided as XXX, YYY (Chi2 positive hits);
# 2.Make blastdb and BLASTP for XXX against YYY;
# 3.run orthomcl for XXX against YYY;
# 4.combine Chi2positivehits from XXX and YYY as you provided and compare to the new (temp) orthomcl groups, find joint groups. (NOT DONE FOR NOW)
###################################################################################################################

use Cwd;

open (CONFIG, '<', $ARGV[0]);

my $SPECIES1REGX;
my $SPECIES2REGX;
my $inputgroup1 =0;
my $inputgroup2 =0;
my $species1regx =0;
my $species2regx =0;
my $species1fasta =0;
my $species2fasta =0;
my %totalGeneTable;

while (<CONFIG>){
                chomp $_;
                @fields = split "\t", $_;
                if (/species 1 chi2/) {
                                open (INPUTGROUP1, '<', $fields[1]); $inputgroup1 = $fields[1]; }
                if (/species 2 chi2/) {
                                open (INPUTGROUP2, '<', $fields[1]); $inputgroup2 = $fields[1];}
                if (/species 1 REGX/) {
                                print "$fields[1]\n";
                                $SPECIES1REGX = $fields[1]; $species1regx = 1;}
                if (/species 2 REGX/) {
                                print "$fields[1]\n";
                                $SPECIES2REGX = $fields[1]; $species2regx = 1;}
                if (/species 1 fasta/) {
                                print "$fields[1]\n";
                                $SPECIES1fasta = $fields[1]; $species1fasta = 1;}
                if (/species 2 fasta/) {
                                print "$fields[1]\n";
                                $SPECIES2fasta = $fields[1]; $species2fasta = 1;}
}


open (QUERY1, '>', "$inputgroup1.query.list");
open (QUERY2, '>',"$inputgroup2.query.list");

while (<INPUTGROUP1>){
                @catched = $_ =~ m/($SPECIES1REGX)/g;
                foreach $item (@catched) {
                                if (defined ($totalGeneTable{$item})) {next;}
                                else { print QUERY1 "$item\n"; 
                                $totalGeneTable{$item} = 1; }
                }
}

while (<INPUTGROUP2>){
                @catched = $_ =~ m/($SPECIES2REGX)/g;
                foreach $item (@catched) {
                                if (defined ($totalGeneTable{$item})) {next;}
                                else { print QUERY2 "$item\n"; 
                                $totalGeneTable{$item} = 1; }
                }
}


$orthomclconfig = qq{dbVendor=mysql 
dbConnectString=dbi:mysql:database=orthomcl_temp;mysql_socket=/work/mysqldata/mysql.sock
dbLogin=orthomcl
dbPassword=orthomcl
similarSequencesTable=similarsequences
orthologTable=ortholog
inParalogTable=inparalog
coOrthologTable=coortholog
interTaxonMatchView=intertaxonmatch
percentMatchCutoff=50
evalueExponentCutoff=-6
oracleIndexTblSpc=NONE
};

$sqlfile = qq {
drop database orthomcl_temp;
create database orthomcl_temp;
grant all privileges on orthomcl_temp.* to "orthomcl"@"localhost";
};

open (ORTHOMCLCONFIG, ">orthomcl.temp.config");
open (SQLCMD, ">temp.sql");

print ORTHOMCLCONFIG $orthomclconfig;
print SQLCMD $sqlfile;

close ORTHOMCLCONFIG;
close SQLCMD;



$blastcmd = qq{ 
        cat $inputgroup1.query.list >> query.list; 
        cat $inputgroup2.query.list >> query.list; 
        faSomeRecords $SPECIES1fasta $inputgroup1.query.list $inputgroup1.fasta;
        faSomeRecords $SPECIES2fasta $inputgroup2.query.list $inputgroup2.fasta;
        orthomclAdjustFasta SPE1 $inputgroup1.fasta 1;
        orthomclAdjustFasta SPE2 $inputgroup2.fasta 1;
        rm -R -f temp.fasta;
        mkdir temp.fasta;
        mv SPE1.fasta temp.fasta;
        mv SPE2.fasta temp.fasta;
        orthomclFilterFasta ./temp.fasta 100 10;
        makeblastdb -in goodProteins.fasta -input_type fasta -dbtype prot -title temp.db -out temp.db;
        blastp -db temp.db -query goodProteins.fasta -evalue 1e-6 -max_target_seqs 50 -num_threads 96 -outfmt 6 -out temp.blastp.out;
};

$orthomclcmd = qq{
        mysql --socket=/work/mysqldata/mysql.sock -uroot --password="$ARGV[1]" < temp.sql; rm temp.sql -f;
        orthomclInstallSchema orthomcl.temp.config; 
        orthomclBlastParser temp.blastp.out ./temp.fasta >> temp.orthomcl.similarsequences.txt;
        orthomclLoadBlast orthomcl.temp.config temp.orthomcl.similarsequences.txt;
        orthomclPairs orthomcl.temp.config orthomcl.temp.pairs.log cleanup=no;
        orthomclDumpPairsFiles orthomcl.temp.config;
        mcl mclInput --abc -I 1.5 -o mclOutput;
        orthomclMclToGroups TEMP 00000 <mclOutput >temp.orthomcl.groups;
};

system $blastcmd;
system $orthomclcmd;
