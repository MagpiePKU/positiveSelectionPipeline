#############################################################################################################
# The program runs PAML in parallel way (currently set to 100, you can set higher in the denoted section)
# run the program as:
# perl ParallelPAML.pl your-regular-expression-for-labelling-trees (such as ENSP[^,\)]{0,1000} or ENSP\d+ is fine)
# Supposingly the directory is filled with previous runs of MUSCLE, stableMuscle, pal2nal, and phyml, containing:
#     groups.txt: a file containing all orthomcl groups
#     workinggroup.annotated: a file generated from the annotation.genelist.pl 
#     *.muscle: MUSCLE runs
#     *.aa.alignment: PAL2NAL runs
#     *.fasta and *.mrna.fasta (protein and NA)
#     *.tree.txt (phyml runs)
#     *.tree.txt.mod1 (modified phyml runs from ./correctTree.sh, also provided in this git. It removes all Bootstrap data)
# The program does:
# 1. labels tree branches as you denoted
# 2. parallelly runs PAML
# 3. performs Chi2 test, catching the LRT-positive groups, output the details
# 
# Note for an unfound (yet) bug it seems that the Chi2 session must be runed as a seperate pl file. Not fixed for now
#############################################################################################################


use Cwd;

open(GROUP, '<', "groups.txt");
my @groups;
my $currentgroup;
my $tree;
my $branch;

system "./correctTree.sh";

while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt.mod1" ){
                push (@groups, $_);
        }
}

close GROUP;


foreach $currentgroup (@groups) {
        open (TREE, '<', "$currentgroup.aa.alignment_phyml_tree.txt.mod1") or die "failed input PHYML tree";
        $tree = <TREE>;
        chomp $tree;
        print <$tree>;
        my @array =$tree =~ /($ARGV[0])/g;
                foreach $branch (@array){
                        print "branch $branch\t";
                        my $modbranch = "$branch #1";
                        print "modbranch $modbranch\n";
                        $tree =~ s/$branch/$modbranch/;
                }
        print "modified! $currentgroup\n";
        open (MODTREE, '>', "$currentgroup.aa.alignment_phyml_tree.txt.mod") or die "failed MODTREE handle";
        print MODTREE "$tree\n";
        close MODTREE;
        close TREE;
}


##############################################
### Running PAML
##############################################




open(GROUP, '<', "groups.txt");
my @groups ="";
my $currentgroup ="";



while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt.mod"){
                push (@groups, $_);
        }
}

foreach $currentgroup (@groups){

                my $token = 0;

                while ($token==0) {
                                $pid = `ps aux|grep codeml|wc -l`;
                                if ($pid <= 100) { ############ NOTE you can change your parallel threads here!
                                                $token = 1;
                                }
                                else {
                                        sleep (10);
                                }
                }

                my $suffix0 = "H0";
                my $suffix1 = "H1";

                my $H0dir = "$currentgroup.classical.H0";
                my $H1dir = "$currentgroup.classical.H1";
                system "mkdir $H0dir; mkdir $H1dir;";

                open(CTL,">$H1dir/codemlModel_H1.ctl");
        my $ctl_text = qq{
                seqfile = ../$currentgroup.aa.alignment
                treefile =  ../$currentgroup.aa.alignment_phyml_tree.txt.mod
                        outfile =  ../$currentgroup.H1.result
                noisy = 0  
                verbose = 0 
                runmode = 0
                seqtype = 1
                CodonFreq = 2 
                clock = 0
                model = 2
                NSsites = 2
                  icode = 0
                  fix_kappa = 0
                  kappa = 2
                fix_omega = 0
                omega = 1
                ncatG = 10
                getSE = 0
                RateAncestor = 0
                method = 0
                cleandata=1
                Small_Diff = .5e-6
                fix_blength = 1  
                };
        print CTL "$ctl_text";
        close CTL;
        
        
    open(CTL,">$H0dir/codemlModel_H0.ctl");
        my $ctl_text = qq{
                seqfile = ../$currentgroup.aa.alignment
                treefile =  ../$currentgroup.aa.alignment_phyml_tree.txt.mod
                outfile =  ../$currentgroup.H0.result
                noisy = 0  
                verbose = 0 
                runmode = 0
                seqtype = 1
                CodonFreq = 2 
                clock = 0
                        model = 2
                NSsites = 2
                icode = 0
                fix_kappa = 0
                kappa = 2
                fix_omega = 1
                omega = 1
                ncatG = 10
                getSE = 0
                RateAncestor = 0
                method = 0
                cleandata=1
                Small_Diff = .5e-6
                fix_blength = 1   
                };
        print CTL "$ctl_text";
        close CTL;
  
        my $original = cwd;
  
        chdir $H1dir;
        system "codeml codemlModel_H1.ctl &";
        chdir $original;
        chdir $H0dir;
        system "codeml codemlModel_H0.ctl &";
        chdir $original;

}

chdir $original;

close GROUP;

##############################################
### Processing Information, granting Chi2 and LRT 
##############################################


open(GROUP, '<', "groups.txt");
open(RAW, '>', "Chi2RawOutput.txt");
open(SIGNIFICANT, '>', "Chi2SignificantGroup.txt");
open(ANNOTATED, '>', "Chi2SignificantGroup.annotated.txt");

@groups = "";

system "./Ln0Ln1.sh;";

while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt.mod" && -s "$_.H0.result.Ln0" && -s "$_.H1.result.Ln1") {
                push (@groups, $_);
        }
}

close GROUP;

foreach my $currentgroup (@groups){
        open (H0, '<', "$currentgroup.H0.result.Ln0") or die "cannot open $currentgroup.H0.result";
        open (H1, '<', "$currentgroup.H1.result.Ln1") or die "cannot open $currentgroup.H1.result" ;
        my $h0=<H0>;
        my $h1=<H1>;
        my $lrt=2*($h1-$h0);
        open (LRT, '>', "$currentgroup.LRT.result") or die "cannot open LRT output";
        print LRT $lrt;
        close LRT;
        close H0;
        close H1;
}

my $wait = system "./LRT2Chi2.sh";

my @significantgroup;

foreach my $currentgroup (@groups){
        open (CHI2,'<',"$currentgroup.Chi2.result");
        my $chi2 = <CHI2>;
        open (H0, '<', "$currentgroup.H0.result.Ln0") or die "cannot open $currentgroup.PRANK.H0.result";
        open (H1, '<', "$currentgroup.H1.result.Ln1") or die "cannot open $currentgroup.PRANK.H1.result" ;
        my $h0=<H0>;
        my $h1=<H1>;
        my $lrt=2*($h1-$h0);
        print RAW ">$currentgroup\tlnL0\t$h0\tlnL1\t$h1\tLRT\t$LRT\tCHI2\t$chi2";
        if ($chi2 lt 0.05){
                print SIGNIFICANT "$currentgroup\n";
                push (@significantgroup, $currentgroup);
                open (ANNOTATION, '<', "./workinggroup.annotated") or die "cannot open annotation";
                while (<ANNOTATION>){
                        if(/$currentgroup/){
                                chomp $_;
                                print ANNOTATED $_,"\t",$chi2;
                                last;
                        }
                }
                close ANNOTATION;
        }
        close CHI2;
        close H0;
        close H1;
}

close RAW;
close SIGNIFICANT;


my $meetBEB = 0;
my $printHeader = 0;

foreach $currentgroup (@significantgroup){
        $meetBEB = 0;
        $printHeader = 0;
        open(H1RESULT,'<',"$currentgroup.H1.result");
        open(H1processed, '>', "$currentgroup.H1.result.hit");
        while (my $lines = <H1RESULT>){
                chomp $lines;
                if ($lines =~ /^Bayes Empirical Bayes/){
                        $meetBEB = 1;
                }
                if ($lines =~ /^The grid/) {
                        $meetBEB = 0;
                }
                if ($meetBEB == 1) {
                        if ($lines =~ /\*/) {
                                if ($printHeader == 0) {
                                         open (ANNOTATION, "./workinggroup.annotated") or die "cannot open annotation";
                while (<ANNOTATION>){
                        if(/$currentgroup/){
                                chomp $_;
                                print H1processed "$_\n";
                                last;
                        }
                }
                close ANNOTATION;
                open (TREE, '<',"$currentgroup.aa.alignment_phyml_tree.txt.mod");
                while (<TREE>) {
                        chomp $_;
                        print H1processed "$_\n";
                }
                close TREE;
                #open (TREE, '<',"$currentgroup.prank.best.anc.dnd");
                #while (<TREE>) {
                #        chomp $_;
                #        print H1processed "$_\n";
                #}
                #close TREE;
           $printHeader = 1;
        }
                                print H1processed "$lines\n";
                        }
                }
        }
        close H1RESULT;
        close H1processed;
}

open(processedHit, '>', "Chi2SignificantGroup.PassedBEB.txt");

foreach $currentgroup (@significantgroup){
        if ( -s "$currentgroup.H1.result.hit") {
                print processedHit "$currentgroup\n";
        }
}

system "rm -f Chi2SignificantGroup.PassedBEB.details.txt; cat *.H1.result.hit >> Chi2SignificantGroup.PassedBEB.details.txt;";

close processedHit;
