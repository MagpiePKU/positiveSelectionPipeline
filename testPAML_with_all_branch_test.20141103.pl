use Cwd;

##############################################
### Readin config file
##############################################

open (CONFIG, '<', $ARGV[0]);

%SPECIESANDREGX;

while (<CONFIG>) {
                chomp $_;
                @array = split "\t", $_;
                $SPECIESANDREGX{$array[0]} = $array[1];
}

print join ("\n", %SPECIESANDREGX);

##############################################
### Making tree 
##############################################

open(GROUP, '<', $ARGV[1]);
my @groups;
my $currentgroup;
my $tree;
my $branch;
my %treehash;

while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt.mod1"){
                push (@groups, $_);
                open (TREE, '<', "$_.aa.alignment_phyml_tree.txt.mod1");
                $treehash{$_}= <TREE>;
                chomp $treehash{$_};
                close TREE;
        }
}

close GROUP;


foreach $currentgroup (@groups) {
        $tree = $treehash{$currentgroup};
    foreach $species (sort keys %SPECIESANDREGX ){
                if (-e "$currentgroup.aa.alignment_phyml_tree.$species.txt.mod"){next;}
                my $token = 0;
                while ($token==0) {
                                $pid = `ps aux|grep codeml|wc -l`;
                                if ($pid < 90) {
                                                $token = 1;
                                }
                }       
        my $temptree = $tree;       
        my @array =$temptree =~ /($SPECIESANDREGX{$species})/g;
                foreach $branch (@array){
                        print "branch $branch\t";
                        my $modbranch = "$branch #1";
                        print "modbranch $modbranch\n";
                        $temptree =~ s/$branch/$modbranch/;
                }
        print "modified! $currentgroup\t$species\n$temptree\n";
        open (MODTREE, '>', "$currentgroup.aa.alignment_phyml_tree.$species.txt.mod") or die "failed MODTREE handle";
        print MODTREE "$temptree\n";
        close MODTREE;


#############################################


##############################################
### Running PAML
##############################################
                my $H0dir = "$currentgroup.$species.H0";
                my $H1dir = "$currentgroup.$species.H1";
                system "mkdir $H0dir& mkdir $H1dir";

                open(CTL1,">$H1dir/codemlModel.$currentgroup.$species.H1.ctl");
        my $ctl_text = qq{
                seqfile = ../$currentgroup.aa.alignment
                treefile =  ../$currentgroup.aa.alignment_phyml_tree.$species.txt.mod
                        outfile =  ../$currentgroup.$species.H1.result
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
        print CTL1 "$ctl_text";
        close CTL1; 


    open(CTL2,">$H0dir/codemlModel.$currentgroup.$species.H0.ctl");
        my $ctl_text = qq{
                seqfile = ../$currentgroup.aa.alignment
                treefile =  ../$currentgroup.aa.alignment_phyml_tree.$species.txt.mod
                outfile =  ../$currentgroup.$species.H0.result
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
        print CTL2 "$ctl_text";
        close CTL2;

        my $original = cwd;

        chdir $H1dir;
        system "codeml codemlModel.$currentgroup.$species.H1.ctl &";
        chdir $original; 
        chdir $H0dir;
        system "codeml codemlModel.$currentgroup.$species.H0.ctl &";
        chdir $original;

        }
}

#######################################
## Running LRT (d.f. = 1)                                                       #
#######################################

open (OUTPUTChi2Significant, '>', "Chi2.speciesdivided.Significant.output");
open (OUTPUTChi2Raw, '>', "Chi2.speciesdivided.Raw.output");

sleep (300);

foreach $currentgroup (@groups) {   
                foreach $species (sort keys %SPECIESANDREGX ){
                                        my $lnLH0 = 0;
                                        my $lnLH1 = 0;
                                        my $LRT = 0;
                                        open (H0RESULT, ,'<', "$currentgroup.$species.H0.result");
                                        open (H1RESULT, ,'<', "$currentgroup.$species.H1.result"); 
                                        while (<H0RESULT>){
                                                        if ($_ =~ /lnL.+-(\d+.\d+)/){
                                                                        $lnLH0 = -1*$1;
                                                                        last;
                                                        }
                                        }
                                        while (<H1RESULT>){
                                                        if ($_ =~ /lnL.+-(\d+.\d+)/){
                                                                        $lnLH1 = -1*$1;
                                                                        last;
                                                        }
                                        }
                                        $LRT= 2*($lnLH1-$lnLH0);
                                        if ($LRT>6.9){
                                                        print OUTPUTChi2Significant "$currentgroup\t$species\tlnL-H0\t$lnLH0\tlnL-H1\t$lnLH1\tLRT\t$LRT\n";
                                        }
                                        print OUTPUTChi2Raw "$currentgroup\t$species\tlnL-H0\t$lnLH0\tlnL-H1\t$lnLH1\tLRT\t$LRT\n";
                }
}

close OUTPUTChi2Significant;
close OUTPUTChi2Raw;

my $meetBEB = 0;
my $printHeader = 0;

foreach $currentgroup (@groups){
foreach $species (sort keys %SPECIESANDREGX ){
        $meetBEB = 0;
        $printHeader = 0;
        open(H1RESULT,'<',"$currentgroup.$species.H1.result");
        open(H1processed, '>', "$currentgroup.$species.H1.result.hit") or die "cannot open hit";
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
                                if ($printHeader == 0 ){
                                print H1processed "$significantgroups{$currentgroup}\n";
           $printHeader = 1;
        }
                                print H1processed "$lines\n";
                        }
                }
        }
        close H1RESULT;
        close H1processed;
}
}

open(processedHit, '>', "Chi2.$species.Significant.output.PassedBEB.annotated.txt") or die "cannot open passbeb output";

foreach $currentgroup (@groups){
foreach $species (sort keys %SPECIESANDREGX ){
        if (-s "$currentgroup.$species.H1.result.hit"){
                open (HIT, "<", "$currentgroup.$species.H1.result.hit");
                while (<HIT>){
                        print processedHit $_;}
                close HIT;
                print processedHit "$currentgroup\n";}
}
}
