open(GROUP, '<', "groups.txt");
open(RAW, '>', "Chi2RawOutput.txt");
open(SIGNIFICANT, '>', "Chi2SignificantGroup.txt");
open(ANNOTATED, '>', "Chi2SignificantGroup.annotated.txt");

my @groups;



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

foreach my $currentgroup (@groups){
        open (CHI2,'<',"$currentgroup.Chi2.result");
        my $chi2 = <CHI2>;
        open (H0, '<', "$currentgroup.H0.result.Ln0") or die "cannot open $currentgroup.H0.result";
        open (H1, '<', "$currentgroup.H1.result.Ln1") or die "cannot open $currentgroup.H1.result" ;
        my $h0=<H0>;
        my $h1=<H1>;
        my $lrt=2*($h1-$h0);
        print RAW ">$currentgroup\tlnL0\t$h0\tlnL1\t$h1\tLRT\t$LRT\tCHI2\t$chi2";
        if ($chi2 lt 0.05){
                print SIGNIFICANT "$currentgroup\n";
                open (ANNOTATION, "/work/zy/output/CDS_analysis/dNds_streamline/mcl_birds_groups.annotated") or die "cannot open annotation";
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

my $meetBEB = 0;
my $printHeader = 0;

foreach $currentgroup (@groups){
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
                                         open (ANNOTATION, "/work/zy/output/CDS_analysis/dNds_streamline/mcl_birds_groups.annotated") or die "cannot open annotation";
                while (<ANNOTATION>){
                        if(/$currentgroup/){
                                chomp $_;
                                print H1processed "$_\n";
                                last;
                        }
                }
                close ANNOTATION;
           $printHeader = 1;
        }
                                print H1processed "$lines\n";
                        }
                }
        }
        close H1RESULT;
        close H1processed;
}

foreach $currentgroup (@groups){
        open(processedHit, '>', "Chi2SignificantGroups.PassedBEB.txt");
        if ( -s "$currentgroup.H1.result.hit") {
                print processedHit "$currentgroup\n";
        }
}
