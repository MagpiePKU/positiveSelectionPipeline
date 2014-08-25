open(GROUP, '<', "groups.txt");
my @groups;
my $currentgroup;


while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt"){
                push (@groups, $_);
        }
}

foreach $currentgroup (@groups){
  open(CTL,">codemlModel_pairwise.ctl");
  my $ctl_text = qq{
    seqfile = $currentgroup.aa.alignment
    treefile =  $currentgroup.aa.alignment_phyml_tree.txt
    outfile =  $currentgroup.pairwise.result
 		noisy = 0  
    verbose = 0 
    runmode = -2
    seqtype = 1
    CodonFreq = 2
    clock = 0
    model = 0
    NSsites = 0
    icode = $icode
    fix_kappa = 0
    kappa = 2
    fix_omega = 0
    omega = 0.1
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 10
    getSE = 0
    RateAncestor = 0
    method = 0
    Small_Diff = .5e-6
	};
  print CTL "$ctl_text";
  close CTL;

  my $wait=system "codeml codemlModel_pairwise.ctl;";
}
