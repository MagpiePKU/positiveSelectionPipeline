open(GROUP, '<', "groups.txt");
my @groups;
my $currentgroup;


while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt.mod"){
                push (@groups, $_);
        }
}

foreach $currentgroup (@groups){
  open(CTL,">codemlModel_H1.ctl");
  my $ctl_text = qq{
    seqfile = $currentgroup.aa.alignment
    treefile =  $currentgroup.aa.alignment_phyml_tree.txt.mod
    outfile =  $currentgroup.H1.result
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
  open(CTL,">codemlModel_H0.ctl");
  my $ctl_text = qq{
    seqfile = $currentgroup.aa.alignment
    treefile =  $currentgroup.aa.alignment_phyml_tree.txt.mod
    outfile =  $currentgroup.H0.result
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
  my $wait=system "codeml codemlModel_H0.ctl ; codeml codemlModel_H1.ctl;";
}
