open(GROUP, '<', "groups.txt");
my @groups;
my $currentgroup;
my $tree;
my $branch;


while (<GROUP>){
        chomp $_;
        if (-s "$_.aa.alignment_phyml_tree.txt" ){
                push (@groups, $_);
        }
}

close GROUP;


foreach $currentgroup (@groups) {
        open (TREE, '<', "$currentgroup.aa.alignment_phyml_tree.txt") or die "failed input PHYML tree";
        $tree = <TREE>;
        chomp $tree;
        print <$tree>;
        $tree =~ m/(evm.model.scaffold[^,\)]{0,1000})/; ## NOTE! change this phrase into whatever you need BEFORE USING
                $branch = $1;
                print "branch $branch\n";
                my $modbranch = "$branch #1";
                print "modbranch $modbranch\n";
                $tree =~ s/$branch/$modbranch/;
                print "modified! $currentgroup\n";
        open (MODTREE, '>', "$currentgroup.aa.alignment_phyml_tree.txt.mod") or die "failed MODTREE handle";
        print MODTREE "$tree\n";
        close MODTREE;
        close TREE;
}


~
