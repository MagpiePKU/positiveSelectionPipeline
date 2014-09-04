# The program switches one name to another between annotation versions
# input: 1. Groups, any orthoMCL-generated file
# input: 2. Name matching list, like:
# xxxx(old name)  yyyy(new name)
# the keyword (evm.xxxx) should be changed to your keyword if you use it


open (GROUPS,'<',$ARGV[0]);
open (LIST, '<', $ARGV[1]);
open (MODIFIEDGROUPS, '>',"$ARGV[0].mod");
my %list;

while (<LIST>){
        chomp $_;
        my @fields = split "\t", $_;
        $list{$fields[0]} = $fields[1];
}

while (<GROUPS>){
        chomp $_;
        while(/(evm.TU.scaffold_\d+.\d+)|(evm.TU.repeat_\d+.\d+)/){ #change the keyword before use
                my $PPIG=$list{$1};
                s/$1/$PPIG/;
        }
        print MODIFIEDGROUPS "$_\n";
}

close GROUPS;
close LIST;
close MODIFIEDGROUPS;
