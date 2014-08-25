open (GROUPS, '<', $ARGV[0]);
open (OUTPUT, '>', "$ARGV[0].annotated");

## note: we are using an example to pull Biomart gene names of Galgal4.7 genome. 
## use Ensembl Biomart, output Entrez GeneID, transcript ID, protein ID, gene name as four columns, as in the sample file ./sample/Biomart_Gallus_full.txt

while(my $group = <GROUPS>){
        chomp $group;
        $group =~ s/GGAL\|//g; ## here we process the orthoMCL outputs. You may modify or may not need these. 
        $group =~ s/PICA\|//g;
        $group =~ s/TGUT\|//g;
        $group =~ s/CBRA\|//g;
        $group =~ s/PHUM\|//g;
        my @fields = split " ", $group;
        my @fields2 = split ":", $group;
        my $groupnumber = $fields2[0];
        my @chickname;
        foreach $item (@fields){
                if ($item =~ /ENSG/) { ## change your item keyword like ENSG according to the species you choose
                        open (GALLUS, '<', "Biomart_Gallus_full.txt");
                        while (my $chickBiomart = <GALLUS>){
                                chomp $chickBiomart;
                                @chickBiomartList = split "\t", $chickBiomart;
                                if ($chickBiomartList[2] eq $item) {
                                        push (@chickname,$chickBiomartList[3]);
                                }
                        }
                        close GALLUS;
                }
        }
        print OUTPUT "$group\t@chickname\n";
}
