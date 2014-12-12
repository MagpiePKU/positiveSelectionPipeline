## 20141213
## Zhang Yi, Peking University
## synapse@pku.edu.cn
## 
## use for computing codon-wise and amino acid-wise local entropy across species using existing alignments
## to do this you need a PAL2NAL alignment file and a config file denoting the unique naming of each species e.g. 
## PAL2NAL:
## SPECIESX_geneA
## ATGATGATG
## SPECIESY_geneA
## GTAGTAGTA
## SPECIESZ_geneA
## GTAGTAATG
## and the config file:
## X (tab) SPECIESX (tab) groupAlpha
## Y (tab) SPECIESY (tab) groupBeta
## Z (tab) SPECIESZ (tab) groupBeta
## The program computes entropy of each site within (A_ and B_entropy) and across ('totalentropy') the group (alpha/beta) of species (or population) 
## The current program will not run if gene name having ambiguities, for example, using P.humilis data together with G.gallus. 


open (INPUT, "<", $ARGV[0]);
open (OUTPUTDNA, ">>", "$ARGV[1].DNA");
open (OUTPUTAA, ">>", "$ARGV[1].PROTEIN");
open (CONFIG, "<", $ARGV[2]);

print OUTPUTDNA "CurrentLength\tMaxTotalCodon\tTotalEntropy\tMaxA_Codon\tAprob{maxAcodon}\tA_Entropy\tMaxB_Codon\tBprob{maxBcodon}\tB_Entropy\tEntropy_difference\n";

print OUTPUTAA "CurrentLength\tMaxTotal_AA\tTotalEntropy\tMaxA_AA\tAprob{maxAcodon}\tA_Entropy\tMaxB_AA\tBprob{maxBcodon}\tB_Entropy\tEntropy_difference\n";

########## Read config 

while (<CONFIG>){
        chomp $_;
        @array = split "\t", $_;
        if ($array[0] =~ /[^\s+]/){
        $configspecies{$array[0]} = $array[1];
        if($array[2]=~"groupA"){
                $groupAspecies{$array[0]} = $array[1];
        }
        if($array[2]=~"groupB"){
                $groupBspecies{$array[0]} = $array[1];
        }
        }

}

########## Read DNA and translate into protein

%usedspecies;
$useline = 1;
%sequencearray;
%proteinarray;
$jump = 0;
$previousspecies = "tempspecies";



while (<INPUT>){
        chomp $_;
        if ($_ =~ /^\s/) {next;} ## not empty line
        $line = $_; ## finish processing line

## check whether is header
        if ($line =~ /\d/) {
                $header = 1;
        } else {$header = 0;}

        if ($line =~/-/) {
                $header = 1; $nextheader = 1; next;} ## temporary method to solve the P.humilis naming issue

        if ($header == 1) {
                foreach $speciesnotation (sort keys(%configspecies)){
                        if ($line =~ /$configspecies{$speciesnotation}/){
                                if (defined $usedspecies{$speciesnotation}){
                                        $nextheader = 1;last;
                                } else {$nextheader = 0;
                                        $usedspecies{$speciesnotation} = 0; 
                                        $usedspecies{$previousspecies} = 1; ## record the finishing status
                                        $currentspecies = $speciesnotation; ## record the current new species
                                        last; ## quit cycle
                                }
                        }
                }
        }

## 

        if ($header == 0 && $nextheader == 0){ ## not header and not jump
                if ($usedspecies{$currentspecies}==0){
                        $sequencearray{$currentspecies} .= $line; ## concatanate
                        $previousspecies = $currentspecies; ## record the current species
                }
        }
}

my %codon2aa = qw(
    TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
    TAC  Y  TAT  Y  TAA  _  TAG  _  TGC  C  TGT  C  TGA  _  TGG  W
    CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
    CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
    ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
    AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
    GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
    GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
    );


foreach $speciesnotation(sort keys (%configspecies)){
        $sequencearray{$speciesnotation} =~ tr/tcga/TCGA/; ## transfer to capital
        $dnasequence = $sequencearray{$speciesnotation};
        chomp $dnasequence;
        while ($dnasequence =~ /(...)/g){
                $codon = $1;
                $proteinarray{$speciesnotation} .= $codon2aa{$1} or die "aa code bad";
        }
}


########## DNA

$lengthofstring = 3;
$maxlength = 0;

foreach $speciesnotation(sort keys (%configspecies)){
        $dnasequence = $sequencearray{$speciesnotation};
        chomp $dnasequence;
        if(length($dnasequence)>$maxlength){
                $maxlength = length($dnasequence);}
}


while ($lengthofstring >= 3){

my %groupAnow;
my %groupBnow;


foreach $speciesnotation (sort keys(%groupAspecies)){ ## cycle through all species
        if ($usedspecies{$speciesnotation} == 1){ ## if read this sequence
                $sequencearray{$speciesnotation} =~ s/^(...)//; ## if still has sequence
                if (length($1)==3){
                        $groupAnow{$speciesnotation} = $1; ## catch sequence
                        $lengthofstring = length($sequencearray{$speciesnotation});## calculate remaining length
                }
                else {last;}
        }
}

foreach $speciesnotation (sort keys(%groupBspecies)){
        if ($usedspecies{$speciesnotation} == 1){ ## if read this sequence
                $sequencearray{$speciesnotation} =~ s/^(...)//; ## if still has sequence
                if (length($1)==3){
                        $groupBnow{$speciesnotation} = $1; ## catch sequence
                        $lengthofstring = length($sequencearray{$speciesnotation});## calculate remaining length
                }
                else {last;}
        }
}


my %Aarray;
my %Barray;
my %totalarray;
my $Atotalnumber = 0;
my $Btotalnumber = 0;
my $totalnumber = 0;


foreach $speciesnotation (sort keys(%groupAspecies)){
        if (defined $groupAnow{$speciesnotation}){
                $Aarray{$groupAnow{$speciesnotation}} += 1;
                $totalarray{$groupAnow{$speciesnotation}} += 1;
                $Atotalnumber += 1;
                $totalnumer += 1;
        }
}

foreach $speciesnotation (sort keys(%groupBspecies)){
        if (defined $groupBnow{$speciesnotation}){
                $Barray{$groupBnow{$speciesnotation}} += 1;
                $totalarray{$groupBnow{$speciesnotation}} += 1;
                $Btotalnumber += 1;
                $totalnumer += 1;
        }
}

@Akeys = keys (%Aarray);
@Bkeys = keys (%Barray);

my $maxAcodon;
my $maxBcodon;
my $maxtotalcodon;

my $previouscodon = 0;
foreach $codon (sort keys (%Aarray)){
                                if ($previouscodon == 0) { 
                                        $previouscodon = $codon;
                                        $maxAcodon = $codon;}
                                if ($Aarray{$codon} > $Aarray{$previouscodon}){
                                                $maxAcodon = $codon;
                                }
}

$previouscodon = 0;
foreach $codon (sort keys (%Barray)){
                                if ($previouscodon == 0) { 
                                        $previouscodon = $codon;
                                        $maxBcodon = $codon;}
                                if ($Barray{$codon} > $Barray{$previouscodon}){
                                                $maxBcodon = $codon;
                                }
}

$previouscodon = 0;
foreach $codon (sort keys (%totalarray)){
                                if ($previouscodon == 0) { 
                                        $previouscodon = $codon;
                                        $maxtotalcodon = $codon;}
                                if ($totalarray{$codon} > $totalarray{$previouscodon}){
                                                $maxtotalcodon = $codon;
                                }
}


my %Aprob;
my %Bprob;
my $totalprob;
my $Aentropy = 0;
my $Bentropy = 0;
my $totalentropy = 0;

if ($Atotalnumber > 0 && $Btotalnumber > 0){
                $totalnumber = $Atotalnumber + $Btotalnumber;}

foreach $codon (sort keys (%Aarray)){
                        $Aprob{$codon} = $Aarray{$codon}/$Atotalnumber;
                        $Aentropy += $Aprob{$codon} * log($Aprob{$codon});
}

foreach $codon (sort keys (%Barray)){
                        $Bprob{$codon} = $Barray{$codon}/$Btotalnumber;
                        $Bentropy += $Bprob{$codon} * log($Bprob{$codon});
}

foreach $codon (sort keys (%totalarray)){
                        $totalprob{$codon} = $totalarray{$codon}/$totalnumber;
                        $totalentropy += $totalprob{$codon} * log($totalprob{$codon});
}

$currentlength = $maxlength - $lengthofstring;

my $entropydifference = abs(abs($Aentropy) - abs($Bentropy));

@write = ($currentlength, $maxtotalcodon, $totalentropy, $maxAcodon, $Aprob{$maxAcodon}, $Aentropy, $maxBcodon, $Bprob{$maxBcodon}, $Bentropy, $entropydifference);

print OUTPUTDNA join ("\t", @write), "\n"; 


}


## above working 





########## Protein

$lengthofstring = 1;

$maxlength = 0;

foreach $speciesnotation(sort keys (%configspecies)){
        $proteinsequence = $proteinarray{$speciesnotation};
        chomp $proteinsequence;
        if(length($proteinsequence)>$maxlength){
                $maxlength = length($dnasequence);}
}


while ($lengthofstring > 0){

my %groupAnow;
my %groupBnow;


foreach $speciesnotation (sort keys(%groupAspecies)){ ## cycle through all species
        if ($usedspecies{$speciesnotation} == 1){ ## if read this sequence
                $proteinarray{$speciesnotation} =~ s/^(.)//; ## if still has sequence
                if (length($1)==1){
                        $groupAnow{$speciesnotation} = $1; ## catch sequence
                        $lengthofstring = length($proteinarray{$speciesnotation});## calculate remaining length
                }
                else {last;}
        }
}

foreach $speciesnotation (sort keys(%groupBspecies)){
        if ($usedspecies{$speciesnotation} == 1){ ## if read this sequence
                $proteinarray{$speciesnotation} =~ s/^(.)//; ## if still has sequence
                if (length($1)==1){
                        $groupBnow{$speciesnotation} = $1; ## catch sequence
                        $lengthofstring = length($proteinarray{$speciesnotation});## calculate remaining length
                }
                else {last;}
        }
}


my %Aarray;
my %Barray;
my %totalarray;
my $Atotalnumber = 0;
my $Btotalnumber = 0;
my $totalnumber = 0;


foreach $speciesnotation (sort keys(%groupAspecies)){
        if (defined $groupAnow{$speciesnotation}){
                $Aarray{$groupAnow{$speciesnotation}} += 1;
                $totalarray{$groupAnow{$speciesnotation}} += 1;
                $Atotalnumber += 1;
                $totalnumer += 1;
        }
}

foreach $speciesnotation (sort keys(%groupBspecies)){
        if (defined $groupBnow{$speciesnotation}){
                $Barray{$groupBnow{$speciesnotation}} += 1;
                $totalarray{$groupBnow{$speciesnotation}} += 1;
                $Btotalnumber += 1;
                $totalnumer += 1;
        }
}

@Akeys = keys (%Aarray);
@Bkeys = keys (%Barray);

my $maxA_amino_acid;
my $maxB_amino_acid;
my $maxtotal_amino_acid;

my $previous_amino_acid = 0;
foreach $amino_acid (sort keys (%Aarray)){
                                if ($previous_amino_acid == 0) { 
                                        $previous_amino_acid = $amino_acid;
                                        $maxA_amino_acid = $amino_acid;}
                                if ($Aarray{$amino_acid} > $Aarray{$previous_amino_acid}){
                                                $maxA_amino_acid = $amino_acid;
                                }
}

$previous_amino_acid = 0;
foreach $amino_acid (sort keys (%Barray)){
                                if ($previous_amino_acid == 0) { 
                                        $previous_amino_acid = $amino_acid;
                                        $maxB_amino_acid = $amino_acid;}
                                if ($Barray{$amino_acid} > $Barray{$previous_amino_acid}){
                                                $maxB_amino_acid = $amino_acid;
                                }
}

$previous_amino_acid = 0;
foreach $amino_acid (sort keys (%totalarray)){
                                if ($previous_amino_acid == 0) { 
                                        $previous_amino_acid = $amino_acid;
                                        $maxtotal_amino_acid = $amino_acid;}
                                if ($totalarray{$amino_acid} > $totalarray{$previous_amino_acid}){
                                                $maxtotal_amino_acid = $amino_acid;
                                }
}


my %Aprob;
my %Bprob;
my $totalprob;
my $Aentropy = 0;
my $Bentropy = 0;
my $totalentropy = 0;

if ($Atotalnumber > 0 && $Btotalnumber > 0){
                $totalnumber = $Atotalnumber + $Btotalnumber;}

foreach $amino_acid (sort keys (%Aarray)){
                        $Aprob{$amino_acid} = $Aarray{$amino_acid}/$Atotalnumber;
                        $Aentropy += $Aprob{$amino_acid} * log($Aprob{$amino_acid});
}

foreach $amino_acid (sort keys (%Barray)){
                        $Bprob{$amino_acid} = $Barray{$amino_acid}/$Btotalnumber;
                        $Bentropy += $Bprob{$amino_acid} * log($Bprob{$amino_acid});
}

foreach $amino_acid (sort keys (%totalarray)){
                        $totalprob{$amino_acid} = $totalarray{$amino_acid}/$totalnumber;
                        $totalentropy += $totalprob{$amino_acid} * log($totalprob{$amino_acid});
}

$currentlength = $maxlength - $lengthofstring;

my $entropydifference = abs(abs($Aentropy) - abs($Bentropy));

@write = ($currentlength, $maxtotal_amino_acid, $totalentropy, $maxA_amino_acid, $Aprob{$maxA_amino_acid}, $Aentropy, $maxB_amino_acid, $Bprob{$maxB_amino_acid}, $Bentropy, $entropydifference);

print OUTPUTAA join ("\t", @write), "\n"; 


}




