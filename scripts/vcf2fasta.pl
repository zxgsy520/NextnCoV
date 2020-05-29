#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;

my %opt=(qual=>30, hom=>"Y", id=>"concensus",mindepth=>5);
GetOptions(\%opt,"ref:s", "vcf:s", "id:s", "qual:30", "GT:s", "depth:s", "mindepth:i");

my $usage=" Description: generate concensus sequence from *.vcf and ref.fna
 Options:
   -*ref      path for reference file
   -*vcf      path for vcf file
   --id       sequence name for output file
   --depth    depth_all.xls
   --mindepth for filtering low depth variations, default=5
   --qual     quality threshold for filtering snps from vcf file, default=30
   --hom      Genotype, hom or not, default=\"Y\"\n";

($opt{ref} && $opt{vcf}) ||die "$usage\n";
my %genotype;
$genotype{"Y"}='1/1';

my %virus; my @seq;
open FA,"<$opt{ref}" || die "ERROR: cannot find reference: $opt{ref}\n";
$/=">";
<FA>;
while(<FA>){
  chomp;
  my @l=split/\n/;
  my $id=(split/\s+/,$l[0])[0];
  shift @l;
  $virus{$id} = join("",@l);
  @seq=split//,join("",@l);
}
close FA;
$/="\n";

my %depth;
if($opt{depth}){
  open DP,"<$opt{depth}";
  while(<DP>){
    chomp;
    my @l=split/\s+/;
    $depth{$l[1]}=$l[2];
  }
  close DP;
}

open VCF,"<$opt{vcf}" || die "ERROR: cannot find vcf file: $opt{vcf}\n";
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
#MN908947	1548	.	G	A	61.455	PASS	pos1=1548;q1=61.455;pos2=1548;q2=61.455	GT:GQ	1/1:61.455
my $vcf_filt;
my $count=0;
while(<VCF>){
  if(/^#/){
    $vcf_filt .= $_;
    next;
  }
  next if ($genotype{$opt{hom}} && ($_=~/ 1\/1/));
  my @l=split/\t/;
  next if ($l[5]<$opt{qual});
  next if(length($l[3])>1 || length($l[4])>1);  #only snps
  next if($opt{depth} && $depth{$l[1]}<=$opt{mindepth});  #filter by mindepth=5
  my $tmp=$l[1];
  $l[1] --;
  $seq[$l[1]] = $l[4];
  $count ++;
  $vcf_filt .= $_;
  warn "[INFO] $l[3] -> $l[4]  $tmp\n";
}
close VCF;
warn "[INFO] Total $count site(s) changed!\n";

open VCF,">snps_Q$opt{qual}.vcf";
print VCF "$vcf_filt";
close VCF;

print ">$opt{id}\n".join("",@seq)."\n";
