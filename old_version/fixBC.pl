#!/usr/bin/perl
use strict;
use warnings;

my %bcs;
my %chl;
my %mit;
my $lib;

open F, "samtools view -h $ARGV[0] | " or die;
while(<F>){
	chomp;
	if($_ =~ /^@/){
		print "$_\n";
		next;
	}
	my @col = split("\t",$_);
	if($_ =~ /XA:Z:/ && $col[4] < 30){
		my @tag = grep(/NM:i:/, @col);
		my @xa = grep(/XA:Z:/, @col);
		my @edits = split(":",$tag[0]);
		my $ed = $edits[2];
		my @mapped = split(";",$xa[0]);
		my $near = 0;
		foreach my $aln (@mapped){
			my @vals = split(",",$aln);
			my $dif = $vals[$#vals] - $ed;
			if($dif < 3){
				$near++;
			}
		}
		if($near > 1){
			next;
		}
	}
	my $bc;
	my $tag;
	my $bcindex;
	my $hasbc = 0;
	for (my $i = 11; $i < @col; $i++){
		if($col[$i] =~ /CB:Z:/){
			$bcindex = $i;
			$bc = $col[$i];
			$hasbc++;
		}
		elsif($col[$i] =~ /RG:Z:/){
			my @id = split(":",$col[$i]);
			my @tag1 = split("_",$id[2]);
			$tag = join("",$tag1[0],$tag1[1],$tag1[2]);
		}
	}
	if($hasbc > 0){
		$bc =~ s/-1/-$tag/;

		# initialize Pt/Mt hashes
		if(!$chl{$bc}){
			$chl{$bc} = 0;
		}
		if(!$mit{$bc}){
			$mit{$bc} = 0;
		}

		# correct array index
		$col[$bcindex] = $bc;
		my $line = join("\t",@col);
		$bcs{$bc}++;
		$lib = $tag;

		# skip Pt/Mt on print
		##updating 042521 use Chr for the Rice
		if($col[2] =~ /ChrPt/){
			$chl{$bc}++;
		}elsif($col[2] =~ /ChrMt/){
			$mit{$bc}++;
		}else{
			print "$line\n";
		}
	}else{
		next;
	}

}
close F;

my $temp = $lib . "_bc_counts.txt";
open (my $t1, '>', $temp) or die;
print $t1 "cellID\ttotal\tnuclear\tlibrary\ttissue\n";

my @ids = sort {$bcs{$b} <=> $bcs{$a}} keys %bcs;
for (my $j = 0; $j < @ids; $j++){
	my $nuc = $bcs{$ids[$j]} - $chl{$ids[$j]} - $mit{$ids[$j]};
	my $tissue;
	if($lib =~ /panicle1/){
		$tissue = 'panicle';
	}elsif($lib =~ /panicle2/){
		$tissue = 'panicle';
	}elsif($lib =~ /seedling1/){
		$tissue = 'seedling';
	}elsif($lib =~ /seedling2/){
		$tissue = 'seedling';
	}elsif($lib =~ /root1/){
		$tissue = 'root';
	}elsif($lib =~ /root2/){
		$tissue = 'root';
	}elsif($lib =~ /Eseed1/){
		$tissue = 'Eseed';
	}elsif($lib =~ /Eseed2/){
		$tissue = 'Eseed';
	}elsif($lib =~ /crownroot1/){
		$tissue = 'crownroot';
	}elsif($lib =~ /crownroot2/){
		$tissue = 'crownroot';
	}elsif($lib =~ /root3/){
		$tissue = 'root';
	}elsif($lib =~ /root4/){
		$tissue = 'root';
	}elsif($lib =~ /seedling3/){
		$tissue = 'seedling';
	}elsif($lib =~ /seedling4/){
		$tissue = 'seedling';
	}elsif($lib =~ /Eseedling1/){
               	$tissue = 'Eseedling';
        }elsif($lib =~ /Eseedling2/){
                $tissue = 'Eseedling';
        }

	print $t1 "$ids[$j]\t$bcs{$ids[$j]}\t$nuc\t$lib\t$tissue\n";
}
close $t1;
