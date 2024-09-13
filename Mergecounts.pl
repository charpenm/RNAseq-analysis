
#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];
CREATE_LOG(\@ARGV);

my ($inpath, $id, $outpath)=@ARGV[1 .. $#ARGV];

# list all the count files
my @array=glob("$inpath/*_genes_htseq_counts.out");

print join ("\t", "NUMBER OF FILES:", scalar @array),"\n";

my %selected;
my %types;

my %counts;
my @headers1;

# generate a hash of the counts for each genes
foreach my $file (@array){
	my @temp=split /\//, $file;
	$temp[$#temp]=~s/_genes_htseq_counts.out//;
	my $name=$temp[$#temp];
	if(!-z $file){
		push @headers1, $name;
		open(IN, "$file")||die"IN $file\n";
		while(<IN>){
			chomp;
			my @split=split /\t/, $_;
			$counts{$split[0]}{$name}=$split[1];
		}
	}
}

open(OUT1, ">$outpath/$id\_expression_matrix_full.out");
print OUT1 join ("\t", "NAME", @headers1),"\n";
foreach my $gene (keys %counts){
	my @res1;
	foreach my $cell (@headers1){
		push @res1, $counts{$gene}{$cell};
	}
	print OUT1 join ("\t", $gene, @res1),"\n";
	@res1=();
}
close OUT1;


#########################################################
# creation of a log file or appending to an existing log file, the command and variables are being captured and printed into the log file with a timestamp
#########################################################
sub CREATE_LOG{
	my ($array)=(@_);
	my $date=localtime();
	if(-e "$project"){
		open ("IN", "$project")||die"IN $project\n";
		my @project_log=<IN>;
		close IN;
		open(LOG, ">$project");
		for my $n (0 .. $#project_log){
			print LOG $project_log[$n];
		}
	}
	else{
		open(LOG, ">$project");
	}
	print LOG "\n\n#########################################################\n#########################################################\n
	$date\n
#########################################################\n#########################################################\n
	\n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}
