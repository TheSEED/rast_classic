# -*- perl -*-
########################################################################
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
########################################################################

use FIG;
$fig = new FIG;

use constant LEFT   =>  0;
use constant RIGHT  =>  1;
use constant LEN    =>  2;
use constant FID    =>  3;


$0 =~ m/([^\/]+)$/;   $this_tool = $1;
$usage = "$this_tool  [-extend_5p=num] org_id  contigs  tbl_1 tbl_2 ... > gap_tbl";

if (not @ARGV) {
    die "\n\n\tusage: $usage\n\n";
}

$trouble = 0;
$extend_5prime_end = 0;

#...Process any flags and switches...
while ($ARGV[0] =~ m/^-/o) {
    if ($ARGV[$i] =~ m/-help/) {
	print STDERR "\n\n\tusage: $usage\n\n";
	exit(0);
    }
    elsif ($ARGV[$i] =~ m/-extend_5p=(\d+)/o) {
	$extend_5prime_end = $1;
	print STDERR "PEG 5-prime end will be extended $extend_5prime_end bp\n" if $ENV{VERBOSE};
    }
    else {
	$trouble = 1;
	warn "Unrecognized flag: $ARGV[$i]\n";
    }
    
    shift @ARGV;
}


#...Process Org-ID argument...
if (not @ARGV) {
    warn "\n\nNo argument list given\n\n";
    die "\n\n\tusage: $usage\n\n";
}
else {
    if ($ARGV[0] !~ m/^\d+\.\d+$/o) {
	$trouble = 1;
	warn "\nOrg-ID $ARGV[0] is malformed\n";
    }
    else {
	$org_id = shift @ARGV;
    }
}


#...Check existence of file arguments...
for ($i=0; $i < @ARGV; ++$i) {
    if (!-e $ARGV[$i]) {
	$trouble = 1;
	warn "\nERROR: File $ARGV[$i] does not exist\n";
    }
}
die "\n\nAborting due to invalid args\n\n   usage: $usage\n\n" if $trouble;

(($contigs_file = shift @ARGV) && (-s $contigs_file)) 
    || die "Contigs file $contigs_file has zero size\n\n\tusage: $usage\n\n";

((@tbls = @ARGV) > 0) || die "\n\tusage: $usage\n\n";

$len_of  = &load_contig_lens($contigs_file);
$regions = &load_regions($len_of, @tbls);
# die Dumper($features);

$gap_num = 0;
foreach $contig (sort keys %$len_of) {
    print STDERR "Processing $contig ...\n" if $ENV{VERBOSE};
    
    $x = $regions->{$contig};
    $x = defined($x) ? $x : [];
    
#   if ($x->[0]  > 1)                  { unshift @$x, [0, 0]; }
#   if ($x->[-1] > $len_of->{$contig}) { push    @$x, [1+$len_of->{$contig}, 1+$len_of->{$contig}]; }
    
    unshift @$x, [0, 0];
    push    @$x, [1+$len_of->{$contig}, 1+$len_of->{$contig}];
    
    for ($i=1; $i < @$x; ++$i) {
	$gap_beg = $x->[$i-1]->[RIGHT] + 1;
	$gap_end = $x->[$i]->[LEFT] - 1;
#	$gap_len = $gap_end - $gap_beg + 1;
	$gap_loc = "$contig\_$gap_beg\_$gap_end";
	
	if (($gap_beg < 1) || ($gap_end > $len_of->{$contig}) || ($gap_beg > $gap_end)) {
	    print STDERR "Skipping out-of-bounds location: $gap_loc\n" if $ENV{VERBOSE};
	    next;
	}
	
	print STDOUT (qq(fig|$org_id.gap.), (++$gap_num), qq(\t$gap_loc\n));
    }
}
print STDERR "\n$this_tool done\n\n" if $ENV{VERBOSE};
exit(0);


sub load_contig_lens
{
    my ($contigs_file) = @_;
    
    my $num_contigs = 0;
    my $num_chars   = 0;
    
    my $len_of = {};
    
    print STDERR "Loading $contigs_file ...\n" if $ENV{VERBOSE};
    open (CONTIGS, "<$contigs_file") or die "could not open $contigs_file to read";
    while (($id, $seqP) = &FIG::read_fasta_record(\*CONTIGS)) {
	++$num_contigs;
	$num_chars += $len_of->{$id} = length($$seqP);
    }
    close(CONTIGS) or die "could not close $contigs_file";
    print STDERR "Loaded $num_contigs contig lengths (total $num_chars chars) from $contigs_file\n\n" if $ENV{VERBOSE};

    return $len_of;
}

sub load_regions
{
    my ($len_of, @tbl_files) = @_;
    my ($fid, $org, $locus, $contig, $beg, $end, $len);
    
    my $regions  = {};
    foreach my $tbl_file (@tbl_files) {
	my $num_fids = 0;
	open(TBL,"<$tbl_file") || die "could not open $tbl_file";
	print STDERR "Loading $tbl_file ...\n" if $ENV{VERBOSE};
	
	while (defined($entry = <TBL>))
	{
	    ++$num_fids;
	    chomp $entry;
	    if ($entry =~ m/^(\S+)\s+(\S+)/o) {
		($fid, $locus) = ($1, $2);
		
		$fid_type = "";
		if ($fid =~ m/^fig\|\d+\.\d+\.(rna|peg|orf)\.\d+$/) {
		    $fid_type = $1;
		}
		else {
		    die "Unrecognized FID type: $fid"; 
		}
		
		($contig, $beg, $end) = $fig->boundaries_of($locus);
		
		if ($extend_5prime_end && ($fid_type eq qq(peg))) {
		    print STDERR "$fid:\t$locus\tbeg = $beg --> " if $ENV{VERBOSE};
		    $beg -= ($beg < $end) ? +$extend_5prime_end : -$extend_5prime_end;
		    $beg = &FIG::max(1, &FIG::min($beg, $len_of->{$contig}));
		    print STDERR "$beg\n" if $ENV{VERBOSE};
		}
		$len = 1 + abs($end-$beg);
		
		unless (defined($regions->{$contig})) { $regions->{$contig} = []; }
		push(@ { $regions->{$contig} }
		     , [ $left = &FIG::min($beg, $end), $right = &FIG::max($beg, $end), (1 + $right - $left), $fid ]
		     );
	    }
	    else
	    {
		print STDERR "Skipping invalid entry $entry\n";
	    }
	}
	print STDERR "Loaded $num_fids features from $tbl_file\n\n" if $ENV{VERBOSE};
	close(TBL);
    }
    
    foreach my $contig (keys(%$regions))
    {
	my $x = [ sort { ($a->[LEFT] <=> $b->[LEFT]) || ($b->[RIGHT] <=> $a->[RIGHT])
			 }  @ { $regions->{$contig} }
		  ];
    
	for (my $i=1; $i < @$x; ++$i)
	{
	    while (($i < @$x) && ($x->[$i-1]->[RIGHT] >= $x->[$i]->[LEFT]))
	    {
		print STDERR "Merging $i:\t$x->[$i-1]->[FID]:$x->[$i-1]->[LEN]\t$x->[$i]->[FID]:$x->[$i]->[LEN]"
		    if ($ENV{VERBOSE} && ($ENV{VERBOSE} > 1));
		
		$x->[$i-1]->[RIGHT]  = &FIG::max( $x->[$i-1]->[RIGHT], $x->[$i]->[RIGHT] );
		my $new_len = 1 + ($x->[$i-1]->[RIGHT] - $x->[$i-1]->[LEFT]);

		print STDERR "\t--- new_len = $new_len\n"
		    if ($ENV{VERBOSE} && ($ENV{VERBOSE} > 1));
		
		$x->[$i-1]->[FID] .= ",$x->[$i]->[FID]";
		$x->[$i-1]->[LEN]  =   $new_len;
		
		splice @$x, $i, 1;
	    }
	}
	
	$regions->{$contig} = $x;
	print STDERR "\nAfter mergers, $contig has ", (scalar @$x), " regions\n\n" if $ENV{VERBOSE}; 
    }

    return $regions;
}
