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

# use File::Basename;   $this_tool_name = basename($0);

use File::Basename;
use strict;
use warnings;
use Data::Dumper;

#use FIG;
use gjoseqlib;


use Pod::Text;
if ((@ARGV == 1) && ($ARGV[0] =~ m/-help/))  {
    pod2text($0);  exit(0);
}

=pod

=head1 NAME

reformat_contigs

=head1 SYNOPSIS
    
reformat_contigs [-help] [-logfile=logfilename] [-v(erbose)] [-split(=len)] [-keep] [-width=line_width] [-[no]renumber] [-min=min_length]  < fasta_in  > fasta_out

reformat_contigs [-help] [-logfile=logfilename] [-v(erbose)] [-split(=len)] [-keep] [-width=line_width] [-[no]renumber] [-min=min_length]    fasta_in    fasta_out

=head1 DESCRIPTION

Reformats a contigs file to a uniform number of characters per line
(default length of 50), checking for non-IUPAC characters.

Default behavior is to zero-pad the digits of IDs that are of the form
'Contig', immediately followed by one or more digits to four digits.

Optionally renumbers contig IDs to format "PrefixNNNN", where 'Prefix'
defaults to "Contig", "NNNN" are zero-padded digits, and the number of digits
defaults to 4.

Optionally splits contigs at runs of ambig chars (default 4 or more ambigs),
and homopolymer runs (default 50 or more of the same character).
If the script is invoked using the two-argument form or explicit
'-input' and '-output' arguments, a file 'scaffold_map' is written
in the same directory as the output file.

Optionally discards contigs that are shorter or longer than some threshold.

=head1 COMMAND-LINE OPTIONS
		
    --split       Split contigs on runs of ambigs longer than len (def=4)
                    (writes scaffold-map to output path if there is an explicit 'fasta_out' arg)
    --keep        Appends original contig ID to new ID if -renumber
    --width       line width (def=50)
    --renumber    Renumbers contigs starting from Contig0001
    --norenumber  Does not renumber contigs  (Default)
    --renumber-map=filename	Save the a mapping of renumbered contig names to filename
    --renumber-digits=ndigits Number of digits used in renumbered contig names
    --renumber-prefix=prefix  Prefix to use instead of Contig in renumbered contig names
    --remove-duplicates   Remove duplicate sequences
    --duplicates-file     File in which to save information about removed duplicates
    --min         Minimum acceptable contig length (def=0)
    --max         Maximum acceptable contig length (def=999999999999)
    --logfile     Name of the optional logfile

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut
    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Main variable declarations...
#-----------------------------------------------------------------------
my $help;
my $verbose;

my $fasta_in;
my $fasta_out;

my $keep        = 0;
my $split       = 0;
my $longrepeat  = 50;
my $width       = 50;

my $renumber    = 0;
my $renumber_map;
my $renumber_digits = 4;
my $renumber_prefix = "Contig";

my $remove_duplicates = 0;
my $duplicates_file;
my $duplicates_fh;
my $duplicates_removed = 0;

my $min_length  = 0;
my $max_length  = 999999999999;
my $max_bad     = 500;

my $short_chars = 0;
my $long_chars  = 0;
my $nx_chars    = 0;

my $logfile_fh;
my $logfilename = undef;


use Getopt::Long;
my $rc = GetOptions('help'         => \$help,
		    'verbose'      => \$verbose,
		    'logfile=s'    => \$logfilename,
		    'input=s'      => \$fasta_in,
		    'output=s'     => \$fasta_out,
		    'keep'         => \$keep,
		    'split:4'      => \$split,
		    'split-homopol=i' => \$longrepeat,
		    'max-badchar=i'   => \$max_bad,
		    'width=i'      => \$width,
		    'min=i'        =>  $min_length,
		    'max=i'        =>  $max_length,
		    'renumber!'         => \$renumber,
		    'renumber-map=s'    => \$renumber_map,
		    'renumber-digits=i' => \$renumber_digits,
		    'renumber-prefix=s' => \$renumber_prefix,
		    'remove-duplicates' => \$remove_duplicates,
		    'duplicates-file=s' => \$duplicates_file,
		    );
if ($verbose) { $ENV{VERBOSE} = 1; }


if (!$rc || $help || (@ARGV != 2 && @ARGV != 0)) {
    seek(DATA, 0, 0);
    while (<DATA>) {
	last if /^=head1 COMMAND-LINE /;
    }
    while (<DATA>) {
	last if (/^=/);
	print $_;
    }
    exit($help ? 0 : 1);
}
my $trouble = 0;



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Grandfather in the "two-argument" invocation mode...
#-----------------------------------------------------------------------
my $basepath;
my $fh_in  = \*STDIN;
my $fh_out = \*STDOUT;
if (@ARGV == 2) {
    my $fasta_in  = shift;
    my $fasta_out = shift;
    
    if (!-s $fasta_in) {
	$trouble = 1;
	warn "Input file \'$fasta_in\' does not exist\n";
    }
    else {
	open($fh_in,  q(<), $fasta_in)  || die "Could not read-open \'$fasta_in\' --- ERROR=$!";
	open($fh_out, q(>), $fasta_out) || die "Could not write-open \'$fasta_out\' --- ERROR=$!";
	$basepath = dirname($fasta_out);
	
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Open various auxilliary files (if requested)
#-----------------------------------------------------------------------
$logfile_fh = \*STDERR;
if (defined($logfilename)) {
    if (not open($logfile_fh, q(>), $logfilename)) {
	$trouble = 1;
	warn "ERROR=$!: Could not write-open logfile=\'$logfilename\'";
    }
}

my $renumber_map_fh;
if ($renumber) {
    if (not open($renumber_map_fh, q(>), $renumber_map)) {
	$trouble = 1;
	warn "ERROR=$!: Cannot write-open renumber-map file \'$renumber_map\'";
    }
}

if ($duplicates_file) {
    if (not open($duplicates_fh, ">", $duplicates_file)) {
	$trouble = 1;
	warn "ERROR=$1: Cannot write-open duplicates file \'$duplicates_file\'";
    }
}

if ($split && $basepath) {
    if (not open(SCAFFOLD_MAP, ">$basepath/scaffold.map")) {
	$trouble = 1;
	warn "ERROR=$1: Could not write-open scaffold-map file \'$basepath/scaffold.map\'";
    }
    elsif (not open(SCAFFOLD_MAP, ">/dev/null")) {
	$trouble = 1;
	warn "Could not write-open scaffold-map filehandle for \'/dev/null\'";
    }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Abort if any problems were detected
#-----------------------------------------------------------------------
if ($trouble || @ARGV) {
    warn qq(There were invalid arguments: ), join(" ", @ARGV), qq(\n\n);
    pod2text($0);
    die "aborted";
}
#=======================================================================



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Main body of code
#-----------------------------------------------------------------------
my $num_contigs   = 0;
my $num_short = 0;
my $num_long  = 0;
my $badchars = 0;
my $num_bad_contigs   = 0;
my $max_contig_id_len = 0;

my %id_by_content;


#while (my($head, $seqP) = &get_fasta_record($fh_in)) {

while (my @record = &gjoseqlib::read_next_fasta_seq($fh_in)) {
    ++$num_contigs;
    my ($contig_id, $comment, $seq) = @record;
    
    $seq = lc($seq);
    my $len = length($seq);
    
    unless (defined($contig_id) && $len) {
	++$num_bad_contigs;
	print $logfile_fh "No contig-ID for record num=$num_contigs, line=$., len=$len" unless (defined($contig_id));
	print $logfile_fh "Zero-length record num=$num_contigs, line=$., contig_id=$contig_id" unless ($len);
	next;
    }
    
    my $prev_id = $id_by_content{$seq};
    
    my $head;
    if ($renumber) {
	my $orig = $contig_id . ($comment ? qq(\t$comment) : q());
	my $contig_id = $renumber_prefix . ("0" x ($renumber_digits - length($num_contigs))) . $num_contigs; 
	if ($keep) {
	    $head = "$contig_id\t$orig";
	}
	else {
	    $head = $contig_id;
	}
	print $renumber_map_fh "$orig\t$head\n" if $renumber;
    }
    else {
	if (not $contig_id) {
	    $trouble = 1;
	    ++$num_bad_contigs;
	    print $logfile_fh "Record $. has leading whitespace or no contig name\n";
	    next;
	}
    }
    
    
    if ($contig_id =~ m/\,/o) {
	++$num_bad_contigs;
	$trouble = 1;
	print STDERR "Record $. has a comma embedded in the contig name\n";
	next;
    }
    
    if ((my $l = length($contig_id)) > $max_contig_id_len) {
	$max_contig_id_len = $l;
    }

    if ($prev_id) {
	if ($remove_duplicates) {
	    print $duplicates_fh "$prev_id\t$contig_id\n" if $duplicates_fh;
	    $duplicates_removed++;
	    next;
	}
    }
    else {
	$id_by_content{$seq} = $contig_id;
    }
    
    my $accept;
    my @ambig_runs = ();
    if ($split) {
	$_ =  $seq;
	while ($_ =~ m/([nbdhvrykmswx]{$split,}|a{$longrepeat,}|c{$longrepeat,}|g{$longrepeat,}|t{$longrepeat,})/gio) {
	    my $run = $1;
	    $nx_chars += length($run);
	    push @ambig_runs, $run;
	}

	my $runs;
	if (@ambig_runs > 1) { $runs = 'runs'; } else { $runs = 'run'; }
	
	if (defined($ENV{VERBOSE}) && @ambig_runs) {
	    if ($ENV{VERBOSE} == 1) {
		print STDERR "$head contains "
		    , (scalar @ambig_runs)
		    , " long $runs of ambiguity characters separating subcontigs, with run-lengths "
		    , join(qq(, ), (map { length($_) } @ambig_runs)), "\n"
		    if (@ambig_runs && $ENV{VERBOSE});
	    }
	    elsif ($ENV{VERBOSE} > 1) {
		print STDERR "$head contains "
		    , (scalar @ambig_runs)
		    , " long $runs of ambiguity characters separating subcontigs: "
		    , join(qq(, ), @ambig_runs), "\n"
		    if (@ambig_runs && $ENV{VERBOSE});
	    }
	}		
    }
    
    print SCAFFOLD_MAP "$contig_id\n" if $split;
    if ($split && @ambig_runs) {
	my $last_pos   = 0;
	my $subcon_num = 0;
	my($subcontig, $subcontig_id);
	my($prefix, $bridge, $suffix) = ("", "", "");
	while ($seq =~ m/[nbdhvrykmswx]{$split,}|a{$longrepeat,}|c{$longrepeat,}|g{$longrepeat,}|t{$longrepeat,}/gio) {
	    ($prefix, $bridge, $suffix) = ($`, $&, $');
	    print STDERR "$bridge\n";
	    
	    $accept = 1;
	    if (length($prefix)) {
		++$subcon_num;
		$subcontig_id = "$contig_id.$subcon_num";
		
		$subcontig = substr($seq, $last_pos, length($prefix) - $last_pos);
		print SCAFFOLD_MAP "$subcontig_id\t", $last_pos, "\n";
		
		if (($len = length($subcontig)) < $min_length) {		    
		    print STDERR "   skipping len=$len $subcontig_id\n" if ($ENV{VERBOSE});
		    ++$num_short;
		    $short_chars += $len;
		    $accept = 0;
		}
		
		if (($len = length($subcontig)) > $max_length) {		    
		    print STDERR "   skipping len=$len $subcontig_id\n" if ($ENV{VERBOSE});
		    ++$num_long;
		    $long_chars += $len;
		    $accept = 0;
		}
		
		print STDERR "   accepting prefix len=$len $subcontig_id\n" if $ENV{VERBOSE};
		&display_id_and_seq($subcontig_id, \$subcontig, $fh_out) if $accept;
	    }
	    $last_pos = pos($seq);
	}

	$accept = 1;
	if ($suffix) {
	    ++$subcon_num;
	    $subcontig_id = "$contig_id.$subcon_num";
	    
	    $subcontig = $suffix;
	    print SCAFFOLD_MAP "$subcontig_id\t", $last_pos,"\n";
	    
	    if (($len = length($subcontig)) < $min_length) {		    
		print STDERR "   skipping len=$len $subcontig_id\n" if ($ENV{VERBOSE});
		++$num_short;
		$short_chars += $len;
		$accept = 0;
	    }
	    
	    if (($len = length($subcontig)) > $max_length) {		    
		print STDERR "   skipping len=$len $subcontig_id\n" if ($ENV{VERBOSE});
		++$num_long;
		$long_chars += $len;
		$accept = 0;
	    }
	    
	    print STDERR "   accepting suffix len=$len $subcontig_id\n" if $ENV{VERBOSE};
	    &display_id_and_seq($subcontig_id, \$subcontig, $fh_out) if $accept;
	}
    }
    else {
	$accept = 1;
	
	if (($len = length($seq)) < $min_length) {		    
	    print STDERR "   skipping len=$len $head\n" if ($ENV{VERBOSE});
	    ++$num_short;
	    $short_chars += $len;
	    $accept = 0;
	}
	
	if (($len = length($seq)) > $max_length) {		    
	    print STDERR "   skipping len=$len $head\n" if ($ENV{VERBOSE});
	    ++$num_long;
	    $long_chars += $len;
	    $accept = 0;
	}
	
	&display_id_and_seq( $contig_id, \$seq, $fh_out ) if $accept;
    }
    print SCAFFOLD_MAP "//\n" if $split;


#...Abort if error-thresholds exceeded...    
    if ($badchars > $max_bad) {
	warn "Aborting reformat: reached $badchars bad characters\n";
	last;
    }
    
    if ($num_bad_contigs > $max_bad) {
	warn "Aborting reformat: reached $num_bad_contigs bad contigs\n";
	last;
    }
}

my ($s, $sa, $sb);
print STDERR "\n" if ($num_bad_contigs);
print STDERR "max id length $max_contig_id_len\n";

$s = ($duplicates_removed == 1) ? qq() : qq(s);
print STDERR "removed $duplicates_removed duplicate$s.\n" if ($duplicates_removed);

$s = ($num_bad_contigs == 1) ? qq() : qq(s);
print STDERR "skipped $num_bad_contigs bad contig$s.\n" if ($num_bad_contigs);

$s = ($badchars == 1) ? qq() : qq(s);
print STDERR "skipped $badchars invalid char$s.\n"   if ($badchars);

$s = ($nx_chars == 1) ? qq() : qq(s);
print STDERR "skipped $nx_chars ambiguity char$s.\n" if ($nx_chars);

$sa = ($short_chars == 1) ? qq() : qq(s);
$sb = ($num_short == 1) ? qq() : qq(s);
print STDERR "skipped $short_chars char$sa in $num_short contig$sb shorter than $min_length bp.\n" if ($num_short);

$sa = ($long_chars == 1) ? qq() : qq(s);
$sb = ($num_long == 1) ? qq() : qq(s);
print STDERR "skipped $long_chars char$sa in $num_long contig$sb longer than $max_length bp.\n"   if ($num_long);
print STDERR "\n" if ($num_bad_contigs);

if (defined($logfilename)) {
    close $logfile_fh;
}

exit($trouble || $num_bad_contigs);

#=======================================================================

sub display_id_and_seq {
    my( $id, $seq, $fh ) = @_;
    my ( $i, $n, $ln );
    
    if (! defined($fh) )  { $fh = \*STDOUT; }
    
    print $fh ">$id\n";
    
    $n = length($$seq);
#   confess "zero-length sequence ???" if ( (! defined($n)) || ($n == 0) );
    for ($i=0; ($i < $n); $i += $width) {
	if (($i + $width) <= $n) {
	    $ln = substr($$seq,$i,$width);
	}
	else {
	    $ln = substr($$seq,$i,($n-$i));
	}
	
	print $fh "$ln\n";
    }
}


my $first_line  = 1;
my $input_eol_marker = $/;
sub get_fasta_record {
    my ( $fh ) = @_;
    my ( $old_eol, $entry, @lines, $head, $seq);
    
    if (not defined($fh))  { $fh = \*STDIN; }
    $old_eol = $/;
    $/ = "$input_eol_marker>";
    
    my @record = ();
    if (defined($entry = <$fh>)) {
	chomp $entry;
	@lines =  split( /$input_eol_marker/, $entry );
	while (@lines and (not defined($lines[0])))  { shift @lines; }
	
	$head  =  shift @lines;
	if ($first_line) {
	    $first_line = 0;
	    if (not ($head  =~ s/^\s*>//)) {
		$trouble = 1;
		warn $head;
		die "ERROR: File does not appear to be in FASTA format\n";
	    }
	}
	else {
	    if ($head  =~ s/^\s*>//) {
		$trouble = 1;
		warn $head;
		die "Spurious beginning-of record mark found in record $.\n";
	    }
	}
	
	foreach my $ln (@lines) {
	    $_  =  $ln;
	    $ln =~ s/\s//g;
	    
	    print STDERR "$head: contains X's\n"    if ($ln =~ s/x/n/ig);
	    print STDERR "$head: contains colons\n" if ($ln =~ s/://g);
	    
	    while ($ln =~ s/([^ACGTUMRWSYKBDHVN]+)/n/i) {
		$trouble = 1;
		$badchars++;
		print STDERR ">$head:\tbad char $1 at ", pos($ln), " at line $.\n";
	    }
	}
	
	$seq   =  join( "", @lines );
	$seq   =~ s/\cM//g;
	$seq   =~ tr/a-z/A-Z/;
	@record = ($head, \$seq);
    }
    
    $/ = $old_eol;
    return @record;
}


sub file_type {
#...No longer used; Save for archival purposes...
    my ($file_name) = @_;
    
    my $file_type;
    if (($file_type = `file '$file_name' | cut -f2 -d:`) && ($file_type =~ m/\S/o)) {
	my $saved_file_type = $file_type;
	
	$file_type =~ s/^\s+//o;   #...trim leading whitespace
	$file_type =~ s/\s+$//o;   #...trim trailing whitespace
	$file_type =~ s/, with very long lines//;
	
	print STDERR "file_type = $file_type\n" if $ENV{VERBOSE};
	
	if    ($file_type =~ m/^ASCII.*text$/) {
	    print STDERR "ASCII text file\n" if $ENV{VERBOSE};
	}
	elsif ($file_type =~ m/^ASCII.*text, with CR line terminators$/) {
	    print STDERR "CR terminated file\n" if $ENV{VERBOSE};
	    $input_eol_marker = "\cM";
	}
	elsif ($file_type =~ m/^ASCII.*text, with CRLF line terminators$/) {
	    print STDERR "CRLF terminated file\n" if $ENV{VERBOSE};
	    $input_eol_marker = "\cM\cJ";
	}
	elsif ($file_type =~ m/^ASCII.*text, with CR, LF line terminators$/) {
	    print STDERR "CR, LF terminated file\n" if $ENV{VERBOSE};
	    $input_eol_marker = "\cM\cJ";
	}
	elsif ($file_type =~ m/^ASCII.*text, with CRLF, LF line terminators$/) {
	    print STDERR "CRLF, LF terminated file\n" if $ENV{VERBOSE};
	    $input_eol_marker = "\cM\cJ\n";
	}
	else {
	    die "Could not handle file-type $saved_file_type";
	}
    }
}

__DATA__
