#
# Create a btree containing the sequences in a fasta file.
#
#

my $usage = "make_fasta_btree fasta-file btree-file len-btree-file";

use strict;
use FIG;
use DB_File;

@ARGV == 3 or die $usage;

my $fasta_file = shift;
my $btree_file = shift;
my $len_btree_file = shift;

open(FA, "<$fasta_file") or die "Cannot open $fasta_file: $!";
#
# Unlink existing file to ensure it is clean.
#
if (-f $btree_file)
{
    unlink($btree_file);
}
if (-f $len_btree_file)
{
    unlink($len_btree_file);
}

my %btree;
my $btree_tie = tie %btree, 'DB_File', $btree_file, O_RDWR | O_CREAT, 0666, $DB_BTREE;

$btree_tie or die "Cannot create btree $btree_file: $!\n";

my %len_btree;
my $len_btree_tie = tie %len_btree, 'DB_File', $len_btree_file, O_RDWR | O_CREAT, 0666, $DB_BTREE;

$len_btree_tie or die "Cannot create btree $len_btree_file: $!\n";

my($seqp, $id, $comment);

while (($id, $seqp, $comment) = &FIG::read_fasta_record(\*FA))
{
    $btree{$id} = $$seqp;
    $len_btree{$id} = length($$seqp);
}
untie %btree;
