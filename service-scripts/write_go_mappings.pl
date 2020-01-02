#
# Write the GO mappings for the data in the given
# genome directory. (Using FIGV to access, originally
# meant for use in RAST).
#

use File::Basename;
use FIG_Config;
use FIGV;
use strict;
use Getopt::Long;

my $pegs;
my $ec_file;
my $peg_fh;
my $override_file;
my $assigns_file;
my $genomes;
my $go_file;

my $rc = GetOptions("pegs" => \$pegs,
		    "genomes" => \$genomes,
		    "override=s" => \$override_file,
		    "assigns=s" => \$assigns_file,
		    "ec=s" => \$ec_file,
		    "go=s" => \$go_file);

($rc && @ARGV == 1) or die "Usage: $0 [-pegs] [-ec ec-file] genome-dir > go.mappings";

my $dir = shift;

my @genomes;

if ($genomes)
{
    my $fh;
    if ($dir eq '-')
    {
	$fh = \*STDIN;
    }
    else
    {
	open($fh, "<", $dir) or die "Cannot open $dir: $!";
    }
    while (<$fh>)
    {
	if (/(\d+\.\d+)/)
	{
	    push(@genomes, $1);
	}
    }
    close($fh);
}
elsif ($pegs)
{
    if ($dir eq '-')
    {
	$peg_fh = \*STDIN;
    }
    else
    {
	open($peg_fh, "<", $dir) or die "Cannot open $dir: $!";
    }
}
else
{
    # genome-dir option
    if (! -d $dir)
    {
	die "Genome directory $dir does not exist";
    }
}

my %override;
if ($override_file)
{
    if (open(FH, "<", $override_file))
    {
	while (<FH>)
	{
	    chomp;
	    my($peg, $fn) = split(/\t/);
	    $override{$peg} = $fn;
	}
	close(FH);
    }
    else
    {
	die "Cannot open override file $override_file: $!";
    }
}

my $go_fh;
if (defined($go_file))
{
    open($go_fh, ">", $go_file) or die "Cannot write $go_file: $!\n";
}
else
{
    $go_fh = \*STDOUT;
}

my $assigns_fh;
if (defined($assigns_file))
{
    open($assigns_fh, ">", $assigns_file) or die "Cannot write $assigns_file: $!\n";
}
    
my $ec_fh;
if (defined($ec_file))
{
    open($ec_fh, ">", $ec_file) or die "Cannot write $ec_file: $!\n";
}

my $fig = -d $dir ? new FIGV($dir) : new FIG();

# get the functional role name -> GO file
open(FH, $FIG_Config::data . "/Ontologies/GO/fr2go") or die "could not open fr2go";
my $fr2go = {};
while (<FH>) {
    chomp;
    my ($fr, $go) = split /\t/;
    $fr2go->{$fr} = [] unless (exists $fr2go->{$fr});
    push @{$fr2go->{$fr}}, $go;
}
close FH;

if ($genomes)
{
    for my $g (@genomes)
    {
	print STDERR "$g\n";
	for my $ent (@{$fig->all_features_detailed_fast($g)})
	{
	    my $fid = $ent->[0];
	    my $func = $override{$fid} || $ent->[6];
	    process($fid, $func);
	}
    }
}
elsif (!$pegs)
{

    my $genome = basename($dir);

    # get the pegs
    foreach my $peg (sort { &FIG::by_fig_id($a,$b) } $fig->pegs_of($genome), $fig->rnas_of($genome))
    {
	my $func = $override{$peg} || $fig->function_of($peg);

	process($peg, $func);
    }
}
else
{
    while (<$peg_fh>)
    {
	chomp;
	my($peg, $func) = split(/\t/);

	process($peg, $func);
    }
}

close($ec_fh) if defined($ec_fh);
	    
sub process
{
    my($peg, $func) = @_;
    
    my %ecs;
    my @gos;

    print $assigns_fh "$peg\t$func\n" if $assigns_fh;
    
    # get EC / GO from role
    if (defined $func) {
	foreach my $role ($fig->roles_of_function($func)) {
	    my ($ec) = ($role =~ /\(EC ((\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-))\)/);
	    $ecs{$ec} = 1 if ($ec);
	    push @gos, @{$fr2go->{$role}} if ($fr2go->{$role});
	}
    }

    my @ecs = keys %ecs;
    if (@ecs && defined($ec_fh))
    {
	print $ec_fh join("\t", $peg, $func, @ecs), "\n";
    }

    return unless @gos;
    
    my %gos = map { $_ => 1 } @gos;
    print $go_fh join("\t", $peg, $func, sort keys %gos), "\n";
}

