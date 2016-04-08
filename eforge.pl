#!/usr/bin/perl

=head1 NAME

eforge.pl - Experimentally derived Functional element Overlap analysis of ReGions from EWAS.

=head1 SYNOPSIS

eforge.pl options (-f file) (-mvp mvplist)

=head1 DESCRIPTION

Analyse a set of MVPs for their overlap with DNase 1 hotspots compared to matched background MVPs. 
Identifies enrichment in DHS by tissue and plots graphs and table to display. Arbitrarily a minumum of 5* MVPs is required.  
Note that if no MVPs are given the script will run on A DEFAULT EWAS* as an example output.

Several outputs are made.

A straight base R graphics pdf chart of the data.

A polychart (https://github.com/Polychart/polychart2) interactive javascript graphic using rCharts (http://ramnathv.github.io/rCharts/).

A dimple (http://dimplejs.org) d3 interactive graphic using rCharts.

A table using the Datatables (https://datatables.net) plug-in for the jQuery Javascript library, again accessed through rCharts.

In each of the graphics the colouring should be consistent. Blue (p value > 0.05), light red or pink (0.05 => p value > 0.01), red or dark red (p value <= 0.01 ) for the 95% and 99% cIs. 
Or whatever other thresholds are specified. 

eForge functions, plotting options and stats are provided by eForge::eForge, eForge::ePlot and eForge::eStats modules.

=head1 OPTIONS

=over

=item B<data>

Dataset to analyse. Either ENCODE data ('encode'), unconsolidated Roadmap Epigenome data ('erc'), consolidated Roadmap Epigenome data ('erc2'), or Blueprint data ('blueprint'). erc by default.

Use -data ? to get a list of available datasets on your local install.

=item B<peaks>

Use peaks instead of hotspots. Peaks are more stringent DNase1 peaks calls representing DNase hypersensitive sites, 
rather than hotspots which are regions of generalised DNase1 sensitivity or open chromatin. Default is to use hotspots.

=item B<bkgd>

Background is set at default to 450k array ('450k'), the Illumina Infinium HumanMethylation450 BeadChip.

For the time being, it is suficient for MVPs to be on the 450k array. Probes within 1kb of each other will undergo filtering.

Use -bkgd ? to get a list of available backgrounds on your local install.

=item B<label>

Supply a label that you want to use for the plotting titles, and filenames.

=item B<f>

Supply the name of a file containing a list of MVPs. 
Format must be given by the -format flag. 
If not supplied the analysis is performed either on mvps provided as probeids (cg or ch probes) in a comma separated list through the mvps option 
or on a set of data from a default ewas* study (provide link for study here*). Note that 5* MVPs are required at a minimum.

=item B<mvps>

Can provide the mvps as probeids in a comma separated list.

=item B<min_mvps>

Specify the minimum number of MVPs to be allowed. Default is 5 now we are using binomial test.

=item B<thresh>

Alter the default binomial p value thresholds. Give a comma separate list of three e.g. 0.05,0.01 for the defaults

=item B<format>

If f is specified, specify the file format as follow:

probeid = list of mvps as probeids each on a separate line. Optionally can add other fields after the probeid which are ignored, 
unless the pvalue filter is specified, in which case eForge assumes that the second field is the minus log10 pvalue

bed  = File given is a bed file of locations (chr\tbeg\tend).  bed format should be 0 based and the chromosome should be given as chrN. 
However we will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same*.

tabix = File contains MVPs in tabix format.

ian = 1-based chr\tbeg\tend\tprobeid\tpval\tminuslog10pval

=item B<filter>

Set a filter on the MVPs based on the -log10 pvalue.  This works for files in the 'ian' or 'probeid' format. 
Give a value as the lower threshold and only MVPs with -log10 pvalues >= to the threshold will be analysed. Default is no filtering.

=item B<bkgrd>

Output background stats for investigation.

=item B<reps>

The number of background matching sets to pick and analyse. Default 1000*.

=item B<proxy>

Apply filter for MVPs in proximity (within 1 kb of another test MVP). With proximity filter specified, eForge will report MVPs removed due to proximity with another MVP in the list and will randomly pick one of the probes among the set of probes that are in proximity (within 1 kb of each other).
To turn off proximity filtering specify -noproxy

=item B<noproxy>

Turn off proximity filtering.

=item B<depletion>

Analyse for DHS depletion pattern instead of the default DHS enrichment analysis. Use when dealing with datasets suspected not to overlap with DHS. Specifying depletion will be indicated on the label (the text "Depletion Analysis" will be added to the file label).

=item B<noplot>

Just make the data file, don't plot.

=item B<help|h|?>

Print a brief help message and exits.

=item B<man|m>

Print this perldoc and exit.

=back

=head1 LICENCE AND COPYRIGHT

eforge.pl Functional analysis of EWAS MVPs

Copyright (C) [2014-2015] EMBL - European Bioinformatics Institute and University College London

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 dated June, 1991 or at your option
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the source tree;
if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

=head1 CONTACT

Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

Javier Herrero, C<< <javier.herrero at ucl.ac.uk> >>

=head1 ACKNOWLEDGEMENTS

This software is based on the FORGE tool developed by Ian Dunham at the EMBL-EBI

Javier Herrero <javier.herrero@ucl.ac.uk>

=cut

use strict;
use 5.010;
use warnings;
use DBI; #database link to sqlite database
use Sort::Naturally;
use Cwd;
use Getopt::Long; #check this module
use File::Basename;
use Config::IniFiles;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use eForge::eStats;
use eForge::ePlot;
use eForge::eForge;
use Data::UUID;
use Statistics::Multtest qw(BY);


my $cwd = getcwd;

my $dbname = "eforge_1.2.db";

my $bkgd; # Default value
my $bkgd_label;
my ($data, $peaks, $label, $file, $format, $min_mvps, $bkgrdstat, $noplot, $reps,
    $help, $man, $thresh, $proxy, $noproxy, $depletion, $filter, $out_dir, @mvplist,
    $web, $autoopen);

GetOptions (
    'data=s'     => \$data,
    'peaks'      => \$peaks,
    'bkgrd'      => \$bkgrdstat,
    'bkgd=s'     => \$bkgd,
    'label=s'    => \$label,
    'f=s'        => \$file,
    'format=s'   => \$format,
    'mvps=s'     => \@mvplist,
    'min_mvps=i' => \$min_mvps,
    'noplot'     => \$noplot,
    'reps=i'     => \$reps,
    'thresh=s'   => \$thresh,
    'proxy=s'    => \$proxy,
    'noproxy'    => \$noproxy,
    'depletion'  => \$depletion,
    'filter=f'   => \$filter,
    'out_dir=s'  => \$out_dir,
    'web=s'        => \$web,
    'autoopen'   => \$autoopen,
    'help|h|?'   => \$help,
    'man|m'      => \$man,

);


pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

if (!$out_dir) {
    my $ug = new Data::UUID;
    $out_dir = $ug->to_hexstring($ug->create());
}

# the minimum number of mvps allowed for test. Set to 5 as we have binomial p
unless (defined $min_mvps) {
    $min_mvps = 5;
}

# Label for plots
unless (defined $label) {
    $label = "No label given";
}
  
## ============================================================================
## Connect to the DB
## ============================================================================
# This reads the config file and sets up the $datadir variable
my $dirname = dirname(__FILE__);
my $cfg = Config::IniFiles->new( -file => "$dirname/eforge.ini" );
my $datadir = $cfg->val('Files', 'datadir');

unless (-s "$datadir/$dbname") {
    die "Database $dbname not found or empty";
}
my $dsn = "dbi:SQLite:dbname=$datadir/$dbname";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;
## ============================================================================


## ============================================================================
## Check the dataset against the info on the DB
## ============================================================================
my $all_datasets = get_all_datasets($dbh);
if (!defined($all_datasets)) {
    die "Empty database: no dataset loaded!\n";
} elsif (!defined($data)) {
    $data = $all_datasets->[0]->{tag};
    print "Using default dataset: [$data] ".$all_datasets->[0]->{name}."\n";
} elsif ($data eq "?") {
    print "Available datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
    exit();
} elsif (!grep {$_ eq $data} map {$_->{tag}} @$all_datasets) {
    die "Dataset $data unknown\nAvailable datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
}
## ============================================================================


## ============================================================================
## Check the array name (A.K.A. background) against DB
## ============================================================================
my $all_arrays = get_all_arrays($dbh);
if (!defined($all_arrays)) {
    die "Empty database: no background loaded!\n";
} elsif (!defined($bkgd)) {
    $bkgd = $all_arrays->[0]->{tag};
    print "Using default background: [$bkgd] ".$all_arrays->[0]->{name}."\n";
    $bkgd_label = $all_arrays->[0]->{name};
} elsif ($bkgd eq "?") {
    print "Available arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
    exit();
} elsif (!grep {$_ eq $bkgd} map {$_->{tag}} @$all_arrays) {
    die "Array $bkgd unknown\nAvailable arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
} else {
    foreach my $this_array (@$all_arrays) {
        if ($this_array->{tag} eq $bkgd) {
            $bkgd_label = $this_array->{name};
            last;
        }
    }
}
## ============================================================================


## ============================================================================
## Check the proxy_filter (A.K.A. filter) against DB
## ============================================================================
# Set proximity filter
unless (defined $noproxy) {
    my $all_proxy_filters = get_all_proxy_filters($dbh);
    if ($all_proxy_filters->{$bkgd}) {
        $proxy = $all_proxy_filters->{$bkgd};
    }
}
## ============================================================================


if (defined $depletion) {
    $label = "$label.depletion";
}

#regexp puts underscores where labels before
(my $lab = $label) =~ s/\s/_/g;
$lab = "$lab.$bkgd.$data";

#format for reading from file
unless (defined $format) {
    $format = 'probeid';
}





# percentile bins for the bkgrd calculations. 
#This is hard coded so there are enough probes to choose from, but could later be altered.
#eForge, unlike Forge, does not use percentile bins
my $per = 10;


# number of sets to analyse for bkgrd.
unless (defined $reps) {
    $reps = 1000;
}


# Define the thresholds to use.
my ($t_marginal, $t_strict);
if (defined $thresh) {
    ($t_marginal, $t_strict) = split(",", $thresh);
    unless (looks_like_number($t_marginal) && looks_like_number($t_strict)){
        die "You must specify numerical p value thresholds in a comma separated list\n";
    }
    unless ((1 >= $t_marginal) && ($t_marginal >= $t_strict) && ($t_strict >= 0)) {
        die "The p-value thresholds should be 1 >= T.marginal >= T.strict >= 0\n";
    }
} else {
    $t_marginal = 0.05; # set binomial p values, bonferroni is applied later based on number of samples (cells)
    $t_strict = 0.01;
}

# mvps need to come either from a file or a list
my $mvps;


# A series of data file formats to accept.

warn "[".scalar(localtime())."] Processing input...\n";

if (defined $file) {
    if (defined $filter) {
        unless ($format eq "ian" or $format eq "probeid") {
            warn "You have specified p value filtering, but this isn't implemented for files of format $format. No filtering will happen."
	    }
    }
    open my $fh, "<", $file or die "cannot open file $file : $!";
    $mvps = process_file($fh, $format, $dbh, $bkgd, $filter);

} elsif (@mvplist) {
    @$mvps = split(/,/,join(',',@mvplist));

} else{
    # Test MVPs from Liu Y et al. Nat Biotechnol 2013  Pulmonary_function.snps.bed (*put EWAS bedfile here)
    # If no options are given it will run on the default set of MVPs
    warn "No probe input given, so running on default set of probes, a set of monocyte tDMPs from Jaffe AE and Irizarry RA, Genome Biol 2014.";
    @$mvps = qw(cg13430807 cg10480329 cg06297318 cg19301114 cg23244761 cg26872907 cg18066690 cg04468741 cg16636767 cg10624395 cg20918393);
}


# Remove redundancy in the input

my %nonredundant;
foreach my $mvp (@$mvps) {
    $nonredundant{$mvp}++;
}


foreach my $mvp (keys %nonredundant) {
    if ($nonredundant{$mvp} > 1) {
        say "$mvp is present " . $nonredundant{$mvp} . " times in the input. Analysing only once."
    }
}

@$mvps = keys %nonredundant;
my @origmvps = @$mvps;


#######!!!!!### proximity filtering starts below:

my ($prox_excluded, $output, $input);
unless(defined $noproxy) {
    $input = scalar @$mvps;
    ($prox_excluded, @$mvps) = prox_filter($mvps, $dbh);
    while (my ($excluded_mvp, $mvp) = each %$prox_excluded) {
        warn "$excluded_mvp excluded for proximity (1 kb) with $mvp\n";
    }

    $output = scalar @$mvps;
}

# Check we have enough MVPs
if (scalar @$mvps < $min_mvps) {
    die "Fewer than $min_mvps MVPs. Analysis not run\n";
}


# get the cell list array and the hash that connects the cells and tissues
my ($cells, $tissues) = get_cells($data, $dbh);

# get the bit strings for the test mvps from the database file
my $rows = get_bits($dbh, $data, $bkgd, $mvps);

# unpack the bitstrings and store the overlaps by cell.
my $test = process_bits($rows, $cells, $data);

# generate stats on the background selection
if (defined $bkgrdstat) {
    bkgrdstat($test, $lab, "test");
}



# Identify probeids that weren't found and warn about them.
#the reference to hash key 'MVPS' is because of use of perl module eforge.pm from the eForge tool
#(in subroutines process_bits etc)

my @missing;
foreach my $probeid (@origmvps) {
    if (defined $proxy) {
        next if exists $$prox_excluded{$probeid};
    }
    unless (exists $$test{'MVPS'}{$probeid}) {
        push @missing, $probeid;
    }
}

if (scalar @missing > 0) {
    warn "The following " . scalar @missing . " MVPs have not been analysed because they were not found on the ".$bkgd_label."\n";
    warn join("\n", @missing) . "\n";
}
if (defined $proxy) {
    if ($output < $input) {
        warn "For $label, $input MVPs provided, " . scalar @$mvps . " retained, " . scalar @missing . " not analysed, "  . scalar(keys %$prox_excluded) . " proximity filtered at 1 kb\n";
    }
}

# only pick background mvps matching mvps that had bitstrings originally.
#reference to hash key 'MVPS' is due to use of eforge.pm module from eForge tool
#(in subroutines process_bits, etc)

my @foundmvps = keys %{$$test{'MVPS'}};
my $mvpcount = scalar @foundmvps;


# identify the feature and cpg island relationship, and then make bkgrd picks (old version just below)
#my $picks = match(\%$test, $bkgd, $datadir, $per, $reps);
warn "[".scalar(localtime())."] Loading the $bkgd background...\n";
my $picks = match(\%$test, $bkgd, $datadir, $reps);

########check below lines:
 
# for bkgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
my %bkgrd; #this hash is going to store the bkgrd overlaps

# Get the bits for the background sets and process
my $backmvps;

warn "[".scalar(localtime())."] Running the analysis with $mvpcount MVPs...\n";
my $num = 0;
foreach my $bkgrd (keys %{$picks}) {
    warn "[".scalar(localtime())."] Repetition $num out of ".$reps."\n" if (++$num%100 == 0);
    #$rows = get_bits(\@{$$picks{$bkgrd}}, $sth);
    $rows = get_bits($dbh, $data, $bkgd, \@{$$picks{$bkgrd}});
    $backmvps += scalar @$rows; #$backmvps is the total number of background probes analysed
    unless (scalar @$rows == scalar @foundmvps) {
        warn "Background " . $bkgrd . " only " . scalar @$rows . " probes out of " . scalar @foundmvps . "\n";
    }
    my $result = process_bits($rows, $cells, $data);
    foreach my $cell (keys %{$$result{'CELLS'}}) {
        push @{$bkgrd{$cell}}, $$result{'CELLS'}{$cell}{'COUNT'}; # accumulate the overlap counts by cell
    }
    if (defined $bkgrdstat) {
        bkgrdstat($result, $lab, $bkgrd);
    }
}

$dbh->disconnect();
warn "[".scalar(localtime())."] All repetitions done.\n";

warn "[".scalar(localtime())."] Calculating p-values...\n";
#Having got the test overlaps and the bkgd overlaps now calculate p values and output 
#the table to be read into R for plotting.


mkdir $out_dir;

open my $bfh, ">", "$out_dir/background.tsv" or die "Cannot open background.tsv";



my @results;
my @pvalues;
###ncmp is a function from Sort::Naturally
foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells){
    # above line sorts by the tissues alphabetically (from $tissues hash values)

    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "MVPs")
    say $bfh join("\t", @{$bkgrd{$cell}});
    my $teststat = ($$test{'CELLS'}{$cell}{'COUNT'} or 0); #number of overlaps for the test MVPs

    # binomial pvalue, probability of success is derived from the background overlaps over the tests for this cell
    # $backmvps is the total number of background mvps analysed
    # $tests is the number of overlaps found over all the background tests
    my $tests;
    foreach (@{$bkgrd{$cell}}) {
        $tests += $_;
    }
    my $p = sprintf("%.6f", $tests/$backmvps);

    # binomial probability for $teststat or more hits out of $mvpcount mvps
    # sum the binomial for each k out of n above $teststat
    my $pbinom;
    if (defined $depletion) {
        foreach my $k (0 .. $teststat) {
            $pbinom += binomial($k, $mvpcount, $p);
        }
    } else {
        foreach my $k ($teststat .. $mvpcount) {
            $pbinom += binomial($k, $mvpcount, $p);
        }
    }
    if ($pbinom >1) {
        $pbinom=1;
    }
    # Store the p-values in natural scale (i.e. before log transformation) for FDR correction
    push(@pvalues, $pbinom);
    $pbinom = sprintf("%.2e", $pbinom);

    # Z score calculation (note: this is here only for legacy reasons. Z-scores assume normal distribution)
    my $zscore = zscore($teststat, $bkgrd{$cell});

    my $mvp_string = "";
    $mvp_string = join(",", @{$$test{'CELLS'}{$cell}{'MVPS'}}) if defined $$test{'CELLS'}{$cell}{'MVPS'};
    # This gives the list of overlapping MVPs for use in the tooltips. If there are a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.

    push(@results, [$zscore, $pbinom, $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'datatype'}, $$tissues{$cell}{'file'}, $mvp_string, $$tissues{$cell}{'acc'}]);
}
close($bfh);

# Correct the p-values for multiple testing using the Benjamini-Yekutieli FDR control method
my $qvalues = BY(\@pvalues);
$qvalues = [map {sprintf("%.2e", $_)} @$qvalues];

# Write the results to a tab-separated file
my $filename = "$lab.chart.tsv";
open my $ofh, ">", "$out_dir/$filename" or die "Cannot open $out_dir/$filename: $!";
print $ofh join("\t", "Zscore", "Pvalue", "Cell", "Tissue", "Datatype", "File", "Probe", "Accession", "Qvalue"), "\n";
for (my $i = 0; $i < @results; $i++) {
    print $ofh join("\t", @{$results[$i]}, $qvalues->[$i]), "\n";
}
close($ofh);


warn "[".scalar(localtime())."] Generating plots...\n";
unless (defined $noplot){
    #Plotting and table routines
    Chart($filename, $lab, $out_dir, $tissues, $cells, $label, $t_marginal, $t_strict, $data); # basic pdf plot
    dChart($filename, $lab, $out_dir, $data, $label, $t_marginal, $t_strict, $web); # rCharts Dimple chart
    table($filename, $lab, $out_dir, $web); # Datatables chart
  }

warn "[".scalar(localtime())."] Done.\n";

if ($autoopen) {
    system("open $out_dir/$lab.table.html");
    system("open $out_dir/$lab.dchart.html");
    system("open $out_dir/$lab.chart.pdf");
}
