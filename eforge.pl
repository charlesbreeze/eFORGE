#!/usr/bin/env perl

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

=item B<--dataset TAG>

Set of functional data to look for enrichment. Either ENCODE data ('encode'), unconsolidated Roadmap
Epigenome data ('erc'), consolidated Roadmap Epigenome data ('erc2'), or Blueprint data ('blueprint').
erc by default.

Use --dataset ? to get a list of available datasets on your local install.

=item B<--array TAG>

Array (FKA background) is set at default to 450k array ('450k'), the Illumina Infinium HumanMethylation450 BeadChip.

For the time being, it is suficient for MVPs to be on the 450k array. Probes within 1kb of each other
will undergo filtering.

Use --array ? to get a list of available backgrounds on your local install.

=item B<--label STRING>

Supply a label that you want to use for the plotting titles, and filenames.

=item B<--f FILENAME>

Supply the name of a file containing a list of MVPs. 
Format must be given by the -format flag. 
If not supplied the analysis is performed either on mvps provided as probeids (cg or ch probes) in a
comma separated list through the mvps option or on a set of data from a default ewas study, namely a
set of monocyte tDMPs from Jaffe AE and Irizarry RA, Genome Biol 2014.

Note that at least 5 MVPs are required at a minimum by default.

=item B<--mvps probe_id,probe_id...>

Can provide the mvps as probeids in a comma separated list.

=item B<--min_mvps INT>

Specify the minimum number of MVPs to be allowed. Default is 5 now we are using binomial test.

=item B<--thresh FLOAT,FLOAT>

Alter the default binomial p value thresholds. Give a comma separate list of three e.g. 0.05,0.01 for the defaults

=item B<--format STRING>

If f is specified, specify the file format as follow:

probeid = list of mvps as probeids each on a separate line. Optionally can add other fields after the probeid which are ignored,
unless the pvalue filter is specified, in which case eForge assumes that the second field is the minus log10 pvalue

bed  = File given is a bed file of locations (chr\tbeg\tend).  bed format should be 0 based and the chromosome should be given as chrN.
However we will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same*.

tabix = File contains MVPs in tabix format.

=item B<--filter FLOAT>

Set a filter on the MVPs based on the -log10 pvalue.  This works for files in the probeid' format.
Give a value as the lower threshold and only MVPs with -log10 pvalues >= to the threshold will be
analysed. Default is no filtering.

=item B<--bkgrd>

Output background stats for investigation.

=item B<--reps INT>

The number of background matching sets to pick and analyse. Default 1000.

=item B<--proxy TAG>

Apply filter for MVPs in proximity (within 1 kb of another test MVP). With proximity filter specified,
eForge will report MVPs removed due to proximity with another MVP in the list and will randomly pick
one of the probes among the set of probes that are in proximity (within 1 kb of each other).

To turn off proximity filtering specify -noproxy

=item B<--noproxy>

Turn off proximity filtering.

=item B<--depletion>

Analyse for depletion pattern instead of the default enrichment analysis. Use when dealing with
datasets suspected not to overlap with DHS (or the relevant functional assay). Specifying depletion
will be indicated on the label (the text "Depletion Analysis" will be added to the file label).

=item B<--noplot>

Just make the data file, don't plot.

=item B<--help|-h|-?>

Print a brief help message and exits.

=item B<--man|-m>

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

my $array; # Default value
my $array_label;
my $format = 'probeid'; # Input format
my $label = 'Unnamed'; # Label for plots
my $reps = 1000;
# set binomial p values, multiple test correction is used
my $thresh; # string for command line option
my $t_marginal = 0.05; # default marginal p-value threshold
my $t_strict = 0.01; # default strict p-value threshold

my $min_num_probes = 5; # the minimum number of probes allowed for test. Set to 5 as we have binomial p

my ($dataset, $filename, $bkgrdstat, $noplot,
    $help, $man, $proxy, $noproxy, $depletion, $filter, $out_dir, $probe_list,
    $web, $autoopen);

GetOptions (
    'dataset=s'  => \$dataset,
    'bkgrd'      => \$bkgrdstat,
    'array|bkgd=s' => \$array,
    'label=s'    => \$label,
    'f=s'        => \$filename,
    'format=s'   => \$format,
    'probes|mvps=s@' => \$probe_list,
    'min_num_probes|min_mvps=i' => \$min_num_probes,
    'noplot'     => \$noplot,
    'reps=i'     => \$reps,
    'thresh=s'   => \$thresh,
    'proxy=s'    => \$proxy,
    'noproxy'    => \$noproxy,
    'depletion'  => \$depletion,
    'filter=f'   => \$filter,
    'out_dir=s'  => \$out_dir,
    'web=s'      => \$web,
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

# Define the thresholds to use.
if ($thresh) {
    ($t_marginal, $t_strict) = parse_pvalue_thresholds($thresh);
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
} elsif (!defined($dataset)) {
    $dataset = $all_datasets->[0]->{tag};
    print "Using default dataset: [$dataset] ".$all_datasets->[0]->{name}."\n";
} elsif ($dataset eq "?") {
    print "Available datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
    exit();
} elsif (!grep {$_ eq $dataset} map {$_->{tag}} @$all_datasets) {
    die "Dataset $dataset unknown\nAvailable datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
}
## ============================================================================


## ============================================================================
## Check the array name (A.K.A. background) against DB
## ============================================================================
my $all_arrays = get_all_arrays($dbh);
if (!defined($all_arrays)) {
    die "Empty database: no background loaded!\n";
} elsif (!defined($array)) {
    $array = $all_arrays->[0]->{tag};
    print "Using default background: [$array] ".$all_arrays->[0]->{name}."\n";
    $array_label = $all_arrays->[0]->{name};
} elsif ($array eq "?") {
    print "Available arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
    exit();
} elsif (!grep {$_ eq $array} map {$_->{tag}} @$all_arrays) {
    die "Array $array unknown\nAvailable arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
} else {
    foreach my $this_array (@$all_arrays) {
        if ($this_array->{tag} eq $array) {
            $array_label = $this_array->{name};
            last;
        }
    }
}
## ============================================================================


## ============================================================================
## Check the proxy_filter (A.K.A. filter) against DB
## ============================================================================
# Set proximity filter
if (defined $noproxy) {
    $proxy = undef;
} else {
    my $all_proxy_filters = get_all_proximity_filters($dbh);
    if ($all_proxy_filters->{$array}) {
        $proxy = $all_proxy_filters->{$array};
    }
}
## ============================================================================


## ============================================================================
## Append main options (depletion on/off; array; dataset) to $label
## ============================================================================
if (defined $depletion) {
    $label = "$label.depletion";
}
(my $lab = $label) =~ s/\s/_/g; # Avoid whitespaces on the label
$lab = "$lab.$array.$dataset";
## ============================================================================


## ============================================================================
## Read and process the input MVPs
## ============================================================================
warn "[".scalar(localtime())."] Processing input...\n";
# This will read the probes from the file if provided, from the probe list otherwise or use the
# example data set as a last resort.
my $mvps = get_input_probes($filename, $probe_list);
my $original_mvps = [@$mvps];
my $num_of_input_mvps = scalar(@$mvps);

# Apply the proximity filter if requested
my ($proximity_excluded);
if(defined $proxy) {
    ($proximity_excluded, $mvps) = proximity_filter($dbh, $array, $mvps);
    while (my ($excluded_mvp, $mvp) = each %$proximity_excluded) {
        warn "$excluded_mvp excluded for $proxy proximity filter with $mvp\n";
    }
}

# $annotated_probes is an arrayref with probe_id, sum, bit, gene_group, cgi_group for each input probe
my $annotated_probes = get_probe_annotations_and_overlap_for_dataset($dbh, $dataset, $array, $mvps);
my $existing_probes = {map {$_->[0] => 1} @$annotated_probes};
$mvps = [keys %$existing_probes];

## Detect and remove the missing probes.
my $num_missing_probes = find_missing_probes($original_mvps, $existing_probes);

# Print summary of filtering and checks:
my $msg = "For $label, $num_of_input_mvps MVPs provided, ". scalar @$mvps.
        " retained: $num_missing_probes were not found";
if (defined $proxy) {
    $msg .= " and " . scalar(keys %$proximity_excluded) . " excluded using $proxy proximity filter";
}
warn $msg, ".\n";

# Check we have enough MVPs left
my $num_of_valid_probes = scalar @$mvps;
if ($num_of_valid_probes < $min_num_probes) {
    die "Fewer than $min_num_probes MVPs. Analysis not run\n";
}
## ============================================================================


# get the cell list array and the hash that connects the cells and tissues
# $samples is a hash whose keys are the $cells (short name for the cell type/lines) and value is
# another hash with 'tissue', 'datatype', 'file' and 'acc' keys.
# IMPORTANT: $cells contains the list of cells in the order defined in the DB. This is critical
# to correctly assign each bit to the right sample.
my ($cells, $samples) = get_samples_from_dataset($dbh, $dataset);

# unpack the bitstrings and store the overlaps by cell.
# $test is a complex hash like:
# $test{'MVPS'}{$probe_id}{'SUM'} (total number of overlaps of this probe with DHS/etc in this dataset)
# $test{'MVPS'}{$probe_id}{'PARAMS'} (gene and CGI annotations for this probe)
# $test{'CELLS'}{$cell}{'COUNT'} (number of input MVPs that overlap with the signal on this sample)
# $test{'CELLS'}{$cell}{'MVPS'} (list of input MVPs that overlap with the signal on this sample)
my $test = process_overlaps($annotated_probes, $cells, $dataset);

# generate stats on the background selection
if (defined $bkgrdstat) {
    bkgrdstat($test, $lab, "test");
}



# only pick background mvps matching mvps that had bitstrings originally.
#reference to hash key 'MVPS' is due to use of eforge.pm module from eForge tool
#(in subroutines process_overlaps, etc)


# identify the feature and cpg island relationship, and then make bkgrd picks (old version just below)
#my $picks = match(\%$test, $array, $datadir, $per, $reps);
warn "[".scalar(localtime())."] Loading the $array background...\n";
my $picks = match(\%$test, $array, $datadir, $reps);

########check below lines:
 
# for bkgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
my %bkgrd; #this hash is going to store the bkgrd overlaps

# Get the bits for the background sets and process
my $backmvps;

warn "[".scalar(localtime())."] Running the analysis with $num_of_valid_probes MVPs...\n";
my $num = 0;
foreach my $bkgrd (keys %{$picks}) {
    warn "[".scalar(localtime())."] Repetition $num out of ".$reps."\n" if (++$num%100 == 0);
    #$annotated_probes = get_probe_annotations_and_overlap_for_dataset(\@{$$picks{$bkgrd}}, $sth);
    $annotated_probes = get_probe_annotations_and_overlap_for_dataset($dbh, $dataset, $array, \@{$$picks{$bkgrd}});
    $backmvps += scalar @$annotated_probes; #$backmvps is the total number of background probes analysed
    unless (scalar @$annotated_probes == $num_of_valid_probes) {
        warn "Background " . $bkgrd . " only " . scalar @$annotated_probes . " probes out of $num_of_valid_probes\n";
    }
    my $result = process_overlaps($annotated_probes, $cells, $dataset);
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

if (!$web) {
    open(BACKGROUND, "| gzip -9 > $out_dir/background.tsv.gz") or die "Cannot open background.tsv";
}



my @results;
my @pvalues;
###ncmp is a function from Sort::Naturally
foreach my $cell (sort {ncmp($$samples{$a}{'tissue'},$$samples{$b}{'tissue'}) || ncmp($a,$b)} @$cells){
    # above line sorts by the tissues alphabetically (from $samples hash values)

    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "MVPs")
    if (!$web) {
        print BACKGROUND join("\t", @{$bkgrd{$cell}}), "\n";
    }
    my $teststat = ($test->{'CELLS'}->{$cell}->{'COUNT'} or 0); #number of overlaps for the test MVPs

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
            $pbinom += binomial($k, $num_of_valid_probes, $p);
        }
    } else {
        foreach my $k ($teststat .. $num_of_valid_probes) {
            $pbinom += binomial($k, $num_of_valid_probes, $p);
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

    push(@results, [$zscore, $pbinom, $shortcell, $$samples{$cell}{'tissue'}, $$samples{$cell}{'datatype'}, $$samples{$cell}{'file'}, $mvp_string, $$samples{$cell}{'acc'}]);
}
if (!$web) {
    close(BACKGROUND);
}

## ============================================================================
## Correct the p-values for multiple testing using the Benjamini-Yekutieli FDR control method
## ============================================================================
my $qvalues = BY(\@pvalues);
$qvalues = [map {sprintf("%.2e", $_)} @$qvalues];
## ============================================================================


## ============================================================================
## Write the results to a tab-separated file
## ============================================================================
my $results_filename = "$lab.chart.tsv.gz";
open(TSV, "| gzip -9 > $out_dir/$results_filename") or die "Cannot open $out_dir/$results_filename: $!";
print TSV join("\t", "Zscore", "Pvalue", "Cell", "Tissue", "Datatype", "File", "Probe", "Accession", "Qvalue"), "\n";
for (my $i = 0; $i < @results; $i++) {
    print TSV join("\t", @{$results[$i]}, $qvalues->[$i]), "\n";
}
close(TSV);
## ============================================================================


## ============================================================================
## Generate plots
## ============================================================================
warn "[".scalar(localtime())."] Generating plots...\n";
unless (defined $noplot){
    #Plotting and table routines
    Chart($results_filename, $lab, $out_dir, $samples, $cells, $label, $t_marginal, $t_strict, $dataset); # basic pdf plot
    dChart($results_filename, $lab, $out_dir, $dataset, $label, $t_marginal, $t_strict, $web); # rCharts Dimple chart
    table($results_filename, $lab, $out_dir, $web); # Datatables chart
  }
## ============================================================================

warn "[".scalar(localtime())."] Done.\n";

if ($autoopen) {
    system("open $out_dir/$lab.table.html");
    system("open $out_dir/$lab.dchart.html");
    system("open $out_dir/$lab.chart.pdf");
}


####################################################################################################
####################################################################################################
##
##  Sub-functions
##
####################################################################################################
####################################################################################################


=head2 parse_pvalue_thresholds

 Arg[1]         : string $thresholds
 Returns        : arrayref of marginal and strict thresholds (floats)
 Example        : ($t_marginal, $t_strict) = parse_pvalue_thesholds("0.05,0.01");
 Description    : This function returns the both marginal and strict p-value thresholds as read from
                  the command line option. The input string should contain both numbers separated by
                  a comma.
 Exceptions     : Dies if $thresholds is empty, does not contain numbers or are not defined between
                  0 and 1 and/or the marginal threshold is not larger or equal to the strict one.

=cut

sub parse_pvalue_thresholds {
    my ($thresh) = @_;
    my ($t_marginal, $t_strict);

    if (!$thresh) {
        die "Cannot read p-value thresholds from an empty string\n";
    }

    ($t_marginal, $t_strict) = split(",", $thresh);
    unless (looks_like_number($t_marginal) && looks_like_number($t_strict)){
        die "You must specify numerical p-value thresholds in a comma separated list\n";
    }
    unless ((1 >= $t_marginal) && ($t_marginal >= $t_strict) && ($t_strict >= 0)) {
        die "The p-value thresholds should be 1 >= T.marginal >= T.strict >= 0\n";
    }
    return ($t_marginal, $t_strict);
}


=head2 get_input_probes

 Arg[1]         : string $filename
 Arg[2]         : arrayref $probe_list
 Returns        : arrayref of probe IDs (string)
 Example        : $mvps = get_input_probes("input.txt", undef);
 Example        : $mvps = get_input_probes(undef, ["cg13430807", "cg10480329,cg06297318,cg19301114"]);
 Example        : $mvps = get_input_probes(undef, undef);
 Description    : This function returns the list of input probe IDs. This can come from either
                  $filename if defined or from $probe_list otherwise. Each element in $probe_list is a
                  string which contains one or more probe IDs separated by commas (see Examples).
                  Falls back to the default data set from Jaffe and Irizarry.
                  The set of probe IDs is checked to remove redundant entries.
 Exceptions     : Dies if the file is not found or cannot be opened for whatever reason.

=cut

sub get_input_probes {
    my ($filename, $probe_list) = @_;
    my $probes;

    if (defined $filename) {
        my $fh;
        if ($filename =~ /\.gz$/) {
            open($fh, "gunzip -c $filename |") or die "cannot open file $filename : $!";
        } elsif ($filename =~ /\.bz2$/) {
            open($fh, "bunzip2 -c $filename |") or die "cannot open file $filename : $!";
        } else {
            open($fh, "$filename") or die "cannot open file $filename : $!";
        }
        $probes = process_file($fh, $format, $dbh, $array, $filter);

    } elsif ($probe_list and @$probe_list) {
        @$probes = split(/,/, join(',', @$probe_list));

    } else{
        # Test MVPs from Liu Y et al. Nat Biotechnol 2013  Pulmonary_function.snps.bed (*put EWAS bedfile here)
        # If no options are given it will run on the default set of MVPs
        warn "No probe input given, so running on default set of probes, a set of monocyte tDMPs from Jaffe AE and Irizarry RA, Genome Biol 2014.";
        @$probes = qw(cg13430807 cg10480329 cg06297318 cg19301114 cg23244761 cg26872907 cg18066690 cg04468741 cg16636767 cg10624395 cg20918393);
    }

    # Remove redundancy in the input
    my %probes_hash;
    foreach my $probe (@$probes) {
        $probes_hash{$probe}++;
    }

    while (my ($probe, $num) = each %probes_hash) {
        if ($num > 1) {
            say "$probe is present $num times in the input. Analysing only once."
        }
    }

    @$probes = keys %probes_hash;

    return($probes);
}


=head2 find_missing_probes

 Arg[1]         : arrayref of strings $original_probe_ids
 Arg[2]         : hashref $existing_probe_ids (keys are probe_ids, values are ignored)
 Returns        : int $num__missing_probes
 Example        : my $num_missing_probes = find_missing_probes(['cg001', 'cg002', 'cg003', 'cg004'],
                    {'cg001' => 1, 'cg003 => 1});
 Description    : Detects and prints the list of missing probes if any.
 Exceptions     :

=cut

sub find_missing_probes {
    my ($original_probes, $existing_probes_hash) = @_;
    my $num_missing_probes = 0;

    my $missing_probes = [];
    foreach my $probe_id (@$original_probes) {
        unless ($existing_probes_hash->{$probe_id}) {
            push @$missing_probes, $probe_id;
        }
    }
    $num_missing_probes = scalar @$missing_probes;

    if ($num_missing_probes > 0) {
        warn "The following $num_missing_probes MVPs have not been analysed because they were not found on the selected array\n";
        warn join("\n", @$missing_probes) . "\n";
    }

    return $num_missing_probes;
}
