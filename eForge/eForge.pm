package eForge::eForge;

=head1 NAME

eForge::eForge - Interface with the DB and various other common functions for eForge

=head1 VERSION

Version 0.01

=head1 LICENCE AND COPYRIGHT

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

=cut

use 5.010;
use strict;
use warnings FATAL => 'all';
use Storable;


my $MAX_SQL_VARIABLES = 999;
our $VERSION = '0.01';
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_all_datasets get_all_arrays get_all_proxy_filters process_file match process_bits get_bits get_cells assign bkgrdstat prox_filter);

=head1 SYNOPSIS

Provide functional modules for eForge

=head1 EXPORT

get_all_datasets
get_all_arrays
get_all_proxy_filters
process_file
match
process_bits
get_bits
get_cells
assign
bkgrdstat
prox_filter

=head1 SUBROUTINES/METHODS

=head2 bkgrdstat
 
Gives the background statistics.

=cut 

sub bkgrdstat{
    my ($test, $lab, $flag) = @_;
    my $bfh;
    my $file = "$lab.bkgrd.stats";
    if ($flag eq "test") {
        open $bfh, ">", $file or die "cannot open $file";
      }
    else{
        open $bfh, ">>", $file or die "cannot open $file";
      }
    my (@feature, @cpg_island_relationship);
    foreach my $probeid (keys %{$$test{'MVPS'}}){
        my ($feature, $cpg_island_relationship) = split "\t", $$test{'MVPS'}{$probeid}{'PARAMS'};
        push @feature, $feature;
        push @cpg_island_relationship, $cpg_island_relationship;
      }
    say $bfh join("\t", $flag, "feature", @feature);
    say $bfh join("\t", $flag, "cpg_island_relationship", @cpg_island_relationship);
  }


=head2 process_file

 Arg[1]         : FILEHANDLE $input_fh
 Arg[2]         : string $format
 Arg[3]         : DB-connection-handle $dbh
 Arg[4]         : string $array
 Arg[5]         : numeric $filter
 Returns        : arrayref of $probe_ids (string)
 Example        : my $probe_ids = process_file($fh, "probeid", $dbh, '450k', undef);
 Description    : Reads the list of probe IDs or locations from the file. If the file contains
                  locations, these will be translated into probe_ids
 Exceptions     :

=cut

sub process_file {
    my ($fh, $format, $dbh, $array, $filter) = @_;
    my $probe_ids = [];

    if ($format =~ /^probe/i) {
        while (<$fh>) {
            next if /^#/;
            chomp;
            my $probe_id;
            if (defined $filter) {
                my $pval;
                ($probe_id, $pval) = split /\s+/, $_;
                next unless $pval >= $filter;
            } else {
                ($probe_id, undef) = split /\s+/, $_; # remove anything that is not supposed to be there :-)
            }
            push @$probe_ids, $probe_id;
        }

    } elsif ($format =~ /^bed/i) {
        if (defined $filter) {
            warn "You have specified p-value filtering, but this isn't implemented for files of format $format. No filtering will happen."
        }
        my $locations = [];
        while (<$fh>) {
            next if /^track/;
            chomp;
            my ($chr, $from, $to) = split "\t", $_;
            next if (!defined($to));
            unless ($chr =~ /^chr/){
                $chr = "chr". $chr;
            }
            push(@$locations, [$chr, $from]);
        }
        $probe_ids = fetch_all_probe_ids($dbh, $array, $locations);

    } elsif ($format =~ /^tabix/i) {
        if (defined $filter) {
            warn "You have specified p-value filtering, but this isn't implemented for files of format $format. No filtering will happen."
        }
        my $locations = [];
        while (<$fh>) {
            chomp;
            my ($chr, $from, $to) = $_ =~ /(.+)\:(\d+)\-(\d+)/;
            push(@$locations, [$chr, $from]);
        }
        $probe_ids = fetch_all_probe_ids($dbh, $array, $locations);
    }
    return $probe_ids;
}

=head2 match

Identifies the bins that each of the probes in a probe hash lies in, and then picks matching probes for the number of reps specified.

=cut

#it is not enough to change match, we also have to change bkgrdstat and process_bits and get_bits
#(we don't use "assign" as we do not use percentile bins so no point in changing "assign")
sub match{
    #we take out the percentile bins ($per):
    #my ($mvps, $array, $datadir, $per, $reps) = @_;
    my ($mvps, $array, $datadir, $reps) = @_;
    #my ($bins, $params, %bins, %params);
    my ($bins, %bins);
    # load up the stored hashes that contain the bins of mvps by feature and cpg island relationship.
    # These are precalculated according to the parameters that are hard coded above.
    # the hash to load is defined by the bkgd option - defaults to '450k'
    #$bins = $datadir . "/mvp_bins";
    #$bins = $datadir . "snp_bins.$per";
    #$params = $datadir . "snp_params.$per";

    if ($array =~ "27k"){
        $bins = $datadir . "/mvp_27k_bins";
        }
    else{
        $bins = $datadir . "/mvp_450k_bins";

        }

    #took params out, do not need params
    #if (-e $bins && -e $params){
    if (-e $bins){
        %bins = %{ retrieve($bins) };
        #took params out, do not need params
        #%params = %{ retrieve($params)};
      }
    else{
        die "Cannot retrieve the file $bins\n";
    }
    
    my (%picks);

    

    foreach my $cg (keys %{$$mvps{'MVPS'}}){
        srand;
        my ($feature, $cpg_island_relationship) = split("\t", join("\t", $$mvps{'MVPS'}{$cg}{'PARAMS'}));
        #$cg is the test mvp, $cgid is the matched mvp.

        #range has to be the number of probes to choose from in that hash subclass
        my $range=scalar @{$bins{$feature}{$cpg_island_relationship}};

        for (my $n = 1; $n <= $reps; $n++) {
            my ($mvp_string, $cgid);
            while (1){
                my $pick = int(rand($range));
                $mvp_string = ${$bins{$feature}{$cpg_island_relationship}}[$pick]; #pick the $pick'th element in the array as the chosen mvp
                #(undef, undef, undef, $cgid) = split /\t/,  $mvp_string; #did not set the splitting on, check whether this can be done
                $cgid=$mvp_string;
                last unless $cgid eq $cg; # must not pick the test mvp itself.
              }
            push @{$picks{$n}}, $cgid; # each $n array is a set of probes matching the test set/ it is allowed to pick the same probe more than once in this background selection
          }
      }
    return \%picks;
  }

#commented previous foreach with three parameters(maf, tss, gc content) can be found below:
#foreach my $cg (keys %{$$mvps{'SNPS'}}){
#        srand;
#        my ($maf, $tss, $gc) = split("\t", join("\t", $$mvps{'SNPS'}{$cg}{'PARAMS'}));
#        #$cg is the test mvp, $cgid is the matched mvp.
#        my ($i, $j, $k) = assign ($gc, $tss, $maf, \%params);
#
#        my $range = scalar @{$bins{$i}{$j}{$k}};
#        for (my $n = 1; $n <= $reps; $n++) {
#            my ($mvp_string, $cgid);
#            while (1){
#                my $pick = int(rand($range));
#                $mvp_string = ${$bins{$i}{$j}{$k}}[$pick]; #pick the $pick'th element in the array as the chosen mvp "
#                (undef, undef, undef, $cgid) = split /\t/,  $mvp_string;
#                last unless $cgid eq $cg; # must not pick the test mvp itself.
#            }
#            push @{$picks{$n}}, $cgid; # each $n array is a set of probes matching the test set/ it is allowed to pick the same probe more than once in this background selection
#        }
#    }
#    return \%picks;




=head2 process_bits

 Arg[1]         : arrayref of arrays $rows ($probe_id, $sum, $bit_string, $feature, $CGI_context)
 Arg[2]         : arrayref of strings $cells (shortname for the cells in the dataset)
 Arg[3]         : string $dataset
 Returns        : hashref $stats
 Example        : my $result = process_bits($rows, $cells, 'erc');
 Description    : Returns a reference to a complex hash with stats from the rows. These are split
                  into 'MVPS' and 'CELLS'. The former contains 'SUM' and 'PARAMS' for each probe ID
                  while the latter contains 'COUNT' and 'MVPS' for each cell.
 Exceptions     : Dies if the number of cells does not match the length of the bit string.

=cut

sub process_bits {
    my ($rows, $cells, $data) = @_;
    my %test;
    my @test_cells;
    my @indexes = 0..(@$cells-1);
    foreach my $row (@{$rows}){
        my ($probeid, $sum, $bit_string, $feature, $cpg_island_relationship) = @$row;
        $test{'MVPS'}{$probeid}{'SUM'} = $sum;
        $test{'MVPS'}{$probeid}{'PARAMS'} = join("\t", $feature, $cpg_island_relationship);
        die "For $data, found ".scalar(@$cells)." cells for ".length($bit_string)." bits\n" if (scalar(@$cells) ne length($bit_string));
        foreach my $index (@indexes) {
            ## $bit_string is a string made of 0s and 1s. If it is a 1 for this position, count and push
            if (substr($bit_string, $index, 1)) {
                $test_cells[$index][0]++;
                push @{$test_cells[$index][1]}, $probeid;
            }
        }
    }
    my $index = 0;
    foreach my $cell (@$cells){
        $test{'CELLS'}{$cell}{'COUNT'} = $test_cells[$index][0] if ($test_cells[$index][0]);
        $test{'CELLS'}{$cell}{'MVPS'} = $test_cells[$index][1] if ($test_cells[$index][1]);
        $index++;
    }

    return \%test;
}


=head2 get_all_datasets

 Arg[1]         : DB-handle $dbh
 Arg[2]         : (optional) string $species_name
 Returns        : arrayref of hashes (tag/name)
 Example        : my $all_datasets = get_all_dataset($dbh);
 Description    : Returns the list of all datasets in the DB. It returns an arrayref of hashes. Each
                  hash has two keys: 'tag' (string ID of the dataset) and 'name' (full name). You
                  can limit the dataset to the ones available for a given species.
 Exceptions     :

=cut

sub get_all_datasets {
    my ($dbh, $species_name) = @_;
    my $datasets;

    my $sql = "SELECT dataset_tag, dataset_name FROM dataset ORDER BY dataset_id";
    my @bind_params = ();
    if ($species_name) {
        $sql .= " WHERE species_name = ?";
        push(@bind_params, $species_name);
    }
    my ($tag, $name);
    my $sth = $dbh->prepare($sql);
    $sth->execute(@bind_params);
    $sth->bind_columns(\$tag, \$name);
    while ($sth->fetch) {
        push(@$datasets, {tag => $tag, name => $name});
    }

    return($datasets);
}


=head2 get_all_arrays

 Arg[1]         : DB-handle $dbh
 Arg[2]         : (optional) string $species_name
 Returns        : arrayref of hashes (tag/name)
 Example        : my $all_arrays = get_all_arrays($dbh);
 Description    : Returns the list of all arrays in the DB. It returns an arrayref of hashes. Each
                  hash has two keys: 'tag' (string ID of the array) and 'name' (full name). You
                  can limit the arrays to the ones available for a given species.
 Exceptions     :

=cut

sub get_all_arrays {
    my ($dbh, $species_name) = @_;
    my $arrays;

    my $sql = "SELECT array_tag, array_name FROM array order by array_id";
    my @bind_params = ();
    if ($species_name) {
        $sql .= " WHERE species_name = ?";
        push(@bind_params, $species_name);
    }
    my ($tag, $name);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$tag, \$name);
    while ($sth->fetch) {
        push(@$arrays, {tag => $tag, name => $name});
    }

    return($arrays);
}


=head2 get_all_proxy_filters

 Arg[1]         : DB-handle $dbh
 Arg[2]         : (optional) string $array
 Returns        : hashref of proxy-filters (string)
 Example        : my $proxy_filters = get_all_proxy_filters($dbh, '450k');
 Description    : Returns the list of proxy filters available in the DB. It returns a hashref whose
                  keys are the bkgd names (i.e. array tags) and values are the name of the proxy
                  filter. At the moment the schema supports just one proxy filter per array.
                  If you provide a $array, the result will be limited to that array. Note that you
                  will still get a hashref with exactly the same structure. Only the content will
                  be limited.
 Exceptions     :

=cut

sub get_all_proxy_filters {
    my ($dbh, $array_tag) = @_;
    my $proxy_filters;

    my $sql = "SELECT array_tag, description FROM proxy_filter_info JOIN array USING (array_id)";
    my @bind_params = ();
    if ($array_tag) {
        $sql .= " WHERE array_tag = ?";
        push(@bind_params, $array_tag);
    }
    my ($this_array_tag, $description);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$this_array_tag, \$description);
    while ($sth->fetch) {
        $proxy_filters->{$this_array_tag} = $description;
    }

    return($proxy_filters);
}

=head2 get_bits

Get the bitstrings for an array of mvps from the sqlite db

=cut

sub get_bits{
    my ($dbh, $dataset_tag, $array_tag, $mvps) = @_;
    my @results;
    
    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$mvps; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$mvps -1 if ($end >= @$mvps);

        my $sql = "SELECT probe_id, sum, bit, gene_group, cgi_group
            FROM probe_annotation
            JOIN probe_bitstring USING (array_id, probe_id)
            JOIN dataset USING (dataset_id)
            JOIN array ON (array.array_id = probe_annotation.array_id)
            WHERE array_tag = ? and dataset_tag = ?
            AND probe_id IN (?". (",?" x ($end - $start)).")";
        my $sth = $dbh->prepare_cached($sql); #get the blocks form the ld table
        $sth->execute($array_tag, $dataset_tag, @$mvps[$start..$end]);
        
        while (my $row = $sth->fetchrow_arrayref()) {
            push @results, [@$row];
        }
        $sth->finish();
    }

    return \@results;# return the bitstring line from the database
}


=head2 fetch_all_probe_ids

 Arg[1]         : DB-handle $dbh
 Arg[2]         : string $array
 Arg[3]         : arrayref of locations [$chr, $pos]
 Returns        : arrayref of probe IDs (string)
 Description    : Given a list of chromosome names and location pairs, the method fetches the name
                  of the probe ID for that location. To do this successfully, this method requires
                  to know the background (i.e. array) you want to translate the locations to.
 Exceptions     :

=cut

sub fetch_all_probe_ids {
    #gets the rsid for a SNP where a location is given
    my ($dbh, $array, $locations) = @_;

    my $sth = $dbh->prepare("SELECT probe_id
                             FROM probe_mapping
                             WHERE  chr_name = ? AND position = ?");
    my $probe_ids = [];

    foreach my $this_loc (@$locations) {
        my ($chr, $pos) = @$this_loc;
        $sth->execute($chr, $pos);
        my $result = $sth->fetchall_arrayref();
        my $probe_id;
        foreach my $row (@{$result}) {
            push(@$probe_ids, $$row[0]);
        }
    }
    $sth->finish();
    return $probe_ids;
}

=head2 prox_filter

Filter MVPs from the MVP list if they are within 1 kb of each other. The rationale is that the first MVP to be identified in a block is chosen, and others are removed.

=cut

sub prox_filter{
    my ($mvps, $dbh) = @_;
    my %prox_excluded; # a hash to store MVPs found in proximity (1 kb) with an MVP in the list
    my @mvps_filtered; # The list of MVPs filtered
    my %mvps;
    foreach my $mvp (@$mvps){
        $mvps{$mvp} = 1;
      }
    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$mvps; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$mvps -1 if ($end >= @$mvps);

        my $sql = "SELECT probe_id,proxy_probes FROM proxy_filter WHERE probe_id IN (?". (",?" x ($end - $start)).")";
        my $sth = $dbh->prepare($sql); #get the blocks form the ld table
        $sth->execute(@$mvps[$start..$end]);
        my $result = $sth->fetchall_arrayref();
        $sth->finish();
        foreach my $row (@{$result}){
            my ($mvp, $block) = @$row;
            next if exists $prox_excluded{$mvp}; # if the mvp is in the proximity filtered set already ignore it
            push @mvps_filtered, $mvp; # if this is the first time it is seen, add it to the filtered mvps, and remove anything in proximity with it
            next if $block =~ /NONE/; # nothing is in proximity
            my (@block) = split (/\|/, $block);
            foreach my $proxmvp (@block){
                if (exists $mvps{$proxmvp}) {
                    $prox_excluded{$proxmvp} = $mvp; #Add to the excluded mvps, if it is in proximity with the current mvp, and it its one of the test mvps.
                  }
              }
          }
      }
    #note that if an MVP doesn't exist in the proximity file it will be rejected regardless, may need to add these back
    return (\%prox_excluded, @mvps_filtered);
}

=head2 get_cells

Read the correct cell list based on data (erc, erc2, blueprint or encode). Also gets the tissue names for the cells.

=cut

sub get_cells {
    my ($dataset_tag, $dbh) = @_;
    my $samples;

    my $sth = $dbh->prepare("SELECT shortcell, tissue, datatype, file, acc FROM dataset JOIN sample USING (dataset_id) WHERE dataset_tag = ? ORDER BY sample_order");
    $sth->execute($dataset_tag);
    my ($shortcell, $tissue, $datatype, $file, $acc);
    $sth->bind_columns(\$shortcell, \$tissue, \$datatype, \$file, \$acc);
    my ($cells, $tissues);
    while ($sth->fetch) {
        $shortcell = "$shortcell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $shortcell;
        $$tissues{$shortcell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$tissues{$shortcell}{'datatype'} = $datatype;
        $$tissues{$shortcell}{'file'} = $file;
        $$tissues{$shortcell}{'acc'} = $acc;
    }
#     use Data::Dumper;
#     print Dumper %$tissues;

    return ($cells, $tissues); # return
}

1;
