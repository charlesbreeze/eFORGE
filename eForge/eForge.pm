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
use Data::Dumper;

my $MAX_SQL_VARIABLES = 99999;
our $VERSION = '0.01';
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_all_datasets get_all_arrays get_all_proximity_filters process_file get_random_matching_picks process_overlaps get_probe_annotations_and_overlap_for_dataset get_samples_from_dataset assign save_probe_annotation_stats proximity_filter);

=head1 SYNOPSIS

Provide functional modules for eForge

=head1 EXPORT

get_all_datasets
get_all_arrays
get_all_proximity_filters
process_file
get_random_matching_picks
process_overlaps
get_probe_annotations_and_overlap_for_dataset
get_samples_from_dataset
assign
save_probe_annotation_stats
proximity_filter

=head1 SUBROUTINES/METHODS


=head2 save_probe_annotation_stats

 Arg[1]         : hashref $overlaps (see get_overlaps)
 Arg[2]         : string $outdir
 Arg[3]         : string $label
 Arg[4]         : string $set_id
 Returns        : 
 Example        : save_probe_annotation_stats($overlaps, ".", "Unnamed", "test");
 Example        : save_probe_annotation_stats($overlaps, ".", "Unnamed", 23);
 Description    : Save stats about the probe annotations on a text file. The $set_id is either
                  "test" (which relates to the input probe_ids, the ones to be tested for
                  enrichment) or the random pick number.
                  The file will contain for each set the whole list of gene features and CpG islands
                  relationships for the probes in that set.
 Exceptions     : Dies if file cannot be opened

=cut

sub save_probe_annotation_stats {
    my ($overlaps, $out_dir, $lab, $flag) = @_;

    my $fh;
    my $file = "$out_dir/$lab.overlaps.stats.txt";
    open(STATS, ">>$file") or die "cannot open $file";
    my (@gene_features, @cpg_island_relationships);
    foreach my $probeid (keys %{$overlaps->{'DMPS'}}){
        my ($this_gene_feature, $this_cpg_island_relationship) =
            split("\t", $overlaps->{'DMPS'}->{$probeid}->{'PARAMS'});
        push @gene_features, $this_gene_feature;
        push @cpg_island_relationships, $this_cpg_island_relationship;
    }
    say STATS join("\t", $flag, "gene_features", @gene_features);
    say STATS join("\t", $flag, "cpg_island_relationships", @cpg_island_relationships);
    close(STATS);
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


=head2 get_random_matching_picks

 Arg[1]         : hashref $overlaps
 Arg[2]         : string $array_tag
 Arg[3]         : string $data_dir
 Arg[4]         : int $num_random_picks
 Returns        : arrayref of arrays of $probe_ids (string)
 Example        : my $random_picks = get_random_matching_picks($overlaps, "450k", ".", 1000);
 Description    : Get several random picks of probes matching the criteria defined in the $overlaps
                  hash. The random picks are selected from a pre-built hash stored in the $data_dir
                  called dmp_450k_bins (or so).
 Exceptions     : Cannot find the bins file for the selected array.

=cut

sub get_random_matching_picks {
    my ($overlaps, $array, $datadir, $num_random_picks) = @_;
    my $picks = [];

    # load up the stored hashes that contain the bins of dmps by feature and cpg island relationship.
    my %bins;
    my $bins_file = $datadir . "/dmp_${array}_bins";
    if (-e $bins_file) {
        %bins = %{ retrieve($bins_file) };
    } else {
        die "Cannot retrieve the file $bins_file\n";
    }
    
    foreach my $probe_id (keys %{$overlaps->{'DMPS'}}) {
        my ($feature, $cpg_island_relationship) = split("\t", join("\t", $overlaps->{'DMPS'}->{$probe_id}->{'PARAMS'}));

        #range has to be the number of probes to choose from in that hash subclass
        my $range = scalar @{$bins{$feature}{$cpg_island_relationship}};

        for (my $n = 0; $n < $num_random_picks; $n++) {
            my $picked_probe_id;
            while (1) {
                my $pick = int(rand($range));
                $picked_probe_id = ${$bins{$feature}{$cpg_island_relationship}}[$pick]; #pick the $pick'th element in the array as the chosen dmp
                last unless $picked_probe_id eq $probe_id; # must not pick the test dmp itself.
            }
            push(@{$picks->[$n]}, $picked_probe_id);
        }
    }

    return $picks;
}


=head2 process_overlaps

 Arg[1]         : arrayref of arrays $annotated_probes ($probe_id, $sum, $bit_string, $feature, $CGI_context)
 Arg[2]         : arrayref of strings $cells (shortname for the cells in the dataset)
 Arg[3]         : string $dataset
 Returns        : hashref $stats
 Example        : my $result = process_overlaps($rows, $cells, 'erc');
 Description    : Returns a reference to a complex hash with stats from the rows. These are split
                  into 'DMPS' and 'CELLS'. The former contains 'SUM' and 'PARAMS' for each probe ID
                  while the latter contains 'COUNT' and 'DMPS' for each cell.
 Exceptions     : Dies if the number of cells does not match the length of the bit string.

=cut

sub process_overlaps {
    my ($rows, $cells, $data) = @_;

    #print Dumper "rows";    
    #print Dumper \$rows;
    
    #print Dumper "cells";    
    #print Dumper \$cells;
    
    my $overlaps;
    my @overlapping_probes_per_cell;
    my @indexes = 0..(@$cells-1);
    foreach my $row (@{$rows}){
        my ($probeid, $sum, $bit_string, $feature, $cpg_island_relationship) = @$row;
        $overlaps->{'DMPS'}->{$probeid}->{'SUM'} = $sum;
        $overlaps->{'DMPS'}->{$probeid}->{'PARAMS'} = join("\t", $feature, $cpg_island_relationship);
        die "For $data, found ".scalar(@$cells)." cells for ".length($bit_string)." bits\n" if (scalar(@$cells) ne length($bit_string));
        
        # profile test 1
        
        #my @bits = map(int(chr), unpack("C*", $bit_string));
        #map{ if ($bits[$_] == 1) { push @{$overlapping_probes_per_cell[$_]}, $probeid } } 0..$#bits;
        
        # profile test 2
        
        #my @bits = unpack("C*", $bit_string);
        #map{ if ($bits[$_] == 49) { push @{$overlapping_probes_per_cell[$_]}, $probeid } } 0..$#indexes;
        
        # profile test 3 (https://tools.altiusinstitute.org/nytprof.033018.008)
        
        my $bsl = length($bit_string);
        my @bits = unpack("C$bsl", $bit_string);
        map{ if ($bits[$_] == 49) { push @{$overlapping_probes_per_cell[$_]}, $probeid } } 0..($bsl-1);
        
        # profile test 4
        
        #my $bsl = length($bit_string);
        #map{ if ((unpack("C$bsl", $bit_string))[$_] == 49) { push @{$overlapping_probes_per_cell[$_]}, $probeid } } 0..($bsl-1);
        
        # profile test 5 (https://tools.altiusinstitute.org/nytprof.033018.009)
        
        #my $bsl = length($bit_string);
        #my @bits = unpack("C$bsl", $bit_string);
        #foreach my $index (0..($bsl-1)) {
        #  if ($bits[$index] == 49) { push @{$overlapping_probes_per_cell[$index]}, $probeid }
        #}
        
        #
        # original (https://tools.altiusinstitute.org/nytprof.033018.007)
        #
        
        #foreach my $index (@indexes) {
        #    ## $bit_string is a string made of 0s and 1s. If it is a 1 for this position, count and push
        #    if (substr($bit_string, $index, 1)) {
        #        push @{$overlapping_probes_per_cell[$index]}, $probeid;
        #    }
        #}
    }
    my $index = 0;
    foreach my $cell (@$cells){
        if ($overlapping_probes_per_cell[$index] and @{$overlapping_probes_per_cell[$index]}) {
            $overlaps->{'CELLS'}->{$cell}->{'COUNT'} = scalar(@{$overlapping_probes_per_cell[$index]});
            $overlaps->{'CELLS'}->{$cell}->{'DMPS'} = $overlapping_probes_per_cell[$index];
        }
        $index++;
    }

    #print Dumper "overlaps";
    #print Dumper \$overlaps;

    return $overlaps;
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

    my $sql = "SELECT dataset_tag, dataset_name FROM dataset ORDER BY dataset_id DESC";
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

    my $sql = "SELECT array_tag, array_name FROM array order by array_id DESC";
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


=head2 get_all_proximity_filters

 Arg[1]         : DB-handle $dbh
 Arg[2]         : (optional) string $array
 Returns        : hashref of proxy-filters (string)
 Example        : my $proximity_filters = get_all_proximity_filters($dbh, '450k');
 Description    : Returns the list of proxy filters available in the DB. It returns a hashref whose
                  keys are the bkgd names (i.e. array tags) and values are the name of the proxy
                  filter. At the moment the schema supports just one proxy filter per array.
                  If you provide a $array, the result will be limited to that array. Note that you
                  will still get a hashref with exactly the same structure. Only the content will
                  be limited.
 Exceptions     :

=cut

sub get_all_proximity_filters {
    my ($dbh, $array_tag) = @_;
    my $proximity_filters;

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
        $proximity_filters->{$this_array_tag} = $description;
    }

    return($proximity_filters);
}


=head2 get_probe_annotations_and_overlap_for_dataset

 Arg[1]         : DB-handle $dbh
 Arg[2]         : string $dataset_tag
 Arg[3]         : string $array
 Arg[4]         : arrayref of strings $probe_ids
 Returns        : arrayref of arrays containing the probe_id, sum, bitstring, gene-group and
                  cgi-group
 Example        : my $annotated_probes = get_probe_annotations_and_overlap_for_dataset($dbh,
                        'erc', '450k', $probe_list);
 Description    : Fetches the gene and CGI related annotation for the set of probes as well as the
                  overlaps with the features defined in the selected dataset
 Exceptions     :

=cut

sub get_probe_annotations_and_overlap_for_dataset {
    my ($dbh, $dataset_tag, $array_tag, $probe_ids) = @_;
    my $results = [];

    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$probe_ids; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$probe_ids - 1 if ($end >= @$probe_ids);

        my $sql = "SELECT probe_id, sum, bit, gene_group, cgi_group
            FROM probe_annotation
            JOIN probe_bitstring USING (array_id, probe_id)
            JOIN dataset USING (dataset_id)
            JOIN array ON (array.array_id = probe_annotation.array_id)
            WHERE array_tag = ? and dataset_tag = ?
            AND probe_id IN (?". (",?" x ($end - $start)).")";

        my $sth = $dbh->prepare_cached($sql);
        $sth->execute($array_tag, $dataset_tag, @$probe_ids[$start..$end]);
        
        # profile test 1 (https://tools.altiusinstitute.org/nytprof.033018.011/)
        
        #my $max_rows = 100;
                
        # profile test 2 (https://tools.altiusinstitute.org/nytprof.033018.012/)
        
        #my $max_rows = 500;
                
        # profile test 3 (https://tools.altiusinstitute.org/nytprof.033018.013/)
        
        my $max_rows = 1000;
                
        # profile test 4 (https://tools.altiusinstitute.org/nytprof.033018.014/)
        
        #my $max_rows = 10000;
        
        my $rows = [];
        while (my $row = shift(@$rows) || shift ( @{$rows = $sth->fetchall_arrayref(undef, $max_rows)||[]} )) {
            push @$results, [@$row];
        }

        #
        # original (https://tools.altiusinstitute.org/nytprof.033018.010/)
        #
        
        #while (my $row = $sth->fetchrow_arrayref()) {
        #    push @$results, [@$row];
        #}
        
        $sth->finish();
    }

    return $results;
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
    my ($dbh, $array_tag, $locations) = @_;
    
    # Get array ID
    my $array_id_sth = $dbh->prepare("SELECT array_id FROM array WHERE array_tag = ?");
    $array_id_sth->execute($array_tag);
    my $array_id_result = $array_id_sth->fetchrow_arrayref();
    my $array_id = @{$array_id_result}[0];
    $array_id_sth->finish();
    
    # Get probe mapping ID
    my $probe_mapping_id_sth = $dbh->prepare("SELECT probe_mapping_id FROM probe_mapping_info WHERE array_id = ?");
    $probe_mapping_id_sth->execute($array_id);
    my $probe_mapping_id_result = $probe_mapping_id_sth->fetchrow_arrayref();
    my $probe_mapping_id = @{$probe_mapping_id_result}[0];
    $probe_mapping_id_sth->finish();

    my $sth = $dbh->prepare("SELECT probe_id
                             FROM probe_mapping
                             WHERE probe_mapping_id = ? AND chr_name = ? AND position = ?");
    my $probe_ids = [];

    foreach my $this_loc (@$locations) {
        my ($chr, $pos) = @$this_loc;
        $sth->execute($probe_mapping_id, $chr, $pos);
        my $result = $sth->fetchall_arrayref();
        my $probe_id;
        foreach my $row (@{$result}) {
            push(@$probe_ids, $$row[0]);
        }
    }
    $sth->finish();
    return $probe_ids;
}

=head2 proximity_filter

Filter DMPs from the DMP list if they are within 1 kb of each other. The rationale is that the first DMP to be identified in a block is chosen, and others are removed.

=cut

sub proximity_filter {
    my ($dbh, $array_tag, $probe_ids, $filter) = @_;
    my %prox_excluded_probes; # a hash to store DMPs found in proximity (1 kb) with a DMP in the list
    my %filtered_probes; # The list of DMPs filtered (i.e. after filtering, the ones to keep)
    my %missing_probes;
    
    # Get array ID
    my $array_id_sth = $dbh->prepare("SELECT array_id FROM array WHERE array_tag = ?");
    $array_id_sth->execute($array_tag);
    my $array_id_result = $array_id_sth->fetchrow_arrayref();
    my $array_id = @{$array_id_result}[0];
    $array_id_sth->finish();

    # Get the full list of probes as a hash (also removes redundancy)
    my %probe_id_hash;
    foreach my $probe_id (@$probe_ids){
        $probe_id_hash{$probe_id} = 1;
    }
    $probe_ids = [keys %probe_id_hash];
    
    # Debugging old sqlite3 constant
    #$MAX_SQL_VARIABLES = 999;

    for (my $loop = 0; $loop * $MAX_SQL_VARIABLES < @$probe_ids; $loop++) {
        my $start = $loop * $MAX_SQL_VARIABLES;
        my $end = ($loop + 1) * $MAX_SQL_VARIABLES - 1;
        $end = @$probe_ids - 1 if ($end >= @$probe_ids);

        # The proxy_filter_id value is array_id
        my $sql = "SELECT probe_id, proxy_probes FROM proxy_filter WHERE proxy_filter_id = $array_id AND probe_id IN (?". (",?" x ($end - $start)).")";
        my $sth = $dbh->prepare($sql); #get the blocks form the ld table
        $sth->execute(@$probe_ids[$start..$end]);
        my $result = $sth->fetchall_arrayref();
        $sth->finish();

        #my @probe_ids = ();
        foreach my $row (@{$result}){
            my ($probe_id, $probe_id_list) = @$row;
            #push @probe_ids, $probe_id;
            # if the probe is in the proximity filtered set already ignore it
            #next if exists $prox_excluded_probes{$probe_id}; 
            if (exists $prox_excluded_probes{$probe_id}) {
                #if ($probe_id eq "cg01914153") {
                #    print Dumper "cg01914153 found";
                #}
                next;
            } 
            # if this is the first time it is seen, add it to the filtered dmps, and remove anything in proximity with it
            $filtered_probes{$probe_id} = 1;
            next if $probe_id_list =~ /NONE/; # nothing is in proximity
            my (@other_probe_ids) = split (/\|/, $probe_id_list);
            foreach my $other_probe_id (@other_probe_ids) {
                if (exists $probe_id_hash{$other_probe_id}) {
                    $prox_excluded_probes{$other_probe_id} = $probe_id; #Add to the excluded dmps, if it is in proximity with the current dmp, and it its one of the test dmps.
                }
            }
        }
        #print Dumper join ',' , @probe_ids;
    }
    #die;
    
    #print Dumper join ',', sort keys %prox_excluded_probes;
    #die;
    
    #print Dumper join ',', sort keys %filtered_probes;
    #die;

    #note that if an DMP doesn't exist in the proximity file it will be rejected regardless, may need to add these back
    return (\%prox_excluded_probes, [keys %filtered_probes]);
}

=head2 get_samples_from_dataset

Read the correct cell list based on data (erc, erc2, blueprint or encode). Also gets the tissue names for the cells.

=cut

sub get_samples_from_dataset {
    my ($dbh, $dataset_tag) = @_;
    my ($cells, $samples);

    my $prepare = "SELECT shortcell, tissue, datatype, file, acc FROM dataset JOIN sample USING (dataset_id) WHERE dataset_tag = ? ORDER BY sample_order";
    my $sth = $dbh->prepare($prepare);
    $sth->execute($dataset_tag);
    my ($shortcell, $tissue, $datatype, $file, $acc);
    $sth->bind_columns(\$shortcell, \$tissue, \$datatype, \$file, \$acc);
    while ($sth->fetch) {
        $shortcell = "$shortcell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $shortcell;
        $$samples{$shortcell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$samples{$shortcell}{'datatype'} = $datatype;
        $$samples{$shortcell}{'file'} = $file;
        $$samples{$shortcell}{'acc'} = $acc;
    }
#     use Data::Dumper;
#     print Dumper %$tissues;

    return ($cells, $samples); # return
}

1;
