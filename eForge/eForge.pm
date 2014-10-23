package eForge::eForge;

use 5.010;
use strict;
use warnings FATAL => 'all';
use Storable;

=head1 NAME

eForge - The great new eForge!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(process_file match process_bits get_bits fetch_rsid fetch_loc get_cells assign bkgrdstat prox_filter);

=head1 SYNOPSIS

Provide functional modules for eForge

=head1 EXPORT

process_file
match
process_bits
get_bits
fetch_rsid
fetch_loc
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

Processes various file formats.

=cut

sub process_file {
    my ($fh, $format, $sth, $filter) = @_;
    my @snps;
    if ($format =~ /probeid/){
      while (<$fh>){
            chomp;
            my $rs;
            if (defined $filter) {
                my $pval;
                ($rs, $pval) = split /\s+/, $_;
                next unless $pval >= $filter;
              }
            else{
                ($rs, undef) = split /\s+/, $_; # remove anything that is not supposed to be there :-)
              }
            my @rsid = split /\:/, $rs;
            my $rsid = pop @rsid; # take the last one for want of a better idea.
            push @snps, $rsid;
          }
      }
    elsif ($format =~ /ian/){
      while (<$fh>){
            my ($chr, $beg, $end, $rsid, $p, $pval) = split "\t", $_;
            if (defined $filter) {
                next unless $pval >= $filter;
              }
            my @rsid = split /\:/, $rsid; # to deal with multiple rsids
            $rsid = pop @rsid; # take the last one for want of a better idea. Can't take all as they are the same thing.
            push @snps, $rsid;
          }
      }
    elsif ($format =~ /vcf/){
      while (<$fh>){
            next if /^#/;
            my ($chr, $beg, $rsid) = split "\t", $_;
            unless ($chr =~ /^chr/){
                $chr = "chr". $chr;
              }
            if ($rsid =~ /^cg\d+/){
                push @snps, $rsid;
              }
            else {
                my $loc = "$chr:$beg-$beg";
                #get the rsid from the db
                $rsid = fetch_rsid($loc, $sth);
                push @snps, $rsid if defined $rsid;
              }
          }
      }
    elsif ($format =~ /bed|tabix/){
      while (<$fh>){
            my $loc;
            if ($format =~/bed/){
                next if /^track/;
                chomp;
                my ($chr, $beg, $end) = split "\t", $_;
                unless ($chr =~ /^chr/){
                    $chr = "chr". $chr;
                  }
                $loc = "$chr:$end-$end";
              }
            elsif ($format =~ /tabix/){
                chomp;
                $loc = $_;
              }
            #get the $rsid from the db
            my $rsid = fetch_rsid($loc, $sth);
            push @snps, $rsid if defined $rsid;
          }
      }
    return @snps;
  }

=head2 match

Identifies the bins that each of the probes in a probe hash lies in, and then picks matching probes for the number of reps specified.

=cut

#it is not enough to change match, we also have to change bkgrdstat and process_bits and get_bits
#(we don't use "assign" as we do not use percentile bins so no point in changing "assign")
sub match{
    #we take out the percentile bins ($per):
    #my ($mvps, $bkgd, $datadir, $per, $reps) = @_;
    my ($mvps, $bkgd, $datadir, $reps) = @_;
    #my ($bins, $params, %bins, %params);
    my ($bins, %bins);
    # load up the stored hashes that contain the bins of mvps by feature and cpg island relationship.
    # These are precalculated according to the parameters that are hard coded above.
    # the hash to load is defined by the bkgd option - defaults to '450k'
    $bins = $datadir . "/mvp_bins";
    #$bins = $datadir . "snp_bins.$per";
    #$params = $datadir . "snp_params.$per";
    
    #took params out, do not need params
    #if (-e $bins && -e $params){
    if (-e $bins){
        %bins = %{ retrieve($bins) };
        #took params out, do not need params
        #%params = %{ retrieve($params)};
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

Processes bitstrings to get a count of overlaps for each cell type.

=cut

sub process_bits{
    my ($rows, $cells, $data) = @_;
    my %test;
    foreach my $row (@{$rows}){
        my ($location, $probeid, $sum, $bit, $feature, $cpg_island_relationship);
        if ($data eq "erc"){
            ($location, $probeid, undef, undef, $bit, $sum, $feature, $cpg_island_relationship) =  @$row;
          }
        else{
            ($location, $probeid, $bit, $sum, undef, undef, $feature, $cpg_island_relationship) = @$row;
          }
        $test{'MVPS'}{$probeid}{'SUM'} = $sum;
        $test{'MVPS'}{$probeid}{'PARAMS'} = join("\t", $feature, $cpg_island_relationship);
        my @bits = split "", $bit;
        my $index = 0;
        foreach my $cell (@$cells){
            $test{'CELLS'}{$cell}{'COUNT'} += $bits[$index];
            push @{$test{'CELLS'}{$cell}{'MVPS'}}, $probeid if $bits[$index] == 1;
            $index++;
          }
      }
    return \%test;
  }


=head2 get_bits

Get the bitstrings for an array of mvps from the sqlite db

=cut

sub get_bits{

    my ($mvps, $dbh) = @_;
    my @results;
    my $sql = "SELECT * FROM bits WHERE probeid IN (?". (",?" x (@$mvps - 1)).")";
    my $sth = $dbh->prepare($sql); #get the blocks form the ld table
    $sth->execute(@$mvps);
    my $result = $sth->fetchall_arrayref();
    $sth->finish();
    foreach my $row (@{$result}){
      push @results, $row;
    }
    return \@results;# return the bitstring line from the database
  }

=head2 fetch_rsid

gets the rsid for a SNP where a location is given.

=cut

sub fetch_rsid{
    #gets the rsid for a SNP where a location is given
    my ($loc, $sth) = @_;
    $sth->execute($loc);
    my $result = $sth->fetchall_arrayref();
    my $rsid;
    foreach my $row (@{$result}){
        $rsid = $$row[0];
      }
    $sth->finish();
    if (defined $rsid &&$rsid =~ /^cg\d+/){
        return $rsid;
      }
    else{
        return "no PROBEID match for $loc";
      }
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
    my $sql = "SELECT probeid,proxy_probes FROM proxy_filter WHERE probeid IN (?". (",?" x (@$mvps - 1)).")";
    my $sth = $dbh->prepare($sql); #get the blocks form the ld table
    $sth->execute(@$mvps);
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
                say "$proxmvp excluded for proximity (1 kb) with $mvp";
              }
          }
      }
    #note that if an MVP doesn't exist in the proximity file it will be rejected regardless, may need to add these back
    return (\%prox_excluded, @mvps_filtered);
}


=head2 get_cells

Read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.

=cut

sub get_cells{
    # read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.
    my ($data, $dbh) = @_;

    my $table = "cells_".$data;

    # Check that the table exists in the DB (note, some magic here that might be SQLite-specific)
    my @tables = grep {/^cells_/} map {$_ =~ s/"//g; $_ =~ s/^main\.//; $_} $dbh->tables();
    if (!grep {/$table/} @tables) {
        die "The database does not contain information for the data background provided.\n";
    }

    my $sth = $dbh->prepare("SELECT shortcell,tissue,file,acc FROM $table");
    $sth->execute();
    my $ver = $sth->fetchall_arrayref();
    $sth->finish();
    my ($cells, $tissues, $acc);
    foreach my $row (@$ver){
        my $cell = shift @$row;
        my $tissue = shift @$row;
        my $file = shift @$row;
        my $acc = shift @$row;
        $cell = "$cell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $cell;
        $$tissues{$cell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$tissues{$cell}{'file'} = $file;
        $$tissues{$cell}{'acc'} = $acc;
      }
    #print Dumper %$tissues;
    return ($cells, $tissues); # return
  }

=head2 assign

Assign any maf, gc, tss values to the percentile bins

=cut

#sub assign{
#    #sub routine to assign any maf, gc, tss values to the percentile bins
#    my ($gc, $tss, $maf, $params) = @_;
#    my ($i, $j, $k);
#    my $n = 1;
#    foreach my $pc (@{$$params{'gc'}}){
#      if ($gc <= $pc) {
#            $i = $n;
#          }  
#        else{
#            $n++;
#          }
#    }
#    $n=1;
#    foreach my $pc (@{$$params{'tss'}}){
#      if ($tss <= $pc) {
#            $j = $n;
#          }
#        else{
#            $n++;
#          }
#    }
#    $n=1;
#    foreach my $pc (@{$$params{'maf'}}){
#      if ($maf <= $pc) {
#            $k = $n;
#          }
#        else{
#            $n++;
#          }
#    }
#    return ($i, $j, $k);
#  }

=head1 AUTHOR

Charles Breeze, C<< <cbreeze at ebi.ac.uk> >>
Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc eForge


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Charles Breeze.

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


=cut

1;
