#! /usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use DBI;

my $db_dir = ".";
my $db_name = "eforge_1.1.db";
my $dataset_name = "ENCODE";
my $dataset_tag = "encode";
my $decode_file = "encode.decode";
my $work_dir = "/tmp";
my $species_name = "Homo sapiens";
my $assembly_name = "GRCh37";
my $bedtools = "bedtools";

my $help;

my $desc = qq{load_dataset.pl [options]

DESCRIPTION:

This script reads a 'decode' file which contains references to a list of samples, each of them
corresponding to a given dataset. Each of these files is a BED or BED-like file that is read by
bedtools to find overlaps with each and every array loaded in the eFORGE database.

Note that you *must* load the arrays first. If you want to include new arrays at a later date, you
will have to reload all the datasets again (i.e. you will have to re-start from scratch).

Required parameters:
 --tag <tag> or --dataset_tag <tag>
    the ID for this database. This needs to be unique.
 --name <name> or --dataset_name <name>
    the name for this database. This can be a longer description. It will be used in the web interface.
 --decode_file <file>
    the file with the information about each sample in this dataset

Optional parameters:
 --db_name <name>
    is the name of the SQLite file [def: $db_name]
 --db_dir <path>
    is the location of the SQLite file [def: $db_dir]
 --work_dir <path>
    is the location where temporary files will be downloaded/created [def: $work_dir]
 --bedtools <name>
    is the name of the bedtools executable you want to use [def: $bedtools]
 --species <species>
    is the name of the species [def: $species_name]
 --assembly <assembly>
    is the name of the assembly [def: $assembly_name]
};

GetOptions(
    "help" => \$help,
    "db_name=s" => \$db_name,
    "db_dir=s" => \$db_dir,
    "tag|dataset_tag=s" => \$dataset_tag,
    "name|dataset_name=s" => \$dataset_name,
    "decode_file=s" => \$decode_file,
    "work_dir=s" => \$work_dir,
    "species=s" => \$species_name,
    "assembly=s" => \$assembly_name,
    );

if ($help) {
    print $desc;
    exit(0);
}

my $dsn = "dbi:SQLite:dbname=$db_dir/$db_name";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;

system("mkdir -p $work_dir");

my $decode_table = get_decode_table($decode_file);

my $dataset_id = load_dataset($dbh, $species_name, $decode_table, $dataset_name, $dataset_tag);

download_bed_files($decode_table, $work_dir);

my $arrays = get_all_arrays_for_species($dbh, $species_name);

foreach my $this_array (@$arrays) {
    my ($this_array_id, $this_array_name) = @$this_array;
    my $sorted_array_bed_file = dump_array_bed_file($dbh, $this_array_id, $assembly_name, $work_dir);
    my $single_overlaps_bed_files = run_single_overlap_bedtools($bedtools, $work_dir, $sorted_array_bed_file, $decode_table, $this_array_id, $this_array_name);
    my $concatenated_overlaps_bed_file = paste_files($single_overlaps_bed_files, $work_dir, $this_array_id);
    my $final_overlaps_bed_file = add_sum_column_to_concatenated_bedfile($concatenated_overlaps_bed_file, $work_dir, $this_array_id);
    load_bitstrings($dbh, $db_name, $final_overlaps_bed_file, $this_array_id, $dataset_id);
}

# 


exit();

my $input_bed_files = get_input_bed_files_from_decode_table($decode_file);

exit(0);

sub load_dataset {
    my ($dbh, $species_name, $decode_table, $dataset_name) = @_;
    my $dataset_id;

    my $sth;
    $sth = $dbh->prepare("INSERT OR IGNORE INTO dataset (dataset_tag, dataset_name, species_name) VALUES (?, ?, ?)");
    $sth->execute($dataset_tag, $dataset_name, $species_name);
    $dataset_id = $dbh->last_insert_id("", "", "", "");
    $sth->finish();
    if ($dataset_id == 0) {
        $dataset_id = $dbh->selectrow_array("SELECT dataset_id FROM dataset WHERE dataset_name = '$dataset_name' AND species_name = '$species_name'");
    }


    $sth = $dbh->prepare("INSERT OR IGNORE INTO sample (dataset_id, sample_order, file, lab, datatype, cell, tissue, shortcell, individual, acc, url)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    my $sample_order = 1;
    foreach my $this_sample (@$decode_table) {
        $sth->execute($dataset_id,
            $sample_order,
            $this_sample->{file},
            $this_sample->{lab},
            $this_sample->{datatype},
            $this_sample->{cell},
            $this_sample->{tissue},
            $this_sample->{shortcell},
            $this_sample->{individual},
            $this_sample->{acc},
            $this_sample->{url});
        $sample_order++;
    }

    return($dataset_id);
}

sub load_bitstrings {
    my ($dbh, $db_name, $bed_file, $array_id, $dataset_id) = @_;

#     my $sql = "INSERT INTO probe_bitstring (array_id, probe_id, dataset_id, sum, bit)
#                 VALUES (?, ?, ?, ?, ?)";
#     my $sth = $dbh->prepare($sql);
    open(BED, $bed_file) or die "Cannot open BED file <$bed_file>\n";
    open(CSV, ">probe_bitstring.csv") or die "Cannot open CVS temporary file <probe_bitstring.csv>\n";
    while(<BED>) {
        chomp;
        my ($chr, $start, $end, $probe_id, $sum, $bitstring) = split("\t", $_);
        print CSV join(",", $array_id, $probe_id, $dataset_id, $sum, $bitstring), "\n";
#         $sth->execute($array_id, $probe_id, $dataset_id, $sum, $bitstring);
    }
    close(BED);
    close(CSV);
    system("echo '.mode csv
.import probe_bitstring.csv probe_bitstring' | sqlite3 $db_name");
#     $sth->finish();
}

=head2 add_sum_column_to_concatenated_bedfile

  Arg[1]        : string $concatenated_overlaps_bed_file (location of the input BED file with 0/1 flags on the 4th column)
  Arg[2]        : string $work_dir (where to put the temporary files)
  Example       : my $final_overlaps_bed_file = add_sum_column_to_concatenated_bedfile($concatenated_overlaps_bed_file, $work_dir);
  Description   : Reads the input BED file ($concatenated_overlaps_bed_file) which contains a series
                  of 0 and 1 flags in the 4th column. This function reads the number of ones in that
                  column and include that value in the output BED file. The output BED file will
                  contain that number in the 4th column and the series of flag in the 5th colum.
  Returns       : string $final_overlaps_bed_file (the location of the resulting BED file)
  Exceptions    : Dies if error when opening the files

=cut

sub add_sum_column_to_concatenated_bedfile {
    my ($concatenated_overlaps_bed_file, $work_dir, $array_id) = @_;
    my $final_overlaps_bed_file = "$work_dir/final_overlaps.array_${array_id}.bed";
    
    open(BED_IN, $concatenated_overlaps_bed_file) or die "Cannot open BED file <$concatenated_overlaps_bed_file>\n";
    open(BED_OUT, ">$final_overlaps_bed_file") or die "Cannot open BED file <$final_overlaps_bed_file>\n";
    while(<BED_IN>) {
        chomp;
        my ($this_chr, $this_start, $this_end, $this_probe, $this_bitstring) = split("\t", $_);
        my $num = $this_bitstring =~ tr/1/1/;
        print BED_OUT join("\t", $this_chr, $this_start, $this_end, $this_probe, $num, $this_bitstring), "\n";
    }
    close(BED_IN);
    close(BED_OUT);
    
    return $final_overlaps_bed_file;
}


=head2 run_single_overlap_bedtools

  Arg[1]        : string $bedtools (either full path or just the binary if in the $PATH)
  Arg[2]        : string $work_dir (where to put the temporary files)
  Arg[3]        : string $sorted_450k_bed_file (location of the BED file with sorted 450K features)
  Arg[4]        : arrayref of hash $decode_table
  Example       : my $single_overlaps_bed_files = run_single_overlap_bedtools($bedtools, $work_dir, $sorted_450k_bed_file, $decode_table);
  Description   : Run bedtools on the sorted 450K features vs all the BED DNAse features, one at a time.
  Returns       : arrayref of string $total_overlaps_bed_files (the locations of the resulting BED files)
  Exceptions    : Dies if error when running bedtools

=cut

sub run_single_overlap_bedtools {
    my ($bedtools, $work_dir, $sorted_array_bed_file, $decode_table, $array_id, $array_name) = @_;
    my $single_overlap_bed_files = [];

    foreach my $this_input_bed_file (map {$_->{"file"}} @$decode_table) {
        print "Overlap between $array_name and $this_input_bed_file...\n";
        my $this_output_bed_file = $this_input_bed_file;
        $this_output_bed_file =~ s/.+\///;
        $this_output_bed_file = "$work_dir/single_overlaps.array_${array_id}.$this_output_bed_file";
        my $runstr = "$bedtools intersect -c -a $sorted_array_bed_file -b $work_dir/$this_input_bed_file > $this_output_bed_file";
        system($runstr) == 0 or die "Error while running bedtools: $?";
        push(@$single_overlap_bed_files, $this_output_bed_file);
    }

    return $single_overlap_bed_files;
}
    

sub get_decode_table {
    my ($decode_file) = @_;
    my $decode_table;

    open(DECODE, $decode_file) or die "Cannot open decode file <$decode_file>\n";
    my $whole_text = join("\n", <DECODE>);
    $whole_text =~ s/[\r\n]+/\n/g;
    my @lines = split("\n", $whole_text);
    close(DECODE);
    
    my @header = split("\t", shift(@lines));
    
    foreach my $this_line (@lines) {
        my @data = split("\t", $this_line);
        my $decode_record  = {};
        for (my $i=0; $i<@header; $i++) {
            $decode_record->{$header[$i]} = $data[$i];
        }
        # Check that file and URL match
        die "File: ".$decode_record->{"file"}."\nURL: ".$decode_record->{"url"}."\nOn line: $this_line"
            if ($decode_record->{"url"} !~ $decode_record->{"file"});
        $decode_record->{"file"} = $decode_record->{"url"};
        $decode_record->{"file"} =~ s/.+\///g;
        push(@$decode_table, $decode_record);
    }
    
    return $decode_table;
}

sub download_bed_files {
    my ($decode_table, $work_dir) = @_;
    
    foreach my $this_decode_entry (@$decode_table) {
        my $url = $this_decode_entry->{"url"};
        my $file = $this_decode_entry->{"file"};
        if (!-e "$work_dir/$file") {
            download_url($url, $work_dir);
        }
    }
}


sub download_url {
    my ($url, $work_dir) = @_;

    print "Downloading $url...\n";
    system("wget", "-q", "-N", "-P", $work_dir, $url) == 0
        or die "Error: wget -q -N -P $work_dir $url\n$!";
}


sub paste_files {
    my ($files, $work_dir, $array_id) = @_;
    my $output_file = "$work_dir/concat_overlaps.array_${array_id}.bed";

    my $max_num_files = 200;
    my $c = 0;
    my $temp_output_file;
    my @original_files = @$files;
    while (@$files > 1) {
        $c++;
        $temp_output_file = "$work_dir/temp.$$.$c.paste_files.txt";
        my $input_files = [splice(@$files, 0, $max_num_files)];
#         print STDERR "Merging ", join(", ", @$input_files), " into $temp_output_file\n";
        merge_files($temp_output_file, $input_files);
        unshift(@$files, $temp_output_file);
    }
    rename($temp_output_file, $output_file);

    unlink glob "$work_dir/temp.$$.*.paste_files.txt";

# 
# 
# 
# 
#     unlink @original_files;
# 
# 
# 
# 

    return $output_file;
}


sub merge_files {
    my ($output_file, $input_files) = @_;
    my @fhs;
    my $c = 0;

    foreach my $this_input_file (@$input_files) {
        open($fhs[$c++], $this_input_file) or die "Cannot open $this_input_file: $!\n";
    }
    open(OUT, ">$output_file") or die "Cannot open $output_file: $!\n";

    my $first_fh = shift(@fhs);

    while (1) {
        my $line = <$first_fh>;
        chomp($line);
        my ($chr, $start, $end, $probe_id, $value) = split("\t", $line);
        print OUT join("\t", $chr, $start, $end, $probe_id, $value);
        foreach my $this_fh (@fhs) {
            my $line = <$this_fh>;
            chomp($line);
            my ($this_chr, $this_start, $this_end, $this_probe_id, $this_value) = split("\t", $line);
            if (($this_chr ne $chr) or ($this_start != $start) or ($this_end != $end) or ($this_probe_id ne $probe_id)) {
                die "Files do not contain the same $chr-$start-$end-$probe_id lines\n";
            }
            print OUT $this_value;
        }
        print OUT "\n";
        if (eof($first_fh)) {
            last;
        }
    }

    foreach my $this_fh (@fhs) {
        close($this_fh);
    }
    close(OUT);
}

sub get_all_arrays_for_species {
    my ($dbh, $species) = @_;
    my $arrays;

    my $sth = $dbh->prepare("SELECT array_id, array_name FROM array WHERE species_name = ?");
    $sth->execute($species);
    $arrays = $sth->fetchall_arrayref();
    $sth->finish();

    return $arrays;
}

sub get_probe_mapping_id {
    my ($dbh, $this_array_id, $assembly_name) = @_;
    my $sql = "SELECT probe_mapping_id
                FROM probe_mapping_info
                JOIN assembly USING (assembly_id)
                WHERE array_id = $this_array_id
                AND assembly_name = '$assembly_name'";
    my $probe_mapping_id = $dbh->selectrow_array($sql);
    return($probe_mapping_id);
}

sub dump_array_bed_file {
    my ($dbh, $this_array_id, $assembly_name, $work_dir) = @_;
    my $sql = "SELECT probe_mapping_id
                FROM probe_mapping_info
                JOIN assembly USING (assembly_id)
                WHERE array_id = $this_array_id
                AND assembly_name = '$assembly_name'";
    my $probe_mapping_id = $dbh->selectrow_array($sql);

    $sql = "SELECT probe_id, chr_name, position FROM probe_mapping WHERE probe_mapping_id = $probe_mapping_id";
    my $probes = $dbh->selectall_arrayref($sql);
    
    my $bed_file = "$work_dir/array_${this_array_id}.bed";
    open(BED, ">$bed_file") or die "Cannot open temporary BED file <$bed_file>\n";
    foreach my $this_probe (sort _sort_probes @$probes) {
        print BED join("\t", $this_probe->[1], $this_probe->[2], $this_probe->[2]+1, $this_probe->[0]), "\n";
    }
    close(BED);

    return $bed_file;
}

sub _sort_probes {
    my $chr_a = $a->[1];
    my $chr_b = $b->[1];
    my $loc_a = $a->[2];
    my $loc_b = $b->[2];
    $chr_a =~ s/chr//;
    $chr_b =~ s/chr//;
    if ($chr_a eq $chr_b) {
        return $loc_a <=> $loc_b;
    } elsif ($chr_a =~ /^\d/ and $chr_b =~ /^\d/) {
        return $chr_a <=> $chr_b;
    } elsif ($chr_a =~ /^\d/) {
        return -1;
    } elsif ($chr_b =~ /^\d/) {
        return 1;
    } else {
        return $chr_a cmp $chr_b;
    }
}

exit();






