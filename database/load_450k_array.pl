#! /usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use File::Spec qw(splitpath);
use DBI;

my $db_dir = ".";
my $db_name = "eforge_1.1.db";
my $array_tag = "450k";
my $array_name = "Illumina Human 450k array";
my $species = "Homo sapiens";
my $proxy_threshold = 1000;
my $illumina450k_csv_file = 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv';
my $work_dir = ".";
my $bedtools = "bedtools";

my $help;

my $desc = qq{load_450k_array.pl [options]

DESCRIPTION:

This script loads the probe locations and annotations of the Illumina Human 450k methylation array
into the eFORGE database.

Note that you *must* load the arrays before loading the datasets. If you want to include new arrays
at a later date, you will have to reload all the datasets again (i.e. you will have to re-start from
scratch).

Optional arguments:
 --db_name <name>
    is the name of the SQLite file [def: $db_name]
 --db_dir <path>
    is the location of the SQLite file [def: $db_dir]
 --work_dir <path>
    is the location where temporary files will be downloaded/created [def: $work_dir]
 --bedtools <name>
    is the name of the bedtools executable you want to use [def: $bedtools]
};

GetOptions(
    "help" => \$help,
    "db_name=s" => \$db_name,
    "db_dir=s" => \$db_dir,
    "work_dir=s" => \$work_dir,
    "bedtools=s" => \$bedtools,
    );

my $dsn = "dbi:SQLite:dbname=$db_dir/$db_name";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;

my ($vol, $apth, $file) = File::Spec->splitpath($illumina450k_csv_file);

if (!-e "$work_dir/$file") {
    download_url($illumina450k_csv_file, $work_dir);
}

my $array_info = parse_450k_file("$work_dir/$file");

load_array($dbh, $db_name, $species, $array_tag, $array_name, $array_info);

load_proxy_filter($dbh, $db_name, $array_tag, "GRCh37", $work_dir, $proxy_threshold);

exit();


sub download_url {
    my ($url, $work_dir) = @_;
    
    print "Downloading $url...\n";
    system("wget -q -N -P $work_dir $url");
}


sub parse_450k_file {
    my ($illumina450k_csv_file) = @_;
    my $array;
    
    open(CSV, $illumina450k_csv_file) or die;
    my $annotation;
    my @gene_annotations = ("TSS200", "TSS1500", "1stExon", "Body", "3'UTR", "5'UTR", "IGR");
    my $cpg_annotations = {
        'N_Shelf' => "Shelf_Shore",
        'S_Shelf' => "Shelf_Shore",
        'N_Shore' => "Shelf_Shore",
        'S_Shore' => "Shelf_Shore",
        'Island' => "Island"};
    while (<CSV>) {
        chomp;
        my @data = split(",", $_);
        my $probe_id = $data[0];
        next if ($probe_id !~ /^cg/ and $probe_id !~ /^ch\./);
        my $probe_chr37 = $data[11];
        my $probe_loc37 = $data[12];
        my $probe_chr36 = $data[14];
        my $probe_loc36 = $data[15];
        my $this_gene_annotation_arrayStr = $data[23];
        my $this_cpg_annotation = $data[25]?$cpg_annotations->{$data[25]}:"NA";

        my $this_gene_annotation;
        if (!defined($this_gene_annotation_arrayStr) or $this_gene_annotation_arrayStr eq "") {
            $this_gene_annotation = "IGR";
        } else {
            my @this_gene_annotations = split(";", $this_gene_annotation_arrayStr);
            foreach my $this_a (@this_gene_annotations) {
                if (grep {$_ eq $this_a} @gene_annotations) {
                    $this_gene_annotation = $this_a;
                    last;
                }
            }
        }

        $annotation->{$this_gene_annotation."-".$this_cpg_annotation}++;

        $array->{$probe_id} = [$probe_chr36, $probe_loc36, $probe_chr37, $probe_loc37, $this_gene_annotation, $this_cpg_annotation];
    }
#     foreach my $this_a (keys %$annotation) { #@gene_annotations) {
#         print $annotation->{$this_a}, "\t", $this_a, "\n";
#     }
    
    return $array;
}


sub load_array_1 {
    my ($dbh, $code_name, $array_info) = @_;

    my $table_name = $code_name;
    $table_name =~ s/\W//g;
    $table_name = "array_$table_name";
    $dbh->do("CREATE TABLE IF NOT EXISTS $table_name (probe_id, location, gene_group, cgi_group)");

    my $sth = $dbh->prepare("INSERT INTO array_$code_name VALUES (?, ?, ?, ?)");
    foreach my $this_probe_id (sort keys $array_info) {
        my ($chr, $loc, $gene_group, $cgi_group) = @{$array_info->{$this_probe_id}};
        my $location = "$chr:$loc-$loc";
        $sth->execute($this_probe_id, $location, $gene_group, $cgi_group);
    }
    $sth->finish();

    return $table_name;
}


sub load_array {
    my ($dbh, $db_name, $species_name, $array_tag, $array_name, $array_info) = @_;

    my $sth;
    $sth = $dbh->prepare("INSERT OR IGNORE INTO array (array_tag, array_name, species_name) VALUES (?, ?, ?)");
    $sth->execute($array_tag, $array_name, $species_name);
    my $array_id = $dbh->last_insert_id("", "", "", "");
    $sth->finish();
    if ($array_id == 0) {
        $array_id = $dbh->selectrow_array("SELECT array_id FROM array WHERE array_name = '$array_name' AND species_name = 'Homo sapiens'");
    }

    $sth = $dbh->prepare("INSERT OR IGNORE INTO assembly (species_name, assembly_name) VALUES (?, ?)");
    $sth->execute($species_name, "GRCh37");
    $sth->execute($species_name, "NCBI36");
    $sth->finish();
    my $human36_assembly_id = $dbh->selectrow_array("SELECT assembly_id FROM assembly WHERE species_name = 'Homo sapiens' AND assembly_name = 'NCBI36'");
    my $human37_assembly_id = $dbh->selectrow_array("SELECT assembly_id FROM assembly WHERE species_name = 'Homo sapiens' AND assembly_name = 'GRCh37'");

    $sth = $dbh->prepare("INSERT OR IGNORE INTO probe_mapping_info (array_id, assembly_id, url) VALUES (?, ?, ?)");
    $sth->execute($array_id, $human36_assembly_id, "");
    my $probe_mapping_human36_id = $dbh->last_insert_id("", "", "", "");
    $sth->execute($array_id, $human37_assembly_id, "");
    my $probe_mapping_human37_id = $dbh->last_insert_id("", "", "", "");
    $sth->finish();

    $sth = $dbh->prepare("INSERT OR IGNORE INTO probe_annotation_info (array_id, gene_reference_name, cgi_reference_name, url) VALUES (?, ?, ?, ?)");
    $sth->execute($array_id, "RefSeq", "UCSC CpG Islands", "");
#     my $probe_annotation_id = $dbh->last_insert_id("", "", "", "");
    $sth->finish();

    open(CSV1, ">probe_mapping.csv") or die;
    open(CSV2, ">probe_annotation.csv") or die;
#     my $sth1 = $dbh->prepare("INSERT OR IGNORE INTO probe_mapping (probe_mapping_id, probe_id, chr_name, position) VALUES (?, ?, ?, ?)");
#     my $sth2 = $dbh->prepare("INSERT OR IGNORE INTO probe_annotation (probe_annotation_id, probe_id, gene_group, cgi_group) VALUES (?, ?, ?, ?)");
    foreach my $this_probe_id (sort keys $array_info) {
        my ($chr36, $loc36, $chr37, $loc37, $gene_group, $cgi_group) = @{$array_info->{$this_probe_id}};
        print CSV1 join(",", $probe_mapping_human36_id, $this_probe_id, "chr$chr36", $loc36), "\n";
        print CSV1 join(",", $probe_mapping_human37_id, $this_probe_id, "chr$chr37", $loc37), "\n";
#         print CSV2 join(",", $probe_annotation_id, $array_id, $this_probe_id, $gene_group, $cgi_group), "\n";
        print CSV2 join(",", $array_id, $this_probe_id, $gene_group, $cgi_group), "\n";
#         $sth1->execute($probe_mapping_human36_id, $this_probe_id, "chr$chr36", $loc36) if ($probe_mapping_human36_id);
#         $sth1->execute($probe_mapping_human37_id, $this_probe_id, "chr$chr37", $loc37) if ($probe_mapping_human37_id);
#         $sth2->execute($probe_annotation_id, $this_probe_id, $gene_group, $cgi_group) if ($probe_annotation_id);
    }
    close(CSV1);
    close(CSV2);
    system("echo '.mode csv
.import probe_mapping.csv probe_mapping
.import probe_annotation.csv probe_annotation' | sqlite3 $db_name");
#     $sth1->finish();
#     $sth2->finish();

}

sub load_proxy_filter {
    my ($dbh, $db_name, $array_tag, $assembly_name, $work_dir, $distance_threshold) = @_;
    
    my $array_id = $dbh->selectrow_arrayref("SELECT array_id FROM array WHERE array_tag = '$array_tag'")->[0];
    my $bed_file = dump_array_bed_file($dbh, $array_id, $assembly_name, $work_dir);
    my $this_output_bed_file = "$work_dir/array_${array_id}.proxy.bed";
    my $runstr = "$bedtools window -w $distance_threshold -a $bed_file -b $bed_file > $this_output_bed_file";
    system($runstr) == 0 or die "Error while running bedtools: $?";
    
    my $sth = $dbh->prepare("INSERT OR IGNORE INTO proxy_filter_info (array_id, description) VALUES (?, ?)");
    $sth->execute($array_id, "Distance-based: 1kb");

    open(BED, $this_output_bed_file) or die;
    my $all_probes;
    my $mapping_probes;
    while (<BED>) {
        chomp;
        my ($chr1, $start1, $end1, $probe1, $chr2, $start2, $end2, $probe2) = split("\t", $_);
        $all_probes->{$probe1} = 1;
        if ($probe1 ne $probe2) {
            $mapping_probes->{$probe1}->{$probe2} = 1;
        }
    }
    close(BED);
    open(CSV, ">proxy_filter.csv") or die;
    foreach my $probe (sort keys %$all_probes) {
        if (defined($mapping_probes->{$probe})) {
            print CSV "$array_id,$probe,", join("|", sort keys %{$mapping_probes->{$probe}}), "\n";
        } else {
            print CSV "$array_id,$probe,NONE\n";
        }
    }
    close(CSV);
    system("echo '.mode csv
.import proxy_filter.csv proxy_filter' | sqlite3 $db_name");

}

sub dump_array_bed_file {
    my ($dbh, $this_array_id, $assembly_name, $work_dir) = @_;
    my $sql = "SELECT probe_mapping_id
                FROM probe_mapping_info
                JOIN assembly USING (assembly_id)
                WHERE array_id = $this_array_id
                AND assembly_name = '$assembly_name'";
    my $probe_mapping_id = $dbh->selectrow_array($sql);
    if (!$probe_mapping_id) {
        die "Cannot find a mapping for array $this_array_id and assembly $assembly_name\n";
    }

    $sql = "SELECT probe_id, chr_name, position FROM probe_mapping WHERE probe_mapping_id = $probe_mapping_id";
    my $probes = $dbh->selectall_arrayref($sql);
    
    my $bed_file = "$work_dir/array_${this_array_id}.bed";
    open(BED, ">$bed_file") or die;
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
