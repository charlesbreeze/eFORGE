#! /usr/bin/env perl
use strict;
use warnings;

use DBI;
use Getopt::Long;

my $db_dir = ".";
my $db_name = "eforge_1.2.db";


my $help;

my $desc = qq{init_db.pl [--db_name $db_name] [--db_dir $db_dir]

where:
 --db_name is the name of the SQLite file [def: $db_name]
 --db_dir is the location of the SQLite file [def: $db_dir]

};

GetOptions(
    "help" => \$help,
    "db_name=s" => \$db_name,
    "db_dir=s" => \$db_dir,
    );

if ($help) {
    print $desc;
    exit(0);
}

my $dsn = "dbi:SQLite:dbname=$db_dir/$db_name";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;


$dbh->do("CREATE TABLE IF NOT EXISTS assembly (
    assembly_id INTEGER PRIMARY KEY AUTOINCREMENT,
    species_name,
    assembly_name,
    UNIQUE (species_name, assembly_name))");


$dbh->do("CREATE TABLE IF NOT EXISTS array (
    array_id INTEGER PRIMARY KEY AUTOINCREMENT,
    array_tag UNIQUE,
    array_name UNIQUE,
    species_name)");
$dbh->do("CREATE TABLE IF NOT EXISTS probe_mapping_info (
    probe_mapping_id INTEGER PRIMARY KEY AUTOINCREMENT,
    array_id INTEGER NOT NULL REFERENCES array(array_id),
    assembly_id INTEGER NOT NULL,
    url,
    UNIQUE (array_id, assembly_id))");
$dbh->do("CREATE TABLE IF NOT EXISTS probe_mapping (
    probe_mapping_id INTEGER NOT NULL REFERENCES probe_mapping_info(probe_mapping_id),
    probe_id,
    chr_name,
    position INTEGER,
    UNIQUE (probe_mapping_id, probe_id))");
$dbh->do("CREATE INDEX probe_mapping_idx1 on probe_mapping (position, chr_name)");


$dbh->do("CREATE TABLE IF NOT EXISTS proxy_filter_info (
    array_id INTEGER NOT NULL REFERENCES array(array_id),
    description,
    UNIQUE(array_id))");
$dbh->do("CREATE TABLE IF NOT EXISTS proxy_filter (
    array_id INTEGER NOT NULL REFERENCES array(array_id),
    probe_id,
    proxy_probes,
    UNIQUE (array_id, probe_id))");


$dbh->do("CREATE TABLE IF NOT EXISTS probe_annotation_info (
    array_id INTEGER NOT NULL REFERENCES array(array_id),
    gene_reference_name NOT NULL,
    cgi_reference_name NOT NULL,
    url,
    UNIQUE(array_id))");
$dbh->do("CREATE TABLE IF NOT EXISTS probe_annotation (
    array_id INTEGER NOT NULL REFERENCES probe_annotation_info(array_id),
    probe_id NOT NULL,
    gene_group,
    cgi_group,
    UNIQUE (array_id, probe_id))");


$dbh->do("CREATE TABLE IF NOT EXISTS dataset (
    dataset_id INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_tag UNIQUE,
    dataset_name UNIQUE,
    species_name)");
$dbh->do("CREATE TABLE IF NOT EXISTS sample (
    dataset_id INTEGER REFERENCES dataset(dataset_id),
    sample_order INTEGER,
    file,
    lab,
    datatype,
    cell,
    tissue,
    shortcell,
    individual,
    acc,
    url)");
$dbh->do("CREATE TABLE IF NOT EXISTS probe_bitstring (
    array_id INTEGER NOT NULL REFERENCES array(array_id),
    probe_id INTEGER NOT NULL,
    dataset_id INTEGER NOT NULL REFERENCES dataset(dataset_id),
    sum INTEGER,
    bit,
    UNIQUE(array_id,probe_id,dataset_id))");
