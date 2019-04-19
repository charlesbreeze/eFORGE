.echo on
.separator "\t"

CREATE TABLE assembly (
       assembly_id INTEGER PRIMARY KEY AUTOINCREMENT, 
       species_name UNIQUE, 
       assembly_name UNIQUE); 

INSERT OR IGNORE INTO assembly (species_name, assembly_name) VALUES ('Homo sapiens','GRCh37');

CREATE TABLE array (
       array_id INTEGER PRIMARY KEY AUTOINCREMENT,
       array_tag UNIQUE,
       array_name UNIQUE,
       species_name);

INSERT OR IGNORE INTO array (array_tag, array_name, species_name) VALUES ('850k','Illumina Human 850k array','Homo sapiens');

CREATE TABLE probe_mapping_info (
       probe_mapping_id INTEGER PRIMARY KEY AUTOINCREMENT,
       array_id INTEGER NOT NULL REFERENCES array(array_id),
       assembly_id INTEGER NOT NULL REFERENCES assembly(assembly_id),
       url,
       UNIQUE(array_id, assembly_id));

INSERT OR IGNORE INTO probe_mapping_info (array_id, assembly_id) VALUES ('1','1');

CREATE TABLE probe_mapping (
       probe_mapping_id INTEGER NOT NULL REFERENCES probe_mapping_info(probe_mapping_id),
       probe_id,
       chr_name,
       position INTEGER,
       UNIQUE (probe_mapping_id, probe_id));

.import probe_mapping.txt probe_mapping

CREATE INDEX probe_mapping_idx1 on probe_mapping (position, chr_name);

CREATE TABLE proxy_filter_info (
       array_id INTEGER NOT NULL REFERENCES array(array_id),
       description,
       UNIQUE(array_id));

INSERT OR IGNORE INTO proxy_filter_info (array_id, description) VALUES (1, 'Distance-based: 1kb');

CREATE TABLE proxy_filter (
       proxy_filter_id INTEGER NOT NULL REFERENCES proxy_filter_info(proxy_filter_id),
       probe_id,
       proxy_probes,
       UNIQUE (proxy_filter_id, probe_id));

.import proxy_filter.txt proxy_filter

CREATE TABLE probe_annotation_info (
       array_id INTEGER NOT NULL REFERENCES array(array_id),
       gene_reference_name NOT NULL,
       cgi_reference_name NOT NULL,
       url,
       UNIQUE(array_id));

INSERT OR IGNORE INTO probe_annotation_info (array_id, gene_reference_name, cgi_reference_name) VALUES (1, 'RefSeq', 'UCSC CpG Islands');

CREATE TABLE probe_annotation (
       array_id INTEGER NOT NULL REFERENCES probe_annotation_info(array_id),
       probe_id NOT NULL,
       gene_group,
       cgi_group,
       UNIQUE (array_id, probe_id));

.import probe_annotation.txt probe_annotation

CREATE TABLE dataset (
       dataset_id INTEGER PRIMARY KEY AUTOINCREMENT,
       dataset_tag UNIQUE,
       dataset_name UNIQUE,
       species_name);

INSERT OR IGNORE INTO dataset (dataset_tag, dataset_name, species_name) VALUES
       ('erc','Roadmap Epigenomics (2012 data) - DHS','Homo sapiens'),
       ('encode','ENCODE - DHS','Homo sapiens'),
       ('erc2-DHS','Consolidated Roadmap Epigenomics - DHS','Homo sapiens'),
       ('erc2-H3K27me3','Consolidated Roadmap Epigenomics - H3K27me3','Homo sapiens'),
       ('erc2-H3K36me3','Consolidated Roadmap Epigenomics - H3K36me3','Homo sapiens'),
       ('erc2-H3K4me3','Consolidated Roadmap Epigenomics - H3K4me3','Homo sapiens'),
       ('erc2-H3K9me3','Consolidated Roadmap Epigenomics - H3K9me3','Homo sapiens'),
       ('erc2-H3K4me1','Consolidated Roadmap Epigenomics - H3K4me1','Homo sapiens'),
       ('erc2-H3-all','Consolidated Roadmap Epigenomics - All H3 marks','Homo sapiens'),
       ('blueprint','Blueprint - DHS','Homo sapiens'),
       ('erc2-chromatin15state-all','Consolidated Roadmap Epigenomics - All chromatin 15-state marks','Homo sapiens');

CREATE TABLE sample (
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
       url);

.import sample.txt sample
.import erc2-chromatin15state-all-files_sample.txt sample

CREATE TABLE probe_bitstring (
       array_id INTEGER NOT NULL REFERENCES array(array_id),
       probe_id INTEGER NOT NULL,
       dataset_id INTEGER NOT NULL REFERENCES dataset(dataset_id),
       sum INTEGER,
       bit,
       UNIQUE(array_id,probe_id,dataset_id));

.import erc/bitstrings.txt probe_bitstring
.import encode/bitstrings.txt probe_bitstring
.import erc2-DHS/bitstrings.txt probe_bitstring
.import erc2-H3K27me3/bitstrings.txt probe_bitstring
.import erc2-H3K36me3/bitstrings.txt probe_bitstring
.import erc2-H3K4me3/bitstrings.txt probe_bitstring
.import erc2-H3K9me3/bitstrings.txt probe_bitstring
.import erc2-H3K4me1/bitstrings.txt probe_bitstring
.import erc2-H3-all/bitstrings.txt probe_bitstring
.import blueprint/bitstrings.txt probe_bitstring
.import erc2-chromatin15state-all/bitstrings.txt probe_bitstring
