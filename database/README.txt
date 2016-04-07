This directory contains the scripts and files necessary to rebuild the eFORGE database from scratch.

In summary, the steps are:
- Create an empty DB
- Load the information about the arrays
- Load the dataset with DHS, Histone peaks, etc
- Move the DB to its final destination

1. CREATE AN EMPTY DATABASE

Please refer to the help of init_db.pl for information about the different options

rm eforge_1.2.db
perl init_db.pl --db_name eforge_1.2.db

2. LOAD THE ARRAYS

Please refer to the help of load_450k_array.pl for information about the different options

perl load_450k_array.pl --work_dir input --db_name eforge_1.2.db

3. LOAD THE DATASETS

Please refer to the help of load_dataset.pl for information about the different options

perl load_dataset.pl --db_name eforge_1.2.db --tag erc --name 'Roadmap Epigenomics (2012 data) - DHS' --decode_file erc.decode --work_dir input/erc/
perl load_dataset.pl --db_name eforge_1.2.db --tag encode --name 'ENCODE - DHS' --decode_file encode.decode --work_dir input/encode/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-DHS --name 'Consolidated Roadmap Epigenomics - DHS' --decode_file erc2.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3K27me3 --name 'Consolidated Roadmap Epigenomics - H3K27me3' --decode_file erc2-H3K27me3.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3K36me3 --name 'Consolidated Roadmap Epigenomics - H3K36me3' --decode_file erc2-H3K36me3.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3K4me3 --name 'Consolidated Roadmap Epigenomics - H3K4me3' --decode_file erc2-H3K4me3.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3K9me3 --name 'Consolidated Roadmap Epigenomics - H3K9me3' --decode_file erc2-H3K9me3.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3K4me1 --name 'Consolidated Roadmap Epigenomics - H3K4me1' --decode_file erc2-H3K4me1.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag erc2-H3-all --name 'Consolidated Roadmap Epigenomics - All H3 marks' --decode_file erc2-H3-all.decode --work_dir input/erc2/
perl load_dataset.pl --db_name eforge_1.2.db --tag blueprint --name 'Blueprint - DHS' --decode_file blueprint.decode --work_dir input/blueprint/

4. MOVE THE DATABASE TO ITS FINAL LOCATION

mv eforge_1.2.db ..
