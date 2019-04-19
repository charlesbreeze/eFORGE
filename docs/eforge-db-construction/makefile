all: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed bitstrings extract eforge_2.0.db

MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed: MethylationEPIC_v-1-0_B4_tabs_hg19_final.txt
	cut -f1,3-5 $^ | awk -v OFS="\t" '{$$3--; $$4--; print $$2,$$3,$$4,$$1;}' | sort-bed - > $@

erc2-chromatin15state-all-files:
	./construct_erc2-chromatin15state-all_files.py
erc2-chromatin15state-all-files_fns.txt:
	./construct_erc2-chromatin15state-all_filelist.py > $@
erc2-chromatin15state-all-files_sample.txt:
	./construct_erc2-chromatin15state-all_samplelist.py > $@

bitstrings: erc encode erc2-DHS erc2-H3K27me3 erc2-H3K36me3 erc2-H3K4me3 erc2-H3K9me3 erc2-H3K4me1 erc2-H3-all blueprint erc2-chromatin15state-all

erc: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc.sh $@ $^ "eFORGE_files/" "UW*merge150.wgt10.zgt2.wig"

encode: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_encode.sh $@ $^ "eFORGE_files/" "*wig"

erc2-DHS: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-DHS.sh $@ $^ "eFORGE_files/" "*DNase.hotspot.fdr0.01.broad.bed"

erc2-H3K27me3: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3K27me3.sh $@ $^ "eFORGE_files/" "*H3K27me3.broadPeak"

erc2-H3K36me3: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3K36me3.sh $@ $^ "eFORGE_files/" "*H3K36me3.broadPeak" 

erc2-H3K4me3: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3K4me3.sh $@ $^ "eFORGE_files/" "*H3K4me3.broadPeak"

erc2-H3K9me3: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3K9me3.sh $@ $^ "eFORGE_files/" "*H3K9me3.broadPeak"

erc2-H3K4me1: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3K4me1.sh $@ $^ "eFORGE_files/" "*H3K4me1.broadPeak"

erc2-H3-all: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-H3-all.sh $@ $^ "eFORGE_files/" "*.broadPeak"

erc2-chromatin15state-all: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_erc2-chromatin15state-all.sh $@ $^

blueprint: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	./bitstrings_blueprint.sh $@ $^ "eFORGE_files/" "*v3*"

extract:
	gunzip -c erc/bitstrings.txt.gz                             | awk -v OFS="\t" '{ print "1", $$1,  "1", $$2, $$3; }' > erc/bitstrings.txt
	gunzip -c encode/bitstrings.txt.gz                          | awk -v OFS="\t" '{ print "1", $$1,  "2", $$2, $$3; }' > encode/bitstrings.txt
	gunzip -c erc2-DHS/bitstrings.txt.gz                        | awk -v OFS="\t" '{ print "1", $$1,  "3", $$2, $$3; }' > erc2-DHS/bitstrings.txt
	gunzip -c erc2-H3K27me3/bitstrings.txt.gz                   | awk -v OFS="\t" '{ print "1", $$1,  "4", $$2, $$3; }' > erc2-H3K27me3/bitstrings.txt
	gunzip -c erc2-H3K36me3/bitstrings.txt.gz                   | awk -v OFS="\t" '{ print "1", $$1,  "5", $$2, $$3; }' > erc2-H3K36me3/bitstrings.txt
	gunzip -c erc2-H3K4me3/bitstrings.txt.gz                    | awk -v OFS="\t" '{ print "1", $$1,  "6", $$2, $$3; }' > erc2-H3K4me3/bitstrings.txt
	gunzip -c erc2-H3K9me3/bitstrings.txt.gz                    | awk -v OFS="\t" '{ print "1", $$1,  "7", $$2, $$3; }' > erc2-H3K9me3/bitstrings.txt
	gunzip -c erc2-H3K4me1/bitstrings.txt.gz                    | awk -v OFS="\t" '{ print "1", $$1,  "8", $$2, $$3; }' > erc2-H3K4me1/bitstrings.txt
	gunzip -c erc2-H3-all/bitstrings.txt.gz                     | awk -v OFS="\t" '{ print "1", $$1,  "9", $$2, $$3; }' > erc2-H3-all/bitstrings.txt
	gunzip -c blueprint/bitstrings.txt.gz                       | awk -v OFS="\t" '{ print "1", $$1, "10", $$2, $$3; }' > blueprint/bitstrings.txt
	gunzip -c erc2-chromatin15state-all/bitstrings.txt.gz       | awk -v OFS="\t" '{ print "1", $$1, "11", $$2, $$3; }' > erc2-chromatin15state-all/bitstrings.txt

probe_mapping.txt: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	awk -v OFS="\t" '{ print "1", $$4, $$1, $$2 }' $^ > $@

proxy_filter.txt: MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	bedmap --echo --echo-map-id-uniq --range 500 --delim '\t' --multidelim '|' $^ | cut -f4- | awk '{ printf("1\t%s\t",$$1); n=split($$2,a,"|"); if(n==1){printf("NONE|");}else{for(i=1;i<=n;i++){if(a[i]!=$$1){printf("%s|",a[i]);}}} printf("\n"); }' | awk '{ print substr($$0,1,length($$0)-1); }' > $@

probe_annotation.txt: MethylationEPIC_v-1-0_B4_tabs_hg19_final.txt
	cut -f1,6,7 $^ | awk -v OFS="\t" '{ print "1", $$0; }' > $@

eforge_2.0.db: probe_mapping.txt proxy_filter.txt probe_annotation.txt
	cat eforge_2.0.sql | sqlite3 $@

clean:
	rm -f MethylationEPIC_v-1-0_B4_tabs_hg19_final.bed
	rm -rf erc
	rm -rf encode
	rm -rf erc2-DHS
	rm -rf erc2-H3K27me3
	rm -rf erc2-H3K36me3
	rm -rf erc2-H3K4me3
	rm -rf erc2-H3K9me3
	rm -rf erc2-H3K4me1
	rm -rf erc2-H3-all
	rm -rf erc2-chromatin15state-all
	rm -rf blueprint
	rm -f probe_mapping.txt
	rm -f proxy_filter.txt
	rm -f probe_annotation.txt
	rm -f eforge_2.0.db
