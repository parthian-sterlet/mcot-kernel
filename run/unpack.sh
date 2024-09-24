chmod a+x command_line_one
chmod a+x command_line_many
chmod a+x command_line_pro
chmod a+x command_line_library
chmod a+x command_line_denovo
chmod a+x demc.pl

cd ../partners
cat h12core_hg38.binary.tar.gz.part* > h12core_hg38.binary.tar.gz
cat h12core_mm10.binary.tar.gz.part* > h12core_mm10.binary.tar.gz
cat dapseq_at10.binary.tar.gz.part* > dapseq_at10.binary.tar.gz
tar -xvzf h12core_hg38.binary.tar.gz
tar -xvzf h12core_mm10.binary.tar.gz
tar -xvzf dapseq_at10.binary.tar.gz

cd ../genomes/at
tar -xzvf ups1500_at10.seq.tar.gz
cd ../dm
tar -xzvf ups1500_dm6.seq.tar.gz
cd ../mm
tar -xzvf ups2kb_mm10.seq.tar.gz
cd ../hs
tar -xzvf ups2kb_hg38.seq.tar.gz
cd ..


