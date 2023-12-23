chmod a+x command_line_one
chmod a+x command_line_many
chmod a+x command_line_pro
chmod a+x command_line_library

cd ../partners
cat h12core_hg38.binary.tar.gz.part* > h12core_hg38.binary.tar.gz
cat h12core_mm10.binary.tar.gz.part* > h12core_mm10.binary.tar.gz
cat dapseq_at10.binary.tar.gz.part* > dapseq_at10.binary.tar.gz
tar -xvzf h12core_hg38.binary.tar.gz
tar -xvzf h12core_mm10.binary.tar.gz
tar -xvzf dapseq_at10.binary.tar.gz

cd ../genomes/at
tar -xvzf ups1500_at10.seq.tar.gz
cd ../dm
tar -xvzf ups1500_dm6.seq.tar.gz
cd ../mm
tar -xvzf ups2kb_mm10.seq.tar.gz
cd ../hs
tar -xvzf ups2kb_hg38.seq.tar.gz
cd ..


