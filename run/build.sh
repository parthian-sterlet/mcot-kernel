chmod a+x command_line_one
chmod a+x command_line_many
chmod a+x command_line_pro

g++ -o mcot_anchor_pro.exe ../src/anchor_pro/mcot_anchor_pro.cpp
g++ -o mcot_anchor.exe ../src/anchor_vs_one_partner/mcot_anchor.cpp
g++ -o mcot.exe ../src/anchor_vs_many_partners/mcot.cpp

chmod a+x mcot_anchor_pro.exe
chmod a+x mcot_anchor.exe
chmod a+x mcot.exe

cd ../include
cat h12core_hg38.binary.tar.gz.part* > h12core_hg38.binary.tar.gz
cat h12core_mm10.binary.tar.gz.part* > h12core_mm10.binary.tar.gz
cat dapseq_at10.binary.tar.gz.part* > dapseq_at10.binary.tar.gz
tar -xvzf h12core_hg38.binary.tar.gz
tar -xvzf h12core_mm10.binary.tar.gz
tar -xvzf dapseq_at10.binary.tar.gz

cd ../genomes/at
tar -xvzf ups1500.seq.tar.gz
cd ../mm
tar -xvzf ups2kb.seq.tar.gz
cd ../hs
tar -xvzf ups2kb.seq.tar.gz
cd ..
cd ..


