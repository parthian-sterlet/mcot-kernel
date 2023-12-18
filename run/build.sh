chmod a+x command_line_one
chmod a+x command_line_many
chmod a+x command_line_pro

g++ -o mcot_anchor_pro.exe ../src/anchor_pro/mcot_anchor_pro.cpp
g++ -o mcot_anchor.exe ../src/anchor_vs_one_partner/mcot_anchor.cpp
g++ -o mcot.exe ../src/anchor_vs_many_partners/mcot.cpp

chmod a+x mcot_anchor_pro.exe
chmod a+x mcot_anchor.exe
chmod a+x mcot.exe
