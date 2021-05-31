#Figures of the manuscript were generated using sequence data, and configuration files provided in the DATA folder and using scripts provided in the Script folder.

#1- Requirements:
#These scripts were tested with python 3.4.3 and required the following modules : 
optparse, os, shutil, subprocess, sys, tempfile, fileinput, operator, time, optparse, math, configparser, numpy,matplotlib, matplotlib.pyplot,and BIO python

#To Draw the figures you need the python scripts located in the Script folder : 
BlastAndDrawDensity.py
DrawStackedDensity.py
PrepareBlastAndDrawDensity.py

# you need also to update the path to your installation of ncbi-blast program (i.e. blastn and tblatn programs) in the locaprogram.conf file (after the "=" symbol): 
#All sequences necessary to draw figures are provided in the Data folder.

#2-To draw the figures on V2 and V4 pseudochromosomes : please send the following commands from the TEST_TOOL folder.

python ../Script/PrepareBlastAndDrawDensity.py -f ../DATA/Final_Figure_TE.fasta -r PHG_V4:../DATA/DH-Pahang.4.3.fasta,PHG_V2:../DATA/musa_acuminata_v2_pseudochromosome.fna -G ../DATA/Group_Final_Figure_TE.tab

#Then send the following commands to draw the figures
python ../Script/DrawStackedDensity.py -f NANICA:NANICA_PHG_V4_density.tab,CL18:CL18_PHG_V4_density.tab,CL33:CL33_PHG_V4_density.tab,5S:5S_PHG_V4_density.tab,45S:45S_PHG_V4_density.tab,MAX:MAX_PHG_V4_density.tab -F y -d gs -n y -g svg -p STACK-Final_Figure_TE_Musa_V4 -c chloro,putative_mito_scaff,putative_mitochondrion,mito1,mito2,mito3,mito4,mito5,mito6,mito7,mito8,mito9,mito10,mito11,mito12 -C  ../DATA/Revised_Figure_TE_Color.conf

python ../Script/DrawStackedDensity.py -f NANICA:NANICA_PHG_V2_density.tab,CL18:CL18_PHG_V2_density.tab,CL33:CL33_PHG_V2_density.tab,5S:5S_PHG_V2_density.tab,45S:45S_PHG_V2_density.tab,MAX:MAX_PHG_V2_density.tab -F y -d gs -n y -g svg -p STACK-Final_Figure_TE_Musa_V2 -c chloro,putative_mito_scaff,putative_mitochondrion,mito1,mito2,mito3,mito4,mito5,mito6,mito7,mito8,mito9,mito10,mito11,mito12 -C ../DATA/Revised_Figure_TE_Color.conf

#two files are generated : STACK-Final_Figure_TE_Musa_V2.svg and STACK-Final_Figure_TE_Musa_V4.svg
#Additonnal artworks from these svg files to generate the final figures have been performed by Franc-Christophe BAURENS using Inkscape v0.91 https://inkscape.org/


#3-To draw the figures 3 : please send the following commands from the TEST_TOOL folder.

python ../Script/PrepareBlastAndDrawDensity.py -f ../DATA/Final_Figure_TE.fasta -r 1Mb_rDNA_Chr01:../DATA/chr01_27061466_28061465.fasta -w 2000 -G ../DATA/Group_Revised_new_colors.tab

#Then send the following commands to draw the figure
python ../Script/DrawStackedDensity.py -f NANICA:NANICA_1Mb_rDNA_Chr01_density.tab,CRM:CRM_1Mb_rDNA_Chr01_density.tab,CL18:CL18_1Mb_rDNA_Chr01_density.tab,CL33:CL33_1Mb_rDNA_Chr01_density.tab,5S:5S_1Mb_rDNA_Chr01_density.tab,45S:45S_1Mb_rDNA_Chr01_density.tab,MAX:MAX_1Mb_rDNA_Chr01_density.tab -F y -d gs -n y -g svg -p STACK_REVISED_Final_TE_Musa_1Mb_rDNA_Chr01 -c chloro,putative_mito_scaff,putative_mitochondrion,mito1,mito2,mito3,mito4,mito5,mito6,mito7,mito8,mito9,mito10,mito11,mito12 -C ../DATA/Revised_Figure_TE_Color.conf

#Dotplots were generated with Gepard 1.40 on a windows computer https://cube.univie.ac.at/gepard
#1Mb sequence from chromosome 01 centronmeric region : chr01_27061466_28061465.fasta
#30kb sequence from chromosome 01 with rDNA Repeats : 30kb_rDNA.fasta 
#Additonnal artworks from svg and Gepard outfiles to generate final figures have been performed by Franc-Christophe BAURENS using Inkscape v0.91 https://inkscape.org/

