#!/usr/bin/env python
#
#  Copyright 2021 CIRAD
# Guillaume Martin Franc-Christophe Baurens
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
# -*- coding: utf-8 -*-

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, operator, time, configparser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Drawing function
def draw_plot_ref_per_chr(X, Y, OUT, TITLE, FILLUNDER, YLIM, COLOR, NEGATIVE): 
	
	fig = plt.figure(figsize=(15, 5))
	ax = plt.subplot2grid((1,1),(0,0))
	
	if NEGATIVE == 'y':
		ax.plot(X, Y, linewidth=2, color=(0.75, 0.31, 0.30))
		ax.plot(X, Y*-1, linewidth=2, color=(0.75, 0.31, 0.30))
		if FILLUNDER == 'y':
			ax.fill_between(X,Y*-1, Y, color=(0.75, 0.31, 0.30))
		ax.set_xlim(0, X[-1])
		if YLIM != None:
			ax.set_ylim(float(YLIM*-1), float(YLIM))
		else:
			ax.set_ylim(max(Y)*-1, max(Y))
	else:
		ax.plot(X, Y, linewidth=2, color=(0.75, 0.31, 0.30))
		if FILLUNDER == 'y':
			ax.fill_between(X,0, Y, color=(0.75, 0.31, 0.30))
		ax.set_xlim(0, X[-1])
		if YLIM != None:
			ax.set_ylim(0, float(YLIM))
		else:
			ax.set_ylim(0, max(Y))
	ax.set_title(TITLE, fontweight='bold', position=(0.5, 1))
	
	
	fig.savefig(OUT)
	plt.close(fig)

def draw_plot_ref_per_ref(DICO, OUT, FILLUNDER, YLIM, DICOCHR, COLOR, NEGATIVE):
	
	# getting chromosomes list
	chr_list = sorted(list(DICO.keys()))
	
	# Calculating chromosome number
	NB = len(chr_list)
	
	# Calculating max chromosome length
	XLIM = 0
	for chr in chr_list:
		XLIM = max(XLIM, DICOCHR[chr])
	
	# Calculating max Y value si on l a pas fixe
	if YLIM == None:
		YLIM = 0
		for chr in chr_list:
			YLIM = max(YLIM, max(DICO[chr][1]))
	else:
		YLIM = int(YLIM)
	
	# Drawing figures
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# Drawing graph
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=14, rowspan=1)
		if NEGATIVE == 'y':
			ax.set_ylim(YLIM*-1, YLIM)
		else:
			ax.set_ylim(0, YLIM)
		ax.set_xlim(0, XLIM)
		ax.tick_params(axis='both', which='major', labelsize=8)
		if NEGATIVE == 'y':
			ax.plot(DICO[chr][0], np.array(DICO[chr][1])*-1, linewidth=2, color=COLOR)
		else:
			pass
		ax.plot(DICO[chr][0], DICO[chr][1], linewidth=2, color=COLOR)
		if FILLUNDER == 'y':
			if NEGATIVE == 'y':
				ax.fill_between(DICO[chr][0],np.array(DICO[chr][1])*-1, DICO[chr][1], color=COLOR)
			else:
				ax.fill_between(DICO[chr][0],0, DICO[chr][1], color=COLOR)
		ax.axvline(x=DICOCHR[chr],linewidth=2, color='black')
		POSSPAN += 1
	
	# Drawing chromosome name
	POSSPAN = 0
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,14), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', ha='left', fontweight='bold')
		
		POSSPAN += 1
	
	fig.savefig(OUT)
	plt.close(fig)

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python PrepareBlastAndDrawDensity.py [options]\nThis program was developped at CIRAD/AGAP institute by Franc-Christophe BAURENS and Guillaume MARTIN\nFor bug report or any questions please contact franc-christophe.baurens@cirad.fr  or  guillaume.martin.cirad.fr" )
	# Wrapper options. 
	parser.add_option( '-f', '--fasta',		dest='fasta',							help='The fasta file to work with')
	parser.add_option( '-r', '--ref',		dest='ref',								help='The fasta files of references to compare. Each fasta should be separated with "," and given as follows: Ref1Name:fasta1,Ref2Name:fasta2,...')
	parser.add_option( '-w', '--window',	dest='window',		default='100000',	help='The window to calculate proportion of the fasta sequence on the reference; default 100kb')
	parser.add_option( '-s', '--slot',		dest='slot',		default='1',		help='The number of slots for blast analysis')
	parser.add_option( '-F', '--FillUnder',	dest='FillUnder',	default='n',		help='Fill the curve under, (y or n)')
	parser.add_option( '-c', '--chrToRm',	dest='chrToRm',		default='',			help='List of chromosome to exclude from drawing separated with ","')
	parser.add_option( '-C', '--color',		dest='color',		default=None,		help='A color file in tabulated format (in RGB). Col1: RefxName, Col2: Red, Col3: Green, Col4: Blue')
	parser.add_option( '-y', '--Ylim',		dest='Ylim',		default=None,		help='Set Y limit for all plots. Recommended to compare genomes with the same graphic scale --> if omitted each graph may have different Yscale (adjusted from each chromosome)')
	parser.add_option( '-P', '--Prot',		dest='Prot',		default='n',		help='Fasta file contain protein sequences i.e. Amino acid coding, (y or n)')
	parser.add_option( '-d', '--draw',		dest='draw',		default='i',		help='Drawing output. possible options: One figure per genome (g), one figure per chromosomes per genome (i) , No drawing (n), options can be combined. n option overwrite others')
	parser.add_option( '-n', '--negative',	dest='negative',	default='n',		help='Draw also negative curve: y or n')
	parser.add_option( '-g', '--graph',		dest='graph',		default='svg',		help='graphic output. possible options: pdf, png, svg')


	(options, args) = parser.parse_args()
	
	ScriptPath = os.path.dirname(sys.argv[0])
	loca_programs = configparser.RawConfigParser()
	if ScriptPath:
		loca_programs.read(ScriptPath+'/locaprogram.conf')
	else:
		loca_programs.read('locaprogram.conf')
	
	ChrToExclude=options.chrToRm.split(',')
	
	WIN = int(options.window)
	SLOT = int(options.slot)
	print("sequence to compare on references = ",options.fasta)
	HEAD = []
	HEAD = options.fasta.split('/')
	fastaname = HEAD[-1]
	HEAD = fastaname.split ('.')
	fastaname = HEAD[0]
	print ("fastaname = "+fastaname)
	file = open(options.fasta,"r")
	count=0
	for line in file:
		data = line.split()
		if ">" in data:
			count+=1
	if count>1:
		print("Sequences in the fasta file are treated together and not reported individually in the graphics")
	
	# Get info in --ref argument
	dicoRef = {}
	REFS = options.ref.split(',')
	i = 0
	for n in REFS:
		i+=1
		Fname = n.split(':')
		dicoRef[Fname[0]] = Fname[1]
		print("Reference sequence"+str(i)+" = ", Fname[0], Fname[1])
	
	# get colors if necessary
	DicoCol = {}
	if options.color == None:
		for ref in dicoRef:
			DicoCol[ref] = (238/255, 16/255, 16/255)
	else:
		file = open(options.color)
		for line in file:
			data = line.split()
			if data:
				DicoCol[data[0]] = (int(data[1])/255, int(data[2])/255, int(data[3])/255)
		file.close()
		for ref in dicoRef:
			if not (ref in DicoCol):
				print(ref+' was not found in color file, an arbitrary color was attributed')
				DicoCol[ref] = (238/255, 16/255, 16/255)
	
	# Running blasts
	for ref in dicoRef:
		if options.Prot == 'y' or options.Prot == 'yes' :
			print('Running protein blast on '+ref)
			cmd = loca_programs.get('Programs','tblastn').split()+["-query", options.fasta, "-db", dicoRef[ref],  "-num_threads", options.slot,"-outfmt", "7", "-evalue", "1E-100", "-out", fastaname+"_"+ref+"_blast.tab"]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			subprocess.Popen.wait(process)
		else:
			print('Running nucleotid blast on '+ref)
			cmd = loca_programs.get('Programs','blastn').split()+["-query", options.fasta, "-db", dicoRef[ref],  "-num_threads", options.slot,"-outfmt", "7", "-dust", "no", "-out", fastaname+"_"+ref+"_blast.tab"]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			subprocess.Popen.wait(process)
	
	# Drawing figures
	for ref in dicoRef:
		print('Calculating density (and drawing) for '+ref)
		outfile = open(fastaname+"_"+ref+"_density.tab","w")
		# Getting ref sequence statistics
		dico_chr = {}
		dico_chr_len = {}
		sequence_dict = SeqIO.index(dicoRef[ref], "fasta")
		for n in sequence_dict:
			dico_chr[n] = [0]*len(str(sequence_dict[n].seq))
			dico_chr_len[n] = len(str(sequence_dict[n].seq))
		del sequence_dict
		
		# Filling the hash ref1
		file = open(fastaname+"_"+ref+"_blast.tab","r")
		for line in file:
			data = line.split()
			if data:
				if "#" not in data: 
					chr = data[1]
					if not (chr in dico_chr):
						sys.exit('The system exited without finishing because sequence '+chr+' has not been found in the fasta')
					if int(data[8]) < int(data[9]):
						debut = int(data[8])-1
						fin = int(data[9])-1
					else:
						debut = int(data[9])-1
						fin = int(data[8])-1
					while debut <= fin:
						dico_chr[chr][debut] = 1
						debut += 1
		
		# Calculating window proportion
		if 'g' in options.draw and not('n' in options.draw):
			dicodraw = {}
		for chr in dico_chr:
			debut = 0
			fin = WIN
			chr_fin = len(dico_chr[chr])-1
			list_pos = [1]
			list_cov = []
			while debut < chr_fin:
				temp_fin = fin
				if chr_fin < fin:
					temp_fin = chr_fin + 1
				sublist = dico_chr[chr][debut:temp_fin]
				cov = (sum(sublist)/float(len(sublist)))*WIN
				outfile.write('\t'.join([chr, str(debut+1), str(temp_fin), str(cov)])+'\n')
				if len(list_cov) == 0:
					list_cov.append(cov)
				list_cov.append(cov)
				list_pos.append(temp_fin)
				debut += WIN
				fin += WIN
				sys.stdout.flush()
			# Draw figures
			if not (chr in ChrToExclude):
				if 'g' in options.draw:
					dicodraw[chr]=[[],[]]
					dicodraw[chr][0] = list_pos[:]
					dicodraw[chr][1] = list_cov[:]
				if 'i' in options.draw and not('n' in options.draw):
					draw_plot_ref_per_chr(np.array(list_pos), np.array(list_cov), ref+"_"+chr+"_"+fastaname+"."+options.graph, chr, options.FillUnder, int(options.Ylim), DicoCol[ref], options.negative)
		if 'g' in options.draw and not('n' in options.draw):
			draw_plot_ref_per_ref(dicodraw, ref+"_"+fastaname+"."+options.graph, options.FillUnder, int(options.Ylim), dico_chr_len, DicoCol[ref], options.negative)
		outfile.close()


if __name__ == "__main__": __main__()


