#!/usr/bin/env python
#
# Copyright 2021 CIRAD
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
#
# -*- coding: utf-8 -*-

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, operator, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python PrepareBlastAndDrawDensity.py [options]\nThis program was developped at CIRAD/AGAP institute by Franc-Christophe BAURENS and Guillaume MARTIN\nFor bug report or any questions please contact franc-christophe.baurens@cirad.fr  or  guillaume.martin.cirad.fr" )
	# Wrapper options. 
	parser.add_option( '-f', '--fasta',		dest='fasta',							help='The multifasta file to work with')
	parser.add_option( '-r', '--ref',		dest='ref',								help='The fasta files of references to compare. Each fasta should be separated with "," and given as follows: Ref1Name:fasta1,Ref2Name:fasta2,...')
	parser.add_option( '-w', '--window',	dest='window',		default='100000',	help='The window (in bp) to calculate proportion of the fasta sequence on the reference; default 100kb')
	parser.add_option( '-G', '--group',		dest='group',		default=None,		help='A file containing all sequences names from --fasta options to use and each associated to a group. Names associated to same group are treated together. Names in fasta not present in the group file are not treated')
	parser.add_option( '-P', '--Prot',		dest='Prot',		default='n',		help='Fasta file contain protein sequences i.e. Amino acid coding, (y or n)')
	parser.add_option( '-s', '--slot',		dest='slot',		default='1',		help='The number of slots for blast analysis')
	parser.add_option( '-D', '--directory',	dest='directory',	default='PBADD',	help='The output directory which will contain output files')

	(options, args) = parser.parse_args()
	
	os.makedirs(options.directory, exist_ok=True)
	pathname = os.path.dirname(sys.argv[0])
	
	# Recording Sequences in fasta group if any
	dicoGroup = {}
	if options.group != None:
		file = open(options.group)
		for line in file:
			data = line.split()
			if data:
				SEQ, GP = data[0:2]
				if not (GP in dicoGroup):
					dicoGroup[GP] = set()
				dicoGroup[GP].add(SEQ)
	else:
		pass
	
	# Creating sub-fasta performing blast and calculating density
	sequence_dict = SeqIO.index(options.fasta, "fasta")
	if dicoGroup:
		for n in dicoGroup:
			outfile = open(options.directory+'/'+n+'.fasta','wt')
			for seq in dicoGroup[n]:
				if not (seq in sequence_dict):
					sys.exit("Oups, the program exited without finishing because the sequence "+n+" was not found in the fasta file\n")
				else:
					SeqIO.write(SeqRecord(sequence_dict[seq].seq, id = seq, description=''),outfile, "fasta")
			outfile.close()
			if pathname: 
				cmd = ["python", pathname+"/BlastAndDrawDensity.py", "-s", options.slot, "-f", options.directory+"/"+n+".fasta", "-r", options.ref, "-d", "n", "-w", options.window]
			else:
				cmd = ["python", "BlastAndDrawDensity.py", "-s", options.slot, "-f", options.directory+"/"+n+".fasta", "-r", options.ref, "-d", "n", "-w", options.window]

			print(' '.join(cmd))
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			subprocess.Popen.wait(process)
	else:
		for seq in sequence_dict:
			outfile = open(options.directory+'/'+seq+'.fasta','wt')
			SeqIO.write(SeqRecord(sequence_dict[seq].seq, id = seq, description=''),outfile, "fasta")
			outfile.close()
			if pathname: 
				cmd = ["python", pathname+"/BlastAndDrawDensity.py", "-s", options.slot, "-f", options.directory+"/"+seq+".fasta", "-r", options.ref, "-d", "n", "-w", options.window]
			else:
				cmd = ["python", "BlastAndDrawDensity.py", "-s", options.slot, "-f", options.directory+"/"+seq+".fasta", "-r", options.ref, "-d", "n", "-w", options.window]
	
			print(' '.join(cmd))
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			subprocess.Popen.wait(process)

	

if __name__ == "__main__": __main__()


