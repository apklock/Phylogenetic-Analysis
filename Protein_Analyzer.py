#! /usr/bin/env python

"""This is my protein analysis program which takes an input of amino acids and gives you
back several key indices and data on your prospective protein or peptide.

This program needs BioPython installed for some of the modules to work correctly.
Also, you will have to have ProteinParamData and IsoelectricPoint compiled on your system.

To call up the program, simply run Protein_Analyzer.
Then specify a protein/peptide sequence in the following manner:
	
	Prot = Protein_Analyzer("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA")
	
	This is an example peptide only...it is the A Chain of the Amyloid Beta protein, which is
	implicated in Alzheimer's disease due to Amyloid plaque creation.  It is a peptide that I
	have done considerable research on and it is very dear to me...in a scientific way, obviously...
 
 Once you have specified your protein/peptide you can call up the various indices/data
 using the following manner:
 
	print(Prot.count_amino_acids())
	print(Prot.get_amino_acids_percent())
	print(Prot.molecular_weight())
	print(Prot.get_amino_acids_MW_percent())
	print(Prot.aromaticity())
	print(Prot.instability_index())
	print(Prot.gravy_index())
	print(Prot.isoelectric_point())
	print(Prot.secondary_structure_fraction())
	
Have fun checking out some protein specs!

When used with the Phylogenetic_Analysis program, the various indices will be called up automatically
Add in this feature!!!!!!!!!!
"""

from __future__ import print_function

import sys
import Protein_Param_Data
import Isoelectric_Point
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData


class ProteinAnalysis(object):

	def __init__(self, prot_sequence):
		if prot_sequence.islower():
			self.sequence = Seq(prot_sequence.upper(), IUPAC.protein)
		else:
			self.sequence = Seq(prot_sequence, IUPAC.protein)
		self.amino_acids_content = None
		self.amino_acids_percent = None
		self.amino_acids_MW_percent = None
		self.length = len(self.sequence)

	def count_amino_acids(self):
		"""Count the number of times each amino acid appears in the protein.  Simple analysis stuff
		"""
		if self.amino_acids_content is None:
			prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
			for aa in prot_dic:
				prot_dic[aa] = self.sequence.count(aa)

			self.amino_acids_content = prot_dic

		return self.amino_acids_content

	def get_amino_acids_percent(self):
		"""Calculate the amino acid percent content of the protein.  Going one step further than
		just simple counts
		"""
		if self.amino_acids_percent is None:
			aa_counts = self.count_amino_acids()

			percentages = {}
			for aa in aa_counts:
				percentages[aa] = aa_counts[aa] / float(self.length) * 100

			self.amino_acids_percent = percentages

		return self.amino_acids_percent

	def molecular_weight(self):
		"""Calculate the total molecular weight of the protein.  Obviously this is an important
		bit of data, as we always need to know how big our protein is.  You can also get 
		the individual AA weights by looking them up in IUPAC
		"""
		iupac_weights = IUPACData.protein_weights
		water = 18.02

		aa_weights = {}
		for i in iupac_weights:
			aa_weights[i] = iupac_weights[i] - water

		total_weight = water
		for aa in self.sequence:
			total_weight += aa_weights[aa]

		return total_weight

	def get_amino_acids_MW_percent(self):
		"""Calculate the amino acid percent content based on MW instead of raw counts.  Because
		the weight of the individual amino acids can vary so much, this is useful for protein cleavage
		analyses to help predict daughter product composition and clipping
		"""
		if self.amino_acids_MW_percent is None:
			aa_mws = self.count_amino_acids()

			iupac_weights = IUPACData.protein_weights
			water = 18.02

			aa_weights = {}
			for i in iupac_weights:
				aa_weights[i] = iupac_weights[i] - water

			total_weight = water
			for aa in self.sequence:
				total_weight += aa_weights[aa]

			mw_percentages = {}
			for aa in aa_mws:
				mw_percentages[aa] = (aa_mws[aa] * aa_weights[aa]) / total_weight * 100

			self.amino_acids_MW_percent = mw_percentages

		return self.amino_acids_MW_percent

	def aromaticity(self):
		"""Calculate the aromaticity of the protein, which is the frequency of Phenylalanine (F), 
		Tyrosine (Y), and Tryptophan (W).  This value is sort of the opposite of the instability_index.
		Aromatic compounds are quite stable so having a higher aromaticity value would indicate 
		stability or, at the very least, will help you spot stable sections of the peptide/protein
		"""
		aromatic_aa = 'FYW'
		aa_percentages = self.get_amino_acids_percent()

		aromaticity = sum(aa_percentages[aa] for aa in aromatic_aa)

		return aromaticity

	def instability_index(self):
		"""Calculate the instability index of the protein.  A value over 40 indicates
		the protein is unstable and will have a short half-life.  This is a very useful value
		in biopharma and other protein based analyses and is especially handy in formulations
		when dealing with shelf life and temperature dependent stability
		"""
		index = Protein_Param_Data.DIWV
		score = 0.0

		for i in range(self.length - 1):
			this, next = self.sequence[i:i+2]
			dipeptide_value = index[this][next]
			score += dipeptide_value

		return (10.0 / self.length) * score

	def gravy_index(self):
		"""Calculate the gravy index of the protein...mmm...gravy...
		Gravy is the grand average of hydropathicity of the protein, with a positive number
		indicating an overall hydrophobic protein and negative indicating an overall hydrophilic
		protein.  Also very useful in biopharma and can be used to characterize the behavior of
		membrane or other proteins
		"""
		total_gravy = sum(Protein_Param_Data.kd[aa] for aa in self.sequence)

		return total_gravy / self.length

	def isoelectric_point(self):
		"""Calculate the isoelectric point of the protein.  This is a very useful piece of
		information for formulation and other biopharma scientists, like myself.  It tells you
		the pH value at which the protein has neither a positive nor negative charge
		"""
		aa_content = self.count_amino_acids()

		ie_point = Isoelectric_Point.IsoelectricPoint(self.sequence, aa_content)
		return ie_point.pi()

	def secondary_structure_fraction(self):
		"""Calculate the fraction of helix, bent and sheet secondary structures
		within the protein sequence.  Secondary structures are exceedingly important in determining
		binding sites and other behaviors of the protein

		Amino acids in helix: V, I, Y, F, W, L
		Amino acids in bent: N, P, G, S
		Amino acids in sheet: E, M, A, L

		**Note** The secondary structures are predictions only
		"""
		aa_percentages = self.get_amino_acids_percent()

		helix = sum(aa_percentages[r] for r in 'VIYFWL')
		bent  = sum(aa_percentages[r] for r in 'NPGS')
		sheet = sum(aa_percentages[r] for r in 'EMAL')

		return helix, bent, sheet


