#! /usr/bin/env python

from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml
import os
import Protein_Analyzer

""" Just needs a bit of polishing
    """
    

AccessionNumbers = raw_input('Please enter all accession numbers here: ')
AN = AccessionNumbers.split(', ')

AccessionType = raw_input('Please indicate accession type (DNA/RNA/Protein): ')
if AccessionType.islower():
    AccessionType = AccessionType.upper()
else:
    AccessionType = AccessionType

if AccessionType == 'DNA':
    gb_dna = genbank.GenBankDna(ids=AN)
    gb_type = gb_dna
    model_type = ['-m', 'GTRCAT', '-N', '250']
    
elif AccessionType == 'RNA':
    gb_rna = genbank.GenBankRna(ids=AN)
    gb_type = gb_rna
    model_type = ['-m', 'GTRCAT', '-N', '250']

else:
    gb_prot = genbank.GenBankProtein(ids=AN)
    gb_type = gb_prot
    model_type = ['-m', 'PROTCATGTR', '-N', '250']

print '\n Here are your query results:\n'

for gb in gb_type:
    print gb
    print ''


os.system('pause')


# data = gb_type.generate_char_matrix(
    # label_components=["accession", "organism"])
# data = muscle.muscle_align(data)
# rr = raxml.RaxmlRunner()
# tree = rr.estimate_tree(data, model_type)
# print tree.as_ascii_plot()


os.system('pause')


if AccessionType == 'DNA':
    print 'You have just made a phylogenetic tree using DNA sequences.'
    print 'That is awesome and YOU are awesome!'
    print 'Thank you and have a fantastic day.'

elif AccessionType == 'RNA':
    print 'You have just made a phylogenetic tree using RNA sequences.'
    print 'That is awesome and YOU are awesome!'
    print 'Thank you and have a fantastic day.'

else:
    for gb in gb_type:
        Prot = Protein_Analyzer.ProteinAnalysis(gb.sequence_text)
        print '\n Here is an analysis of' 
        print gb
        print '\n Amino acid counts:'
        print(Prot.count_amino_acids())
        print '\n Amino acid percentages based on count:'
        print(Prot.get_amino_acids_percent())
        print '\n Total protein weight:'
        print(Prot.molecular_weight())
        print '\n Amino acid percentages based on molecular weight:'
        print(Prot.get_amino_acids_MW_percent())
        print '\n Aromaticity index:'
        print(Prot.aromaticity())
        print '\n Instability index:'
        print(Prot.instability_index())
        print '\n Gravy index:'
        print(Prot.gravy_index())
        print '\n Isoelectric point:'
        print(Prot.isoelectric_point())
        print '\n Fractions of predicted secondary structures:'
        print(Prot.secondary_structure_fraction())

