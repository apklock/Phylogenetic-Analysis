#! /usr/bin/env python

from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml
import os
import Protein_Analyzer

""" Need to connect the various modules to each other.
    Protein_Analyzer takes an AA sequence as its' input, can get that from GenBankProtein.
    gb.sequence_text gets sequence as a string!!! FINALLY!!!!
    """
    

AccessionNumbers = raw_input('Please enter all accession numbers here: ')
AN = AccessionNumbers.split(', ')

AccessionType = raw_input('Please indicate accession type (DNA/RNA/Protein): ')


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

print ''
print 'Here are your query results:\n'

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
        print(Prot.count_amino_acids())
        print(Prot.get_amino_acids_percent())
        print(Prot.molecular_weight())
        print(Prot.get_amino_acids_MW_percent())
        print(Prot.aromaticity())
        print(Prot.instability_index())
        print(Prot.gravy_index())
        print(Prot.isoelectric_point())
        print(Prot.secondary_structure_fraction())

