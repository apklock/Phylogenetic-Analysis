#! /usr/bin/env python

from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml
import os

""" Need to connect the various modules to each other.
    Protein_Analyzer takes an AA sequence as its' input, can get that from GenBankProtein.
    """
    

AccessionNumbers = raw_input('Please enter all accession numbers here: ')
AN = AccessionNumbers.split(', ')

AccessionType = raw_input('Please indicate accession type (DNA/RNA/Protein): ')


if AccessionType is 'DNA':
    gb_dna = genbank.GenBankDna(ids=AN)
    gb_type = gb_dna
    model_type = ['-m', 'GTRCAT', '-N', '250']
    
elif AccessionType is 'RNA':
    gb_rna = genbank.GenBankRna(ids=AN)
    gb_type = gb_rna
    model_type = ['-m', 'GTRCAT', '-N', '250']

else:
    gb_prot = genbank.GenBankProtein(ids=AN)
    gb_type = gb_prot
    model_type = ['-m', 'PROTCATGTR', '-N', '250']


for gb in gb_type:
    print gb


os.system('pause')


data = gb_type.generate_char_matrix(
    label_components=["accession", "organism"])
data = muscle.muscle_align(data)
rr = raxml.RaxmlRunner()
tree = rr.estimate_tree(data, model_type)
print tree.as_ascii_plot()


os.system('pause')