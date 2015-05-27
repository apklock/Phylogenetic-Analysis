#! /usr/bin/env python

from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml

"""Add in more info for DNA/RNA/Protein? Like in the print gb code (gb.rec, etc.)
    Add in continue points?
    Then need to connect the various modules to each other.
    Protein_Analyzer takes an AA sequence as its' input, can get that from GenBankProtein.
    """
    

AccessionNumbers = raw_input('Please enter all accession numbers here: ')
AN = AccessionNumbers.split(', ')

AccessionType = raw_input('Please indicate accession species (DNA/RNA/Protein): ')


if AccessionType is 'DNA':
    gb_dna = genbank.GenBankDna(ids=AN)
    gb_type = gb_dna
    
elif AccessionType is 'RNA':
    gb_rna = genbank.GenBankRna(ids=AN)
    gb_type = gb_rna

else:
    gb_prot = genbank.GenBankProtein(ids=AN)
    gb_type = gb_prot


for gb in gb_type:
    print gb


data = gb.generate_char_matrix(
    label_components=["accession", "organism"])
data = muscle.muscle_align(data)
rr = raxml.RaxmlRunner()
tree = rr.estimate_tree(data, ['-N', '250'])
print tree.as_ascii_plot()

