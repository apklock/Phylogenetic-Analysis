#! /usr/bin/env python

from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml

"""Need to add in entry format for the user and way to parse the input into the specific categories.
    Then need to figure out how to cut the accession numbers into id_range and prefix (dict?).
    Then need to connect the various modules to each other.
    Protein_Analyzer takes an AA sequence as its' input, can get that from GenBankProtein.
    """
    

if """entry is DNA"""
    gb = genbank.GenBankDna(id_range=(), prefix="")
elif """entry is RNA"""
    gb = genbank.GenBankRna(id_range=(), prefix="")
else """entry is Protein"""
    gb = genbank.GenBankProtein(id_range=(), prefix="")


data = gb.generate_char_matrix(
    label_components=["accession", "organism"])
data = muscle.muscle_align(data)
rr = raxml.RaxmlRunner()
tree = rr.estimate_tree(data, ['-N', '250'])
print tree.as_ascii_plot()

