#! /usr/bin/env python
 
from dendropy.interop import genbank
from dendropy.interop import muscle
from dendropy.interop import raxml
gb = genbank.GenBankProtein(id_range=(56652, 56759), prefix="AAX")
data = gb.generate_char_matrix(
        label_components=["accession", "organism"])
data = muscle.muscle_align(data)
rr = raxml.RaxmlRunner()
tree = rr.estimate_tree(data, ['-N', '250'])
print tree.as_ascii_plot()