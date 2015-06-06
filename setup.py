from setuptools import setup

setup(
    name='Phylogenetic Analysis',
    version='1.0',
    author='Andrew Klock',
    author_email='apklock@gmail.com',
    packages=['Phylogenetic-Analysis'],
    scripts=['Phylogenetic_Analysis.py', 'Protein_Analyzer.py', 'Isoelectric_Point.py', 'Protein_Param_Data.py'],
    url='https://github.com/apklock/Phylogenetic-Analysis',
    license='LICENSE.txt',
    description='An awesome package that creates phylogenetic trees and analyzes proteins',
    long_description=open('README.txt').read(),
    install_requires=[
        "biopython >= 1.61",
        "dendropy",
        "numpy",
        "muscle",
        "raxml",
        "scipy",
    ],
)