"""This module contains hydrophobicity indices according to Kyte and Doolittle and 
DIWV indices for instability calculation according to Guruprasad et al. 1990
"""
   
kd = {'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5, 'Q':-3.5, 'E':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5, 'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2 } 

DIWV = {'A': {'A': 1.0, 'C': 44.94, 'E': 1.0, 'D': -7.49, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': -7.49, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0}, 
	'C': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 20.26, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 33.60, 'K': 1.0, 'M': 33.60, 'L': 20.26, 'N': 1.0, 'Q': -6.54, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 33.60, 'W': 24.68, 'V': -6.54, 'Y': 1.0}, 
	'E': {'A': 1.0, 'C': 44.94, 'E': 33.60, 'D': 20.26, 'G': 1.0, 'F': 1.0, 'I': 20.26, 'H': -6.54, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': 1.0, 'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0}, 
	'D': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0, 'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 1.0, 'S': 20.26, 'R': -6.54, 'T': -14.03, 'W': 1.0, 'V': 1.0, 'Y': 1.0}, 
	'G': {'A': -7.49, 'C': 1.0, 'E': -6.54, 'D': 1.0, 'G': 13.34, 'F': 1.0, 'I': -7.49, 'H': 1.0, 'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': -7.49, 'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': 13.34, 'V': 1.0, 'Y': -7.49}, 
	'F': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 13.34, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -14.03, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 33.601}, 
	'I': {'A': 1.0, 'C': 1.0, 'E': 44.94, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 13.34, 'K': -7.49, 'M': 1.0, 'L': 20.26, 'N': 1.0, 'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0}, 
	'H': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -9.37, 'F': -9.37, 'I': 44.94, 'H': 1.0, 'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 24.68, 'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': -6.54, 'W': -1.88, 'V': 1.0, 'Y': 44.94}, 
	'K': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -7.49, 'F': 1.0, 'I': -7.49, 'H': 1.0, 'K': 1.0, 'M': 33.60, 'L': -7.49, 'N': 1.0, 'Q': 24.64, 'P': -6.54, 'S': 1.0, 'R': 33.60, 'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0}, 
	'M': {'A': 13.34, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 58.28, 'K': 1.0, 'M': -1.88, 'L': 1.0, 'N': 1.0, 'Q': -6.54, 'P': 44.94, 'S': 44.94, 'R': -6.54, 'T': -1.88, 'W': 1.0, 'V': 1.0, 'Y': 24.68}, 
	'L': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 33.60, 'P': 20.26, 'S': 1.0, 'R': 20.26, 'T': 1.0, 'W': 24.68, 'V': 1.0, 'Y': 1.0}, 
	'N': {'A': 1.0, 'C': -1.88, 'E': 1.0, 'D': 1.0, 'G': -14.03, 'F': -14.03, 'I': 44.94, 'H': 1.0, 'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': -6.54, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 1.0}, 
	'Q': {'A': 1.0, 'C': -6.54, 'E': 20.26, 'D': 20.26, 'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': -6.54, 'Y': -6.54}, 
	'P': {'A': 20.26, 'C': -6.54, 'E': 18.38, 'D': -6.54, 'G': 1.0, 'F': 20.26, 'I': 1.0, 'H': 1.0, 'K': 1.0, 'M': -6.54, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': -6.54, 'T': 1.0, 'W': -1.88, 'V': 20.26, 'Y': 1.0}, 
	'S': {'A': 1.0, 'C': 33.60, 'E': 20.26, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 44.94, 'S': 20.26, 'R': 20.26, 'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0}, 
	'R': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 20.26, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 13.34, 'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 58.28, 'T': 1.0, 'W': 58.28, 'V': 1.0, 'Y': -6.54}, 
	'T': {'A': 1.0, 'C': 1.0, 'E': 20.26, 'D': 1.0, 'G': -7.49, 'F': 13.34, 'I': 1.0, 'H': 1.0, 'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': -14.03, 'Q': -6.54, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0}, 
	'W': {'A': -14.03, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -9.37, 'F': 1.0, 'I': 1.0, 'H': 24.68, 'K': 1.0, 'M': 24.68, 'L': 13.34, 'N': 13.34, 'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': -14.03, 'W': 1.0, 'V': -7.49, 'Y': 1.0}, 
	'V': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': -14.03, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -1.88, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': 1.0, 'V': 1.0, 'Y': -6.54}, 
	'Y': {'A': 24.68, 'C': 1.0, 'E': -6.54, 'D': 24.68, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 13.34, 'K': 1.0, 'M': 44.94, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 13.34, 'S': 1.0, 'R': -15.91, 'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 13.34}, 
	} 