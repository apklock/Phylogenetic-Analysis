"""Calculate isoelectric points of polypeptides using methods from Bjellqvist et al. 1993 and 1994
""" 

positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98} 
negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0} 
pKcterminal = {'D': 4.55, 'E': 4.75} 
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7} 
charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y') 


class IsoelectricPoint(object):
	def __init__(self, ProteinSequence, AminoAcidsContent):
		self.sequence = ProteinSequence
		self.charged_aas_content = self._select_charged(AminoAcidsContent)

	def _select_charged(self, AminoAcidsContent):
		charged = {}
		for aa in charged_aas:
			charged[aa] = float(AminoAcidsContent[aa])
		charged['Nterm'] = 1.0
		charged['Cterm'] = 1.0
		return charged

	def _chargeR(self, pH, pos_pKs, neg_pKs):
		PositiveCharge = 0.0
		for aa, pK in pos_pKs.items():
			CR = 10**(pK-pH)
			partial_charge = CR/(CR+1.0)
			PositiveCharge += self.charged_aas_content[aa] * partial_charge

		NegativeCharge = 0.0
		for aa, pK in neg_pKs.items():
			CR = 10**(pH-pK)
			partial_charge = CR/(CR+1.0)
			NegativeCharge += self.charged_aas_content[aa] * partial_charge

		return PositiveCharge - NegativeCharge

	def pi(self):
		pos_pKs = dict(positive_pKs)
		neg_pKs = dict(negative_pKs)
		nterm = self.sequence[0]
		cterm = self.sequence[-1]
		if nterm in pKnterminal:
			pos_pKs['Nterm'] = pKnterminal[nterm]
		if cterm in pKcterminal:
			neg_pKs['Cterm'] = pKcterminal[cterm]

		pH = 7.0
		Charge = self._chargeR(pH, pos_pKs, neg_pKs)
		if Charge > 0.0:
			pH1 = pH
			Charge1 = Charge
			while Charge1 > 0.0:
				pH = pH1 + 1.0
				Charge = self._chargeR(pH, pos_pKs, neg_pKs)
				if Charge > 0.0:
					pH1 = pH
					Charge1 = Charge
				else:
					pH2 = pH
					Charge2 = Charge
					break
		else:
			pH2 = pH
			Charge2 = Charge
			while Charge2 < 0.0:
				pH = pH2 - 1.0
				Charge = self._chargeR(pH, pos_pKs, neg_pKs)
				if Charge < 0.0:
					pH2 = pH
					Charge2 = Charge
				else:
					pH1 = pH
					Charge1 = Charge
					break

		while pH2 - pH1 > 0.0001 and Charge != 0.0:
			pH = (pH1 + pH2) / 2.0
			Charge = self._chargeR(pH, pos_pKs, neg_pKs)
			if Charge > 0.0:
				pH1 = pH
				Charge1 = Charge
			else:
				pH2 = pH
				Charge2 = Charge

		return pH