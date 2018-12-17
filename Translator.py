'''
Useful for validation: http://www.attotron.com/cybertory/analysis/trans.htm
'''

class Translator:
    rna_codons = {
        "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
        "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
        "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
        "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
        "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
        "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
        "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "UAA": "STOP", "CAA": "Q", "AAA": "K", "GAA": "E",
        "UAG": "STOP", "CAG": "Q", "AAG": "K", "GAG": "E",
        "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
        "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "UGA": "STOP", "CGA": "R", "AGA": "R", "GGA": "G",
        "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }

    def __init__(self, rna_sequence: str) -> None:
        allowed_symbols = set(['A', 'C', 'G', 'T', 'U'])
        symbols = set(rna_sequence.upper())
        assert symbols <= allowed_symbols, 'Sequence contains invalid nucleotide symbols!'

        # If find(arg) returns -1 -> arg not present in string
        start_codon_pos: int = rna_sequence.find('AUG')
        stop_codon_pos: int = max(rna_sequence.find('UAA'), rna_sequence.find('UAG'), rna_sequence.find('UGA'))

        assert start_codon_pos > -1, 'Sequence does not contain a start codon!'
        assert stop_codon_pos > -1, 'Sequence does not contain any of stop codons!'
        # Notice that AUGAUG (duplicated start codons) can be read as aUGA (stop codon)
        # If that situation is impossible, then we should replace 1 with 3
        assert start_codon_pos + 1 <= stop_codon_pos, 'Start codon must be placed before stop codon!'

        self.rna_seq = rna_sequence

    @property
    def to_protein(self) -> str:
        '''Terminates when stop codon is found'''
        protein: str = ''

        # 1. Ignore everything before start codon
        start_codon_pos = self.rna_seq.find('AUG')

        # 2. Translate codons to amino acids until no stop codon was found
        for nucl_pos in range(start_codon_pos, len(self.rna_seq), 3):
            codon = self.rna_seq[nucl_pos:nucl_pos + 3]
            if self.rna_codons[codon] == 'STOP':
                break
            protein += self.rna_codons[codon]

        return protein
