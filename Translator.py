class Translator:
    rna_codon = {
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

    @staticmethod
    def rna_sequence_to_protein(rna_sequence: str) -> str:
        protein = ''
        for i in range(0, len(rna_sequence) - (3 + len(rna_sequence) % 3), 3):
            if Translator.rna_codon[rna_sequence[i:i + 3]] == "STOP":
                break
            protein += Translator.rna_codon[rna_sequence[i:i + 3]]
        return protein
