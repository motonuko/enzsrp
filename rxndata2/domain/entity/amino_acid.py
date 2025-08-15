amino_acids = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp',
    'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly',
    'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
    'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',

    'U': 'Sec', 'O': 'Pyl'  # Not standard amino acids
}


def one_to_three(code):
    if code in amino_acids:
        return amino_acids[code]
    else:
        raise ValueError("Invalid single-letter amino acid code")


def three_to_one(code):
    reverse_dict = {v: k for k, v in amino_acids.items()}
    if code in reverse_dict:
        return reverse_dict[code]
    else:
        raise ValueError("Invalid three-letter amino acid code")
