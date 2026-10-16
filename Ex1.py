import sys
from collections import Counter

rna_codon_table = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def read_fasta(file_path):
    """Reads a FASTA file and returns the continuous genome sequence."""
    sequence = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):
                    sequence.append(line)
        return "".join(sequence).upper()
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)

def get_top_codons(sequence, top_n=10):
    """Uses a sliding window of 3 positions to count triplet frequencies."""
    triplets = [sequence[i:i+3] for i in range(len(sequence) - 2)]
    
    valid_triplets = [t for t in triplets if set(t).issubset({'A', 'C', 'G', 'T', 'U'})]
    
    freq = Counter(valid_triplets)
    return freq.most_common(top_n)

def analyze_genome(file_path):
    """Main application logic to read, count, translate, and display results."""
    print(f"Reading genome from: {file_path}...\n")
    dna_sequence = read_fasta(file_path)
    
    if not dna_sequence:
        print("No valid sequence data found.")
        return

    top_10_triplets = get_top_codons(dna_sequence, 10)

    print(f"{'DNA Triplet':<15} | {'RNA Codon':<15} | {'Amino Acid':<15} | {'Frequency'}")
    print("-" * 65)

    for dna_triplet, count in top_10_triplets:
        rna_codon = dna_triplet.replace('T', 'U')
        
        amino_acid = rna_codon_table.get(rna_codon, "Unknown")
        
        print(f"{dna_triplet:<15} | {rna_codon:<15} | {amino_acid:<15} | {count}")

if __name__ == "__main__":
    fasta_file = "genome.fasta" 
    analyze_genome(fasta_file)


# The aminoacids found in the test genome are:
#   Phenylalanine, Leucine, and Lysine (Essential Amino Acids): The human body cannot synthesize these, so they must be ingested directly from complete protein sources like beef, poultry, pork, fish, eggs, dairy, and soy products.

#   Arginine: Highly concentrated in turkey, pork, chicken, pumpkin seeds, peanuts, and chickpeas.

#   Glutamine: Found abundantly in beef, chicken, fish, eggs, milk, cabbage, and spinach.

#   Cysteine: Available in pork, beef, chicken, tuna, lentils, oatmeal, and sunflower seeds.

#   Alanine: Easily sourced from all major meats, seafood, eggs, dairy, and beans.

#   Glycine: Most prominent in collagen-rich connective tissues, bone broth, gelatin, poultry skin, and legumes.