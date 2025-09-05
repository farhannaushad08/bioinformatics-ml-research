# Pairwise sequence alignment using Biopython with BLOSUM62
# Calculates number of matches and percentage similarity

from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

# Load substitution matrix
Blosum62 = substitution_matrices.load("BLOSUM62")

# Read FASTA sequences
seq1 = SeqIO.read("data/sample_sequences/seq1.fasta", "fasta")  # Gorilla
seq2 = SeqIO.read("data/sample_sequences/seq2.fasta", "fasta")  # Physeter catodon
seq3 = SeqIO.read("data/sample_sequences/seq3.fasta", "fasta")  # Grus americana
seq4 = SeqIO.read("data/sample_sequences/seq4.fasta", "fasta")  # Homo sapiens

# ---- Pairwise Global Alignments ----
alignments = [
    ("seq1 vs seq2", pairwise2.align.globalds(seq1.seq, seq2.seq, Blosum62, -10, -0.5)),
    ("seq1 vs seq3", pairwise2.align.globalds(seq1.seq, seq3.seq, Blosum62, -10, -0.5)),
    ("seq1 vs seq4", pairwise2.align.globalds(seq1.seq, seq4.seq, Blosum62, -10, -0.5)),
    ("seq2 vs seq3", pairwise2.align.globalds(seq2.seq, seq3.seq, Blosum62, -10, -0.5)),
    ("seq2 vs seq4", pairwise2.align.globalds(seq2.seq, seq4.seq, Blosum62, -10, -0.5)),
    ("seq3 vs seq4", pairwise2.align.globalds(seq3.seq, seq4.seq, Blosum62, -10, -0.5)),
]

# ---- Calculate and Print Results ----
def calculate_similarity(aln):
    seq_a, seq_b = aln[0], aln[1]
    matches = sum(1 for a, b in zip(seq_a, seq_b) if a == b)
    percentage = round(matches / len(seq_a) * 100, 2)
    return matches, percentage

print("\n=== Global Alignment with BLOSUM62 ===")
for label, aln in alignments:
    best = aln[0]
    matches, percentage = calculate_similarity(best)
    print(f"\n{label}")
    print(pairwise2.format_alignment(*best))
    print(f"Matches: {matches}")
    print(f"Similarity: {percentage}%")
