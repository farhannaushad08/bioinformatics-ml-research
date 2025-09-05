# Comprehensive pairwise sequence alignment with multiple matrices

from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

# Load substitution matrices
blosum62 = substitution_matrices.load("BLOSUM62")
pam250 = substitution_matrices.load("PAM250")

# Read FASTA sequences
seq1 = SeqIO.read("data/sample_sequences/seq1.fasta", "fasta")
seq2 = SeqIO.read("data/sample_sequences/seq2.fasta", "fasta")
seq3 = SeqIO.read("data/sample_sequences/seq3.fasta", "fasta")
seq4 = SeqIO.read("data/sample_sequences/seq4.fasta", "fasta")

# Helper to print alignment
def print_alignment(title, aln):
    print(f"\n{title}")
    print(pairwise2.format_alignment(*aln[0]))

# ---- Global Alignment with BLOSUM62 ----
print("\n=== Global Alignment with BLOSUM62 ===")
print_alignment("seq1 vs seq2", pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5))
print_alignment("seq1 vs seq3", pairwise2.align.globalds(seq1.seq, seq3.seq, blosum62, -10, -0.5))
print_alignment("seq1 vs seq4", pairwise2.align.globalds(seq1.seq, seq4.seq, blosum62, -10, -0.5))
print_alignment("seq2 vs seq3", pairwise2.align.globalds(seq2.seq, seq3.seq, blosum62, -10, -0.5))
print_alignment("seq2 vs seq4", pairwise2.align.globalds(seq2.seq, seq4.seq, blosum62, -10, -0.5))
print_alignment("seq3 vs seq4", pairwise2.align.globalds(seq3.seq, seq4.seq, blosum62, -10, -0.5))

# ---- Global Alignment with PAM250 ----
print("\n=== Global Alignment with PAM250 ===")
print_alignment("seq1 vs seq2", pairwise2.align.globalds(seq1.seq, seq2.seq, pam250, -10, -0.5))
print_alignment("seq1 vs seq3", pairwise2.align.globalds(seq1.seq, seq3.seq, pam250, -10, -0.5))
print_alignment("seq1 vs seq4", pairwise2.align.globalds(seq1.seq, seq4.seq, pam250, -10, -0.5))

# ---- Local Alignment with BLOSUM62 ----
print("\n=== Local Alignment with BLOSUM62 ===")
print_alignment("seq2 vs seq3", pairwise2.align.localds(seq2.seq, seq3.seq, blosum62, -10, -0.5))
print_alignment("seq3 vs seq4", pairwise2.align.localds(seq3.seq, seq4.seq, blosum62, -10, -0.5))

# ---- Local Alignment with PAM250 ----
print("\n=== Local Alignment with PAM250 ===")
print_alignment("seq2 vs seq4", pairwise2.align.localds(seq2.seq, seq4.seq, pam250, -10, -0.5))
