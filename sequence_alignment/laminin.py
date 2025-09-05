# Pairwise sequence alignment using Biopython with PAM250 and BLOSUM62

from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

# Load substitution matrices
Pam250 = substitution_matrices.load("PAM250")
Blosum62 = substitution_matrices.load("BLOSUM62")

# Read FASTA sequences
seq1 = SeqIO.read("data/sample_sequences/seq1.fasta", "fasta")
seq2 = SeqIO.read("data/sample_sequences/seq2.fasta", "fasta")
seq3 = SeqIO.read("data/sample_sequences/seq3.fasta", "fasta")
seq4 = SeqIO.read("data/sample_sequences/seq4.fasta", "fasta")

# ---- Global Alignment with BLOSUM62 ----
alignments = [
    ("seq1 vs seq2", pairwise2.align.globalds(seq1.seq, seq2.seq, Blosum62, -11, -1)),
    ("seq1 vs seq3", pairwise2.align.globalds(seq1.seq, seq3.seq, Blosum62, -11, -1)),
    ("seq1 vs seq4", pairwise2.align.globalds(seq1.seq, seq4.seq, Blosum62, -11, -1)),
    ("seq2 vs seq3", pairwise2.align.globalds(seq2.seq, seq3.seq, Blosum62, -11, -1)),
    ("seq2 vs seq4", pairwise2.align.globalds(seq2.seq, seq4.seq, Blosum62, -11, -1)),
    ("seq3 vs seq4", pairwise2.align.globalds(seq3.seq, seq4.seq, Blosum62, -11, -1)),
]

print("\n=== Global Alignment with BLOSUM62 ===")
for label, aln in alignments:
    print(f"\n{label}")
    print(pairwise2.format_alignment(*aln[0]))

# ---- Global Alignment with PAM250 ----
alignments_pam = [
    ("seq1 vs seq2", pairwise2.align.globalds(seq1.seq, seq2.seq, Pam250, -10, -0.5)),
    ("seq1 vs seq3", pairwise2.align.globalds(seq1.seq, seq3.seq, Pam250, -10, -0.5)),
    ("seq1 vs seq4", pairwise2.align.globalds(seq1.seq, seq4.seq, Pam250, -10, -0.5)),
]

print("\n=== Global Alignment with PAM250 ===")
for label, aln in alignments_pam:
    print(f"\n{label}")
    print(pairwise2.format_alignment(*aln[0]))
