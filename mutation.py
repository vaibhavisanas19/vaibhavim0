import streamlit as st
from Bio import pairwise2
from Bio import SeqIO
from io import StringIO
import pandas as pd

# Streamlit App Title
st.title("🧬 Pairwise Sequence Alignment & Mutation Detection")

# Input FASTA sequences
st.subheader("Enter FASTA Sequences")
fasta_input = st.text_area(
    "Paste sequences in FASTA format:",
    """>Seq1
ATGCTAGCTAGGCT
>Seq2
ATGCGAGCTTGGCT""",
    height=150,
)


# Function to parse FASTA sequences
def parse_fasta(fasta_text):
    fasta_io = StringIO(fasta_text.strip())
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records if len(records) >= 2 else None


# Function for pairwise alignment
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments[0] if alignments else None


# Function to detect mutations (Fixed lowercase parameters)
def detect_mutations(seq_a, seq_b):
    mutation_list = []
    for i, (base1, base2) in enumerate(zip(seq_a, seq_b), start=1):
        if base1 != base2:
            if base1 == "-":
                mutation_list.append(f"Insertion at position {i}: {base2}")
            elif base2 == "-":
                mutation_list.append(f"Deletion at position {i}: {base1}")
            else:
                mutation_list.append(f"Substitution at position {i}: {base1} → {base2}")
    return mutation_list


# Run alignment and mutation detection when the button is clicked
if st.button("Analyze Sequences"):
    parsed_records = parse_fasta(fasta_input)

    if not parsed_records:
        st.error("❌ Please provide at least two valid FASTA sequences!")
    else:
        # Extract the first two sequences
        sequence1 = str(parsed_records[0].seq)
        sequence2 = str(parsed_records[1].seq)

        # Perform alignment
        best_alignment = align_sequences(sequence1, sequence2)

        if best_alignment:
            st.success("✅ Alignment successful!")

            # Display alignment
            st.subheader("🔬 Pairwise Alignment")
            aligned_seq_a = best_alignment.seqA  # Renamed variables to lowercase
            aligned_seq_b = best_alignment.seqB
            st.text(f"{aligned_seq_a}\n{''.join('|' if a == b else ' ' for a, b in zip(aligned_seq_a, aligned_seq_b))}\n{aligned_seq_b}")

            # Detect mutations
            st.subheader("🧬 Detected Mutations")
            mutation_results = detect_mutations(aligned_seq_a, aligned_seq_b)  # Fixed function call

            if mutation_results:
                df_mutations = pd.DataFrame({"Mutation": mutation_results})
                st.dataframe(df_mutations)
            else:
                st.info("✅ No mutations detected. Sequences are identical!")
        else:
            st.error("❌ No valid alignments found. Please check the input sequences.")
