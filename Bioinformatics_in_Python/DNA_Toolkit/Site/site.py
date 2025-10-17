import sys
import os
import streamlit as st

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now you can import your local packages
import bio_seq
import bio_structs
import utilities

st.header("Genome Toolit Plug-in Tool", divider = "violet")

seq = bio_seq.bio_seq()
seq.generate_rnd_seq(length = 100, seq_type = "DNA")

st.markdown("Randomly generated DNA string:")
st.markdown(f":red[Sequence:] {getattr(seq, "seq")}")
st.markdown(f":red[Sequence length:] {getattr(seq, "seq_len")}")
st.markdown(f":red[Sequence type:] {getattr(seq, "seq_type")}")
st.markdown(f":red[Sequence label:] {getattr(seq, "label")}")

st.subheader("Sequence modification and analysis:")
st.badge("Import your own FASTA.txt file below",icon ="🧬")

uploaded_file = st.file_uploader("Choose a file", type="txt")
type_of_seq = st.radio("What is the sequence type?", ["DNA","RNA","Protein"])

st.write(type(uploaded_file))
uploaded_file = utilities.parse_fasta_streamlit(uploaded_file, seq_type=type_of_seq)


