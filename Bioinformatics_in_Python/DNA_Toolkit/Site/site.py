import sys
import os
import streamlit as st
import io

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

uploaded_file = st.file_uploader("Choose a file", type="txt",label_visibility="visible", key="my_upload")
type_of_seq = st.radio("What is the sequence type?", ["DNA","RNA","Protein"])

if uploaded_file is not None:
    # wrap the UploadedFile (a bytes-file) as a text file
    uploaded_file.seek(0)
    adjusted_file = uploaded_file.getvalue().decode("utf-8")
    adjusted_file = io.StringIO(adjusted_file)


st.write(adjusted_file)

sequences= utilities.parse_fasta_streamlit(adjusted_file, seq_type=type_of_seq)

st.write(sequences)
