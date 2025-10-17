import sys
import os
import streamlit as st

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now you can import your local packages
import bio_seq
import bio_structs
import utilities

st.write(sys.executable)
st.write(sys.path)

st.header("Genome Toolit Plug-in Tool", divider = "violet")

seq = bio_seq.bio_seq()
seq.generate_rnd_seq(length = 100, seq_type = "DNA")
st.write(seq)