import sys
import os
import streamlit as st
import io
import pandas as pd

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now you can import your local packages
import bio_seq
import bio_structs
import utilities

st.header("Genome Toolit Plug-in Tool", divider = "violet")

rnd_seq = bio_seq.bio_seq()
rnd_seq.generate_rnd_seq(length = 100, seq_type = "DNA")

st.markdown("Randomly generated DNA string:")
st.markdown(f":red[Sequence:] {getattr(rnd_seq, "seq")}")
st.markdown(f":red[Sequence length:] {getattr(rnd_seq, "seq_len")}")
st.markdown(f":red[Sequence type:] {getattr(rnd_seq, "seq_type")}")
st.markdown(f":red[Sequence label:] {getattr(rnd_seq, "label")}")

st.subheader("Sequence modification and analysis:")
st.markdown("*For now only single sequence functionality*")

#upload custom seq:
verified = False
message = False

while verified == False:
    seq_upload = st.text_input("Enter sequence below:")
    seq_upload = seq_upload.strip()
    seq_upload = seq_upload.replace(" ", "")

    seq_upload_label = st.text_input("Enter sequence label below:")
    seq_upload_type = st.radio("Enter sequence type below:",["DNA","RNA","Protein"])

    #verify custom seq validity:
   
    custom_sequence = bio_seq.bio_seq(seq = seq_upload, seq_type= seq_upload_type, label = seq_upload_label)
    st.write(f'label: {getattr(custom_sequence, "label")}')
    st.write(f'type: {getattr(custom_sequence, "seq_type")}')
    st.write(f'seq: {getattr(custom_sequence, "seq")[0:100]} ...')
    
    
    #verify if user wants to continue with the seq:
    if len(seq_upload) != 0:

        check = st.radio("Is this sequence correct?",["Yes","No"])

    if check == "Yes":
        verified = True
        message = True

    else:
        st.write("Try again")
        break

if message == True:
    st.subheader("Check completed, continue below!")

custom_sequence_nuc_freq = custom_sequence.countNucFrequency()
custom_sequence_nuc_freq_df = pd.Series(custom_sequence_nuc_freq, name="Count").to_frame()
# Optionally reset index to have a column named 'nucleotide'
custom_sequence_nuc_freq_df = custom_sequence_nuc_freq_df.rename_axis("Nucleotide").reset_index()
# Result: columns ['nucleotide', 'count']
st.write("Nucleotide frequency:")
st.dataframe(custom_sequence_nuc_freq_df)

custom_sequence_gc_content = custom_sequence.gc_content_round()
st.write("GC Content:")
st.badge(custom_sequence_gc_content)

custom_sequence_gc_content_subseq = custom_sequence.gc_content_subsec(k=100)
custom_sequence_gc_content_subseq = pd.Series(custom_sequence_gc_content_subseq, name="Percentage").to_frame()
custom_sequence_gc_content_subseq = custom_sequence_gc_content_subseq.rename_axis("Subseq").reset_index()
st.write("GC Content per subseq (k=100):")
st.dataframe(custom_sequence_gc_content_subseq)

custom_sequence_translate = custom_sequence.translate_seq()
custom_sequence_translate = "".join(custom_sequence_translate)
st.write("Sequence translation:")
st.write(custom_sequence_translate)

custom_sequence_proteins = bio_seq.bio_seq(seq=custom_sequence)
custom_sequence_proteins = custom_sequence_proteins.all_proteins_from_ORF(ordered=True)


st.write("All proteins from ORFs:")
st.dataframe(custom_sequence_proteins)


    # st.badge("Import your own FASTA.txt file below",icon ="🧬")

    # uploaded_file = st.file_uploader("Choose a file", type="txt",label_visibility="visible", key="my_upload")
    # type_of_seq = st.radio("What is the sequence type?", ["DNA","RNA","Protein"])

    # if uploaded_file is not None:
    #     # wrap the UploadedFile (a bytes-file) as a text file
    #     uploaded_file.seek(0)
    #     adjusted_file = uploaded_file.getvalue().decode("utf-8")
    #     adjusted_file = io.StringIO(adjusted_file)


    # sequences= utilities.parse_fasta_streamlit(adjusted_file, seq_type=type_of_seq, create_variables=True)
    # #sequences_info = bio_seq.bio_seq.get_seq_info(sequences)

    # st.write(bio_seq.bio_seq.get_seq_info(seq_1))
    # st.write(print(seq_1))

    # for i in range(len(sequences_info)):
    #     st.write(sequences_info[i])
