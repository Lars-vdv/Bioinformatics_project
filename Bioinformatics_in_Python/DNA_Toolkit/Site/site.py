import sys
import streamlit as st
import importlib.util

st.write(sys.executable)
st.write(sys.path)

# import bio_seq 
#import utilities
#import bio_structs
#import bio_seq

# spec = importlib.util.spec_from_file_location("bio_seq", r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\bio_seq.py")
# foo = importlib.util.module_from_spec(spec)
# sys.modules["bio_seq"] = foo
# spec.loader.exec_module(foo)
# foo.bio_seq()

# #import bio_structs
# spec2 = importlib.util.spec_from_file_location("bio_structs", r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\bio_structs.py")
# foo2 = importlib.util.module_from_spec(spec2)
# sys.modules["bio_structs"] = foo2
# spec2.loader.exec_module(foo2)
# foo2.bio_seq()



st.header("Genome Toolit Plug-in Tool", divider = "violet")
