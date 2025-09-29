#!/usr/bin/env python3
# grantham_app.py

import sys
import traceback
import importlib.metadata as metadata
import streamlit as st

# Page config
st.set_page_config(page_title="Grantham Score Calculator", layout="centered")
st.title("🔬 Grantham Score Calculator")
st.write("Compute the Grantham distance between two amino acids.")

# Debug info (hidden by default)
with st.expander("Debug / Environment"):
    st.write("Python:", sys.version)
    st.write("Streamlit:", st.__version__)
    if st.button("Show installed streamlit package"):
        try:
            st.write("streamlit", metadata.version("streamlit"))
        except metadata.PackageNotFoundError:
            st.error("Streamlit not found in this environment!")

# Grantham matrix (symmetric)
grantham_matrix = {
    ('A','C'):195, ('A','D'):126, ('A','E'):107, ('A','F'):113, ('A','G'):60,  ('A','H'):86,
    ('A','I'):94,  ('A','K'):106, ('A','L'):96,  ('A','M'):84,  ('A','N'):111, ('A','P'):27,
    ('A','Q'):91,  ('A','R'):112, ('A','S'):99,  ('A','T'):58,  ('A','V'):64,  ('A','W'):148,
    ('A','Y'):112,
    ('C','D'):154, ('C','E'):170, ('C','F'):205, ('C','G'):159, ('C','H'):174, ('C','I'):198,
    ('C','K'):202, ('C','L'):198, ('C','M'):196, ('C','N'):139, ('C','P'):169, ('C','Q'):154,
    ('C','R'):180, ('C','S'):112, ('C','T'):149, ('C','V'):192, ('C','W'):215, ('C','Y'):194,
    ('D','E'):45,  ('D','F'):177, ('D','G'):94,  ('D','H'):81,  ('D','I'):168, ('D','K'):101,
    ('D','L'):172, ('D','M'):160, ('D','N'):23,  ('D','P'):108, ('D','Q'):61,  ('D','R'):96,
    ('D','S'):65,  ('D','T'):85,  ('D','V'):152, ('D','W'):181, ('D','Y'):160,
    ('E','F'):140, ('E','G'):98,  ('E','H'):40,  ('E','I'):134, ('E','K'):56,  ('E','L'):138,
    ('E','M'):126, ('E','N'):42,  ('E','P'):93,  ('E','Q'):29,  ('E','R'):54,  ('E','S'):80,
    ('E','T'):65,  ('E','V'):121, ('E','W'):152, ('E','Y'):122,
    ('F','G'):153, ('F','H'):100, ('F','I'):21,  ('F','K'):102, ('F','L'):22,  ('F','M'):28,
    ('F','N'):158, ('F','P'):114, ('F','Q'):116, ('F','R'):97,  ('F','S'):155, ('F','T'):103,
    ('F','V'):50,  ('F','W'):40,  ('F','Y'):22,
    ('G','H'):98,  ('G','I'):135, ('G','K'):127, ('G','L'):138, ('G','M'):127, ('G','N'):80,
    ('G','P'):42,  ('G','Q'):87,  ('G','R'):125, ('G','S'):56,  ('G','T'):59,  ('G','V'):109,
    ('G','W'):184, ('G','Y'):147,
    ('H','I'):94,  ('H','K'):32,  ('H','L'):99,  ('H','M'):87,  ('H','N'):68,  ('H','P'):77,
    ('H','Q'):24,  ('H','R'):29,  ('H','S'):89,  ('H','T'):47,  ('H','V'):84,  ('H','W'):115,
    ('H','Y'):83,
    ('I','K'):102, ('I','L'):5,   ('I','M'):10,  ('I','N'):149, ('I','P'):95,  ('I','Q'):109,
    ('I','R'):97,  ('I','S'):142, ('I','T'):89,  ('I','V'):29,  ('I','W'):61,  ('I','Y'):33,
    ('K','L'):107, ('K','M'):95,  ('K','N'):94,  ('K','P'):103, ('K','Q'):53,  ('K','R'):26,
    ('K','S'):121, ('K','T'):78,  ('K','V'):97,  ('K','W'):110, ('K','Y'):85,
    ('L','M'):15,  ('L','N'):153, ('L','P'):98,  ('L','Q'):113, ('L','R'):102, ('L','S'):145,
    ('L','T'):92,  ('L','V'):32,  ('L','W'):61,  ('L','Y'):36,
    ('M','N'):142, ('M','P'):87,  ('M','Q'):101, ('M','R'):91,  ('M','S'):135, ('M','T'):81,
    ('M','V'):21,  ('M','W'):67,  ('M','Y'):36,
    ('N','P'):91,  ('N','Q'):46,  ('N','R'):86,  ('N','S'):46,  ('N','T'):65,  ('N','V'):133,
    ('N','W'):174, ('N','Y'):143,
    ('P','Q'):76,  ('P','R'):103, ('P','S'):74,  ('P','T'):38,  ('P','V'):68,  ('P','W'):147,
    ('P','Y'):110,
    ('Q','R'):43,  ('Q','S'):68,  ('Q','T'):42,  ('Q','V'):96,  ('Q','W'):130, ('Q','Y'):99,
    ('R','S'):110, ('R','T'):71,  ('R','V'):96,  ('R','W'):101, ('R','Y'):77,
    ('S','T'):58,  ('S','V'):64,  ('S','W'):177, ('S','Y'):144,
    ('T','V'):69,  ('T','W'):128, ('T','Y'):92,
    ('V','W'):88,  ('V','Y'):55,
    ('W','Y'):37
}

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Map one-letter codes to three-letter codes
aa_lookup = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"
}

# Create symmetric Grantham matrix
symmetric_grantham = {}
for (a1, a2), val in grantham_matrix.items():
    symmetric_grantham[(a1, a2)] = val
    symmetric_grantham[(a2, a1)] = val
for aa in amino_acids:
    symmetric_grantham[(aa, aa)] = 0

def grantham_score(wt, mu):
    return symmetric_grantham.get((wt, mu))

def get_one_letter_from_3letter(three_letter):
    reverse_lookup = {v: k for k, v in aa_lookup.items()}
    return reverse_lookup[three_letter]

# Dropdown options (three-letter codes)
aa_options_3letter = list(aa_lookup.values())

# UI columns
col1, col2 = st.columns(2)
with col1:
    wildtype_sel = st.selectbox("Select Wildtype Amino Acid", aa_options_3letter)
with col2:
    mutant_sel = st.selectbox("Select Mutant Amino Acid", aa_options_3letter)

wildtype = get_one_letter_from_3letter(wildtype_sel)
mutant = get_one_letter_from_3letter(mutant_sel)

# Calculate button
if st.button("Calculate"):
    try:
        score = grantham_score(wildtype, mutant)
        if score is None:
            st.warning(f"No Grantham score found for {wildtype} → {mutant}.")
        else:
            st.success(f"Grantham score between **{wildtype_sel} → {mutant_sel}** is: **{score}**")
    except Exception as e:
        st.error("An unexpected error occurred; see details below.")
        st.exception(e)
        st.text(traceback.format_exc())
