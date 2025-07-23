import pandas as pd
from Bio import SeqIO
import random

# === CONFIG ===
CSV_PATH = "Genomovirus_Proteins.csv"         # Metadata file
FASTA_PATH = "Genomovirus_Proteins.faa"       # Full FASTA file
OUTPUT_FASTA = "rep_subset_by_genus.faa"      # Output file

# === STEP 1: Load and filter metadata ===
df = pd.read_csv(CSV_PATH)

# Keep only rows with "rep" in Protein name (case-insensitive)
df_filtered = df[df['Protein'].str.contains('rep', case=False, na=False)]

# Sample up to 10 accessions per Genus
sampled_accessions = (
    df_filtered.groupby('Genus')['Accession']
    .apply(lambda x: x.sample(n=min(10, len(x)), random_state=42))
    .explode()
    .tolist()
)
sampled_set = set(sampled_accessions)

# === STEP 2: Filter FASTA ===
selected_records = []
for record in SeqIO.parse(FASTA_PATH, "fasta"):
    accession = record.id.split('|')[0].strip()
    if accession in sampled_set and 'rep' in record.description.lower():
        selected_records.append(record)

# === STEP 3: Save output ===
SeqIO.write(selected_records, OUTPUT_FASTA, "fasta")

print(f"✅ Saved {len(selected_records)} Rep sequences to '{OUTPUT_FASTA}'")
