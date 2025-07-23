import pandas as pd
from Bio import SeqIO
import random

# === CONFIG ===
CSV_PATH = "Genomovirus_Proteins.csv"         # Metadata file
FASTA_PATH = "Genomovirus_Proteins.faa"       # Full FASTA file
OUTPUT_FASTA = "rep_subset_by_genus2.faa"      # Output file

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

# Make lookup tables for country and genus
metadata_lookup = df_filtered.set_index('Accession')[['Country', 'Genus']].to_dict(orient='index')

# === STEP 2: Filter FASTA and rename headers ===
selected_records = []
for record in SeqIO.parse(FASTA_PATH, "fasta"):
    accession = record.id.split('|')[0].strip()
    if accession in sampled_set and 'rep' in record.description.lower():
        if accession in metadata_lookup:
            # Fix NaN-safe handling
            country = metadata_lookup[accession]['Country']
            genus = metadata_lookup[accession]['Genus']
            country = str(country).strip().replace(" ", "_") if pd.notna(country) else "Unknown"
            genus = str(genus).strip().replace(" ", "_") if pd.notna(genus) else "Unknown"
            
            new_id = f"{accession}_{country}_{genus}"
            record.id = new_id
            record.name = new_id
            record.description = ""
            selected_records.append(record)

# === STEP 3: Save output ===
SeqIO.write(selected_records, OUTPUT_FASTA, "fasta")
print(f"✅ Saved {len(selected_records)} renamed Rep sequences to '{OUTPUT_FASTA}'")

