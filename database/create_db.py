import sqlite3

conn = sqlite3.connect("database/ndm_fixed.db")
cursor = conn.cursor()

# Nucleotide sequences ke liye alag table
cursor.execute("""
CREATE TABLE IF NOT EXISTS ndm_nucleotide (
    variant_id   TEXT PRIMARY KEY,
    description  TEXT,
    sequence     TEXT
)
""")

# Protein sequences ke liye alag table
cursor.execute("""
CREATE TABLE IF NOT EXISTS ndm_protein (
    variant_id   TEXT PRIMARY KEY,
    description  TEXT,
    sequence     TEXT
)
""")

# Metadata ke liye table (country, host, publication etc.)
cursor.execute("""
CREATE TABLE IF NOT EXISTS ndm_variants (
    variant_id        TEXT PRIMARY KEY,
    protein_name      TEXT,
    dna_name          TEXT,
    description       TEXT,
    country           TEXT,
    host_bacteria     TEXT,
    publication       TEXT,
    sequence          TEXT
)
""")

conn.commit()
conn.close()
print("NDM Database Created Successfully (2 sequence tables + 1 metadata table)")
