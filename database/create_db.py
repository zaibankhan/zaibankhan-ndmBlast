import sqlite3

conn = sqlite3.connect("database/ndm.db")

cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS ndm_variants (
    variant_id TEXT PRIMARY KEY,
    protein_name TEXT,
    dna_name TEXT,
    description TEXT,
    country TEXT,
    host_bacteria TEXT,
    publication TEXT,
    sequence TEXT
)
""")

conn.commit()
conn.close()

print("NDM Database Created Successfully")