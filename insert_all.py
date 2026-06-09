from Bio import SeqIO
import sqlite3

conn = sqlite3.connect("database/ndm.db")
cursor = conn.cursor()

# PROTEIN DATA INSERT
for record in SeqIO.parse("database/PROTEIN_ndmfinal.txt", "fasta"):

    try:
        cursor.execute("""
        INSERT OR IGNORE INTO ndm_variants
        VALUES (?, '', '', ?, '', '', 'https://pubmed.ncbi.nlm.nih.gov/', ?)
        """, (
            record.id,
            record.description,
            str(record.seq)
        ))
    except Exception as e:
        print(e)

# NUCLEOTIDE DATA INSERT
for record in SeqIO.parse("database/NUCLEOTIDE_ndmfinal.txt", "fasta"):

    try:
        cursor.execute("""
        INSERT OR IGNORE INTO ndm_variants
        VALUES (?, '', '', ?, '', '', 'https://pubmed.ncbi.nlm.nih.gov/', ?)
        """, (
            record.id,
            record.description,
            str(record.seq)
        ))
    except Exception as e:
        print(e)


conn.commit()
conn.close()

print("Done")