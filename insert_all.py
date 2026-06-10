from Bio import SeqIO
import sqlite3

conn = sqlite3.connect("database/ndm_fixed.db")
cursor = conn.cursor()

# NUCLEOTIDE DATA — ndm_nucleotide table mein
n_count = 0
for record in SeqIO.parse("database/NUCLEOTIDE_ndmfinal.txt", "fasta"):
    try:
        cursor.execute(
            "INSERT OR REPLACE INTO ndm_nucleotide VALUES (?, ?, ?)",
            (record.id, record.description, str(record.seq).upper())
        )
        n_count += 1
    except Exception as e:
        print("NUC error:", e)

print(f"Nucleotide records inserted: {n_count}")

# PROTEIN DATA — ndm_protein table mein
p_count = 0
for record in SeqIO.parse("database/PROTEIN_ndmfinal.txt", "fasta"):
    try:
        cursor.execute(
            "INSERT OR REPLACE INTO ndm_protein VALUES (?, ?, ?)",
            (record.id, record.description, str(record.seq).upper())
        )
        p_count += 1
    except Exception as e:
        print("PRO error:", e)

print(f"Protein records inserted: {p_count}")

conn.commit()
conn.close()
print("Done! ndm_fixed.db ready.")
