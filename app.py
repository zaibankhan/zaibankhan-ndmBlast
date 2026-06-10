from flask import Flask, render_template, request, session
from Bio import SeqIO
import re
import sqlite3
from difflib import SequenceMatcher

app = Flask(__name__)
app.secret_key = "ndm123"

DATABASE_NUCLEOTIDE = "database/NUCLEOTIDE_ndmfinal.txt"
DATABASE_PROTEIN = "database/PROTEIN_ndmfinal.txt"
DATABASE_FIXED = "database/ndm_fixed.db"


# FASTA DATABASE LOAD KARNE KE LIYE
def load_database(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append({
            "id": record.id,
            "description": record.description,
            "sequence": str(record.seq).upper()
        })
    return records


# SIMILARITY CALCULATE KARNE KE LIYE — SequenceMatcher use karo
def calculate_similarity(query, target):
    if not query or not target:
        return 0
    matcher = SequenceMatcher(None, query, target)
    similarity = matcher.ratio() * 100
    return round(similarity, 2)


# E-VALUE — Biology-standard approximate formula
def calculate_evalue(similarity, query_length, db_size=1000000):
    if similarity >= 100:
        return 0.0
    identity_fraction = similarity / 100.0
    evalue = db_size * query_length * (1 - identity_fraction) ** query_length
    if evalue > 10:
        evalue = 10.0
    return round(evalue, 6)


# DATABASE SEARCH
def search_database(query_sequence, database_file):
    database = load_database(database_file)
    results = []
    for entry in database:
        similarity = calculate_similarity(query_sequence, entry["sequence"])
        evalue = calculate_evalue(similarity, len(query_sequence))
        results.append({
            "name": entry["id"],
            "description": entry["description"],
            "similarity": similarity,
            "evalue": evalue,
            "sequence": entry["sequence"],
        })
    results.sort(key=lambda x: x["similarity"], reverse=True)
    return results[:5]


@app.route("/")
def home():
    return render_template("HomeInfo.html")


@app.route("/dn")
def dn():
    return render_template("DNinformation.html")


@app.route("/blastinfo")
def blastinfo():
    return render_template("BlastInfo.html")


@app.route("/blastn")
def blastn():
    return render_template("BlastN.html")


@app.route("/blastp")
def blastp():
    return render_template("BlastP.html")


# BLAST N
@app.route("/runblastn", methods=["POST"])
def runblastn():
    sequence = request.form["sequence"]
    sequence = sequence.upper()
    sequence = re.sub(r"[^ATGC]", "", sequence)
    session["query_sequence"] = sequence
    results = search_database(sequence, DATABASE_NUCLEOTIDE)
    session["blastn_results"] = results
    session["last_blast"] = "n"
    print("RESULTS =", results)
    print("FIRST NAME =", results[0]["name"] if results else "No results")
    return render_template("BlastN.html", results=results)


# BLAST P
@app.route("/runblastp", methods=["POST"])
def runblastp():
    sequence = request.form["sequence"]
    sequence = sequence.upper()
    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYV]", "", sequence)
    session["query_sequence"] = sequence
    results = search_database(sequence, DATABASE_PROTEIN)
    session["blastp_results"] = results
    session["last_blast"] = "p"
    print("RESULTS =", results)
    return render_template("BlastP.html", results=results)


@app.route("/blastp_results")
def blastp_results():
    results = session.get("blastp_results", [])
    return render_template("BlastP.html", results=results)


@app.route("/blastn_results")
def blastn_results():
    results = session.get("blastn_results", [])
    return render_template("BlastN.html", results=results)


@app.route("/details/<name>")
def details(name):

    blast_type = session.get("last_blast", "p")

    conn = sqlite3.connect(DATABASE_FIXED)
    cursor = conn.cursor()

    # *** KEY FIX: blast_type ke hisaab se SAHI TABLE se sequence lo ***
    # BlastN → ndm_nucleotide table (DNA sequences)
    # BlastP → ndm_protein table (Protein sequences)
    if blast_type == "n":
        table = "ndm_nucleotide"
    else:
        table = "ndm_protein"

    cursor.execute(
        f"SELECT variant_id, description, sequence FROM {table} WHERE variant_id=?",
        (name,)
    )
    seq_row = cursor.fetchone()

    # Metadata (country, host, publication etc.) ndm_variants se lo agar ho
    meta_row = None
    try:
        cursor.execute(
            "SELECT * FROM ndm_variants WHERE variant_id=?",
            (name,)
        )
        meta_row = cursor.fetchone()
    except Exception:
        pass

    conn.close()

    print("NAME =", name, "| blast_type =", blast_type, "| table =", table)
    print("seq_row found:", seq_row is not None)

    if seq_row:

        subject_seq = seq_row[2]  # Sahi table se sahi sequence

        record = {
            "id": seq_row[0],
            "description": seq_row[1],
            "protein_name": meta_row[1] if meta_row and len(meta_row) > 1 else "",
            "dna_name": meta_row[2] if meta_row and len(meta_row) > 2 else "",
            "country": meta_row[4] if meta_row and len(meta_row) > 4 else "",
            "host_bacteria": meta_row[5] if meta_row and len(meta_row) > 5 else "",
            "publication": meta_row[6] if meta_row and len(meta_row) > 6 else "",
            "sequence": subject_seq,
            "detection_technique": meta_row[8] if meta_row and len(meta_row) > 8 else "N/A"
        }

        query_seq = session.get("query_sequence", "")

        query_length = len(query_seq)
        subject_length = len(subject_seq)
        max_len = max(query_length, subject_length)

        # Query aur subject ko equal length karo — gap character '-' se
        query_padded = query_seq.ljust(max_len, '-')
        subject_padded = subject_seq.ljust(max_len, '-')

        # Match line — poori length pe banao
        match_line = ""
        for q, s in zip(query_padded, subject_padded):
            if q == s and q != '-':
                match_line += "|"
            else:
                match_line += " "

        # Alignment ko 60-character chunks mein tod do readability ke liye
        chunk_size = 60
        alignment_blocks = []
        for i in range(0, max_len, chunk_size):
            alignment_blocks.append({
                "query": query_padded[i:i+chunk_size],
                "match": match_line[i:i+chunk_size],
                "subject": subject_padded[i:i+chunk_size],
                "start": i + 1
            })

        return render_template(
            "details.html",
            record=record,
            query_seq=query_seq,
            subject_seq=subject_seq,
            match_line=match_line,
            query_length=query_length,
            subject_length=subject_length,
            blast_type=blast_type,
            alignment_blocks=alignment_blocks
        )

    return "Record Not Found"


if __name__ == "__main__":
    app.run(debug=True)
