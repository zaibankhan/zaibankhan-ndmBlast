from flask import Flask, render_template, request, session, redirect, url_for
from Bio import SeqIO
import re
import sqlite3
import math


# ── NDM VARIANT STRUCTURE LINKS ───────────────────────────────────────────────
NDM_ALPHAFOLD_LINKS = {
    "NDM-1":  "https://alphafold.ebi.ac.uk/entry/AF-C7C422-F1",
    "NDM-2":  "https://alphafold.ebi.ac.uk/entry/AF-F2YZ26-F1",
    "NDM-3":  "https://alphafold.ebi.ac.uk/entry/AF-A0A2L0ARV7-F1",
    "NDM-4":  "https://alphafold.ebi.ac.uk/entry/AF-A0A8D6APR8-F1",
    "NDM-5":  "https://alphafold.ebi.ac.uk/entry/AF-A0A222U9D1-F1",
    "NDM-6":  "https://alphafold.ebi.ac.uk/entry/AF-A0A290DR99-F1",
    "NDM-7":  "https://alphafold.ebi.ac.uk/entry/AF-A0A1V0M4U5-F1",
    "NDM-8":  "https://alphafold.ebi.ac.uk/entry/AF-M1VE66-F1",
    "NDM-9":  "https://alphafold.ebi.ac.uk/entry/AF-A0A6G6ANE4-F1",
    "NDM-10": "https://alphafold.ebi.ac.uk/entry/AF-S5ZIP8-F1",
    "NDM-11": "https://alphafold.ebi.ac.uk/entry/AF-A0A7T1X4S0-F1",
    "NDM-12": "https://alphafold.ebi.ac.uk/search/text/blaNDM-12",
    "NDM-13": "https://alphafold.ebi.ac.uk/entry/AF-A0A0A8J940-F1",
    "NDM-14": "https://alphafold.ebi.ac.uk/entry/AF-A0A0C5H135-F1",
    "NDM-15": "https://alphafold.ebi.ac.uk/entry/AF-A0A0F6ZNP0-F1",
    "NDM-16": "https://alphafold.ebi.ac.uk/entry/AF-A0A286QUN1-F1",
    "NDM-17": "https://alphafold.ebi.ac.uk/entry/AF-A0A2P1H1P0-F1",
    "NDM-18": "https://alphafold.ebi.ac.uk/entry/AF-A0A1P8VH15-F1",
    "NDM-19": "https://alphafold.ebi.ac.uk/search/text/NDM-19",
    "NDM-21": "https://alphafold.ebi.ac.uk/entry/AF-A0A291NY14-F1",
    "NDM-22": "https://alphafold.ebi.ac.uk/entry/AF-A0A2S1T3V3-F1",
    "NDM-23": "https://alphafold.ebi.ac.uk/entry/AF-A0A2Z4BV56-F1",
    "NDM-24": "https://alphafold.ebi.ac.uk/entry/AF-A0A2Z4BV07-F1",
    "NDM-25": "https://alphafold.ebi.ac.uk/entry/AF-A0A5K6VNM3-F1",
    "NDM-26": "https://alphafold.ebi.ac.uk/entry/AF-A0A5K6W925-F1",
    "NDM-27": "https://alphafold.ebi.ac.uk/entry/AF-A0A3G3C0Q6-F1",
    "NDM-28": "https://alphafold.ebi.ac.uk/entry/AF-A0A410SN60-F1",
    "NDM-29": "https://alphafold.ebi.ac.uk/entry/AF-A0A5Q0MV96-F1",
    "NDM-30": "https://alphafold.ebi.ac.uk/search/text/NDM-30",
    "NDM-31": "https://alphafold.ebi.ac.uk/entry/AF-A0A7U3SV85-F1",
    "NDM-34": "https://alphafold.ebi.ac.uk/entry/AF-A0A8E7DAJ6-F1",
    "NDM-35": "https://alphafold.ebi.ac.uk/entry/AF-A0A8E7DD97-F1",
    "NDM-38": "https://alphafold.ebi.ac.uk/entry/AF-A0A8F1D5U2-F1",
    "NDM-39": "https://alphafold.ebi.ac.uk/entry/AF-A0A8G1A6Q8-F1",
    "NDM-40": "https://alphafold.ebi.ac.uk/entry/AF-A0A8G1A6Z8-F1",
}

NDM_PDB_LINKS = {
    "NDM-1": "https://www.rcsb.org/structure/9O2W",
    "NDM-4": "https://www.rcsb.org/structure/8SK2",
}

app = Flask(__name__)
app.secret_key = "ndm123"

DATABASE_NUCLEOTIDE = "database/NDM_genes.fasta"
DATABASE_PROTEIN    = "database/NDM_proteins.fasta"
DATABASE_FIXED      = "database/ndm_new.db"

# Protein-unique chars — ye DNA mein kabhi nahi aate
PROTEIN_ONLY_CHARS = set("RDEQHILKMFPWYV")


def load_database(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append({
            "id":          record.id,
            "description": record.description,
            "sequence":    str(record.seq).upper()
        })
    return records


def calculate_similarity(query, target):
    if not query or not target:
        return 0.0
    matches = sum(1 for a, b in zip(query, target) if a == b)
    length  = max(len(query), len(target))
    return round((matches / length) * 100, 2)


def calculate_evalue(similarity, query_length, db_size=97):
    if similarity >= 100.0:
        return 0.0
    mismatches = round((1.0 - similarity / 100.0) * query_length)
    if mismatches == 0:
        return 0.0
    try:
        evalue = (db_size * math.exp(mismatches * 0.05)) / (query_length * 10)
        return round(min(evalue, 10.0), 6)
    except Exception:
        return 10.0


def search_database(query_sequence, database_file):
    database = load_database(database_file)
    results  = []
    for entry in database:
        similarity = calculate_similarity(query_sequence, entry["sequence"])
        evalue     = calculate_evalue(similarity, len(query_sequence))
        results.append({
            "name":        entry["id"],
            "description": entry["description"],
            "similarity":  similarity,
            "evalue":      evalue,
            "sequence":    entry["sequence"],
        })
    results.sort(key=lambda x: x["similarity"], reverse=True)
    return results[:5]


def build_alignment(query_seq, subject_seq):
    max_len        = max(len(query_seq), len(subject_seq))
    query_padded   = query_seq.ljust(max_len, "-")
    subject_padded = subject_seq.ljust(max_len, "-")

    match_line = ""
    for q, s in zip(query_padded, subject_padded):
        if q == "-" and s == "-":
            match_line += " "
        elif q == "-" or s == "-":
            match_line += "-"
        elif q == s:
            match_line += "|"
        else:
            match_line += "*"

    chunk_size = 60
    blocks = []
    for i in range(0, max_len, chunk_size):
        q_chunk = query_padded[i:i+chunk_size]
        m_chunk = match_line[i:i+chunk_size]
        s_chunk = subject_padded[i:i+chunk_size]
        if q_chunk.strip("-") or s_chunk.strip("-"):
            blocks.append({
                "query":   q_chunk,
                "match":   m_chunk,
                "subject": s_chunk,
                "start":   i + 1,
            })
    return blocks, match_line


def is_protein_sequence(raw):
    """Return True agar protein-unique characters milein"""
    clean = re.sub(r"[^A-Z]", "", raw.upper())
    return any(c in PROTEIN_ONLY_CHARS for c in clean)


def is_dna_sequence(raw):
    """Return True agar 90%+ chars ATGCN hain aur koi protein char nahi"""
    clean = re.sub(r"[^A-Z]", "", raw.upper())
    if not clean:
        return False
    has_protein = any(c in PROTEIN_ONLY_CHARS for c in clean)
    dna_ratio   = len(re.sub(r"[^ATGCN]", "", clean)) / len(clean)
    return (not has_protein) and (dna_ratio > 0.90)


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
    # Fresh page load — purana error clear karo
    session["blast_error"]    = ""
    session["blastn_results"] = []
    return render_template("BlastN.html", results=[])

@app.route("/blastp")
def blastp():
    # Fresh page load — purana error clear karo
    session["blast_error"]    = ""
    session["blastp_results"] = []
    return render_template("BlastP.html", results=[])


@app.route("/runblastn", methods=["POST"])
def runblastn():
    raw = request.form["sequence"].strip()

    session["last_blast"]     = "n"
    session["blastn_results"] = []
    session["blast_error"]    = ""

    # Protein sequence BlastN mein dali?
    if is_protein_sequence(raw):
        session["blastn_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ Protein sequence detect hui! BlastN sirf DNA sequences (A, T, G, C) accept karta hai. Protein ke liye BlastP use karo."
        return redirect(url_for("blastn_results"))

    sequence = re.sub(r"[^ATGCN]", "", raw.upper())

    if len(sequence) < 50:
        session["blastn_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ Sequence bahut chhoti hai! Minimum 50 nucleotides required."
        return redirect(url_for("blastn_results"))

    session["blastn_query"]   = sequence
    session["query_sequence"] = sequence

    results = search_database(sequence, DATABASE_NUCLEOTIDE)
    session["blastn_results"] = results

    print("BlastN top hit:", results[0]["name"] if results else "No results")
    return redirect(url_for("blastn_results"))


@app.route("/runblastp", methods=["POST"])
def runblastp():
    raw = request.form["sequence"].strip()

    session["last_blast"]     = "p"
    session["blastp_results"] = []
    session["blast_error"]    = ""

    # DNA sequence BlastP mein dali?
    if is_dna_sequence(raw):
        session["blastp_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ DNA sequence detect hui! BlastP sirf Protein sequences (amino acids) accept karta hai. DNA ke liye BlastN use karo."
        return redirect(url_for("blastp_results"))

    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYVX]", "", raw.upper())

    if len(sequence) < 20:
        session["blastp_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ Sequence bahut chhoti hai! Minimum 20 amino acids required."
        return redirect(url_for("blastp_results"))

    session["blastp_query"]   = sequence
    session["query_sequence"] = sequence

    results = search_database(sequence, DATABASE_PROTEIN)
    session["blastp_results"] = results

    print("BlastP top hit:", results[0]["name"] if results else "No results")
    return redirect(url_for("blastp_results"))


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

    if blast_type == "n":
        query_seq   = session.get("blastn_query", "")
        all_results = session.get("blastn_results", [])
    else:
        query_seq   = session.get("blastp_query", "")
        all_results = session.get("blastp_results", [])

    conn = sqlite3.connect(DATABASE_FIXED)
    cursor = conn.cursor()
    if blast_type == "n":
        cursor.execute("SELECT variant_id, description, sequence FROM ndm_nucleotide WHERE variant_id=?", (name,))
    else:
        cursor.execute("SELECT variant_id, description, sequence FROM ndm_protein WHERE variant_id=?", (name,))
    seq_row = cursor.fetchone()

    meta_row = None
    try:
        cursor.execute("SELECT * FROM ndm_variants WHERE variant_id=?", (name,))
        meta_row = cursor.fetchone()
    except Exception:
        pass
    conn.close()

    subject_seq = seq_row[2] if seq_row else ""

    record = {
        "id":                  name,
        "description":         seq_row[1] if seq_row else "",
        "protein_name":        meta_row[1] if meta_row and len(meta_row) > 1 else "",
        "dna_name":            meta_row[2] if meta_row and len(meta_row) > 2 else "",
        "country":             meta_row[4] if meta_row and len(meta_row) > 4 else "",
        "host_bacteria":       meta_row[5] if meta_row and len(meta_row) > 5 else "",
        "publication":         meta_row[6] if meta_row and len(meta_row) > 6 else "",
        "sequence":            subject_seq,
        "detection_technique": meta_row[8] if meta_row and len(meta_row) > 8 else "N/A"
    }

    hit = next((r for r in all_results if r["name"] == name), None)
    alignment_blocks, match_line = build_alignment(query_seq, subject_seq)

    alphafold_url = NDM_ALPHAFOLD_LINKS.get(name, "")
    pdb_url       = NDM_PDB_LINKS.get(name, "")

    return render_template(
        "details.html",
        record=record,
        query_seq=query_seq,
        subject_seq=subject_seq,
        match_line=match_line,
        query_length=len(query_seq),
        subject_length=len(subject_seq),
        blast_type=blast_type,
        alignment_blocks=alignment_blocks,
        hit=hit,
        alphafold_url=alphafold_url,
        pdb_url=pdb_url,
    )


if __name__ == "__main__":
    app.run(debug=True)
