from flask import Flask, render_template, request, session, redirect, url_for
from Bio import SeqIO
import re
import sqlite3
import math

app = Flask(__name__)
app.secret_key = "ndm123"

DATABASE_NUCLEOTIDE = "database/NUCLEOTIDE_ndmfinal.txt"
DATABASE_PROTEIN    = "database/PROTEIN_ndmfinal.txt"
DATABASE_FIXED      = "database/ndm_fixed.db"


# ── DATABASE LOAD ─────────────────────────────────────────────────────────────
def load_database(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append({
            "id":          record.id,
            "description": record.description,
            "sequence":    str(record.seq).upper()
        })
    return records


# ── SIMILARITY — Direct positional comparison (most accurate for NDM variants)
# Aapka data same-length sequences hai (810-828bp), isliye direct comparison best hai
def calculate_similarity(query, target):
    if not query or not target:
        return 0.0

    # Matched positions / longer sequence length
    matches = sum(1 for a, b in zip(query, target) if a == b)
    length  = max(len(query), len(target))
    return round((matches / length) * 100, 2)


# ── E-VALUE — Mismatch based meaningful formula ──────────────────────────────
def calculate_evalue(similarity, query_length, db_size=67):
    """
    similarity 100%  → E-value 0.0       (perfect match)
    similarity 99%   → E-value ~0.0001   (very significant)
    similarity 80%   → E-value ~0.001    (significant)
    similarity 42%   → E-value ~0.75     (low significance)
    similarity <30%  → E-value ~10       (not significant)
    """
    if similarity >= 100.0:
        return 0.0
    if similarity <= 0.0:
        return 10.0

    mismatch_fraction = 1.0 - (similarity / 100.0)
    norm_len = query_length / 100.0

    try:
        evalue = db_size * (mismatch_fraction ** norm_len)
        return round(min(evalue, 10.0), 4)
    except Exception:
        return 10.0


# ── SEARCH DATABASE ───────────────────────────────────────────────────────────
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


# ── ALIGNMENT BUILDER ─────────────────────────────────────────────────────────
def build_alignment(query_seq, subject_seq):
    """
    Poori query aur poora subject align karo.
    Chhoti sequence ko '-' se pad karo.
    """
    max_len         = max(len(query_seq), len(subject_seq))
    query_padded    = query_seq.ljust(max_len, "-")
    subject_padded  = subject_seq.ljust(max_len, "-")

    match_line = ""
    for q, s in zip(query_padded, subject_padded):
        if q == "-" or s == "-":
            match_line += " "   # gap
        elif q == s:
            match_line += "|"   # exact match
        else:
            match_line += "."   # mismatch

    # 60-char chunks mein tod do
    chunk_size = 60
    blocks     = []
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


# ── ROUTES ────────────────────────────────────────────────────────────────────
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


@app.route("/runblastn", methods=["POST"])
def runblastn():
    sequence = request.form["sequence"].upper()
    sequence = re.sub(r"[^ATGCN]", "", sequence)

    session["blastn_query"]   = sequence
    session["query_sequence"] = sequence
    session["last_blast"]     = "n"

    results = search_database(sequence, DATABASE_NUCLEOTIDE)
    session["blastn_results"] = results

    print("BlastN top hit:", results[0]["name"], results[0]["similarity"], "%" if results else "No results")
    return redirect(url_for("blastn_results"))


@app.route("/runblastp", methods=["POST"])
def runblastp():
    sequence = request.form["sequence"].upper()
    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYVX]", "", sequence)

    session["blastp_query"]   = sequence
    session["query_sequence"] = sequence
    session["last_blast"]     = "p"

    results = search_database(sequence, DATABASE_PROTEIN)
    session["blastp_results"] = results

    print("BlastP top hit:", results[0]["name"], results[0]["similarity"], "%" if results else "No results")
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

    # DB se sahi sequence lo
    conn = sqlite3.connect(DATABASE_FIXED)
    cursor = conn.cursor()
    if blast_type == "n":
        cursor.execute(
            "SELECT variant_id, description, sequence FROM ndm_nucleotide WHERE variant_id=?",
            (name,)
        )
    else:
        cursor.execute(
            "SELECT variant_id, description, sequence FROM ndm_protein WHERE variant_id=?",
            (name,)
        )
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

    # Result stats
    hit = next((r for r in all_results if r["name"] == name), None)

    # Alignment build karo
    alignment_blocks, match_line = build_alignment(query_seq, subject_seq)

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
    )


if __name__ == "__main__":
    app.run(debug=True)
