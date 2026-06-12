from flask import Flask, render_template, request, session, redirect, url_for
from Bio import SeqIO
from Bio.Blast import NCBIXML
import re
import sqlite3
import tempfile
import os
import subprocess

app = Flask(__name__)
app.secret_key = "ndm123"

BLAST_NUCL_DB = "database/blast_db/ndm_nucl"
BLAST_PROT_DB = "database/blast_db/ndm_prot"
DATABASE_FIXED = "database/ndm_fixed.db"


def run_blast(query_sequence, blast_type="n", top_n=5):
    """
    subprocess se local BLAST run karo — Bio.Blast.Applications nahi chahiye.
    Sirf local database use hogi — internet se kuch nahi.
    """
    results = []

    # Query temp FASTA file mein likho
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                     delete=False) as qf:
        qf.write(f">query\n{query_sequence}\n")
        query_file = qf.name

    out_file = query_file + "_blast.xml"

    try:
        if blast_type == "n":
            cmd = [
                "blastn",
                "-query", query_file,
                "-db", BLAST_NUCL_DB,
                "-outfmt", "5",
                "-out", out_file,
                "-max_target_seqs", str(top_n),
                "-evalue", "10",
                "-word_size", "7",
                "-dust", "no",
            ]
        else:
            cmd = [
                "blastp",
                "-query", query_file,
                "-db", BLAST_PROT_DB,
                "-outfmt", "5",
                "-out", out_file,
                "-max_target_seqs", str(top_n),
                "-evalue", "10",
                "-word_size", "2",
            ]

        # BLAST run karo locally
        proc = subprocess.run(cmd, capture_output=True, text=True)

        if proc.returncode != 0:
            print("BLAST stderr:", proc.stderr)
            return results

        # XML parse karo
        with open(out_file) as result_handle:
            blast_records = list(NCBIXML.parse(result_handle))

        if blast_records and blast_records[0].alignments:
            for alignment in blast_records[0].alignments[:top_n]:
                hsp = alignment.hsps[0]

                identity = round((hsp.identities / hsp.align_length) * 100, 2)
                evalue   = hsp.expect

                subject_id = alignment.hit_id
                if "|" in subject_id:
                    subject_id = subject_id.split("|")[-1]
                if not subject_id or subject_id == "N/A":
                    subject_id = alignment.hit_def.split()[0]

                results.append({
                    "name":        subject_id,
                    "description": alignment.hit_def,
                    "similarity":  identity,
                    "evalue":      evalue,
                    "sequence":    hsp.sbjct.replace("-", ""),
                    "query_hsp":   hsp.query,
                    "subject_hsp": hsp.sbjct,
                    "match_hsp":   hsp.match,
                    "q_start":     hsp.query_start,
                    "s_start":     hsp.sbjct_start,
                    "align_len":   hsp.align_length,
                })

    except Exception as e:
        print("BLAST error:", e)

    finally:
        if os.path.exists(query_file): os.remove(query_file)
        if os.path.exists(out_file):   os.remove(out_file)

    results.sort(key=lambda x: x["similarity"], reverse=True)
    return results


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

    results = run_blast(sequence, blast_type="n")
    session["blastn_results"] = results

    print("BlastN hits:", len(results))
    if results:
        print("Top:", results[0]["name"], results[0]["similarity"], "%", "E:", results[0]["evalue"])

    return redirect(url_for("blastn_results"))


@app.route("/runblastp", methods=["POST"])
def runblastp():
    sequence = request.form["sequence"].upper()
    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYVX]", "", sequence)

    session["blastp_query"]   = sequence
    session["query_sequence"] = sequence
    session["last_blast"]     = "p"

    results = run_blast(sequence, blast_type="p")
    session["blastp_results"] = results

    print("BlastP hits:", len(results))
    if results:
        print("Top:", results[0]["name"], results[0]["similarity"], "%", "E:", results[0]["evalue"])

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
        all_results = session.get("blastn_results", [])
        query_seq   = session.get("blastn_query", "")
    else:
        all_results = session.get("blastp_results", [])
        query_seq   = session.get("blastp_query", "")

    hit = next((r for r in all_results if r["name"] == name), None)

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

    subject_seq = seq_row[2] if seq_row else (hit["sequence"] if hit else "")

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

    if hit and "query_hsp" in hit:
        query_aligned   = hit["query_hsp"]
        subject_aligned = hit["subject_hsp"]
        match_str       = hit["match_hsp"]
        q_start         = hit["q_start"]
        s_start         = hit["s_start"]
    else:
        max_len         = max(len(query_seq), len(subject_seq))
        query_aligned   = query_seq.ljust(max_len, "-")
        subject_aligned = subject_seq.ljust(max_len, "-")
        match_str       = "".join(
            "|" if a == b and a != "-" else " "
            for a, b in zip(query_aligned, subject_aligned)
        )
        q_start = s_start = 1

    chunk_size = 60
    alignment_blocks = []
    align_len = max(len(query_aligned), len(subject_aligned))
    for i in range(0, align_len, chunk_size):
        q_chunk = query_aligned[i:i+chunk_size]
        m_chunk = match_str[i:i+chunk_size]
        s_chunk = subject_aligned[i:i+chunk_size]
        if q_chunk.strip("-") or s_chunk.strip("-"):
            alignment_blocks.append({
                "query":   q_chunk,
                "match":   m_chunk,
                "subject": s_chunk,
                "q_pos":   q_start + i,
                "s_pos":   s_start + i,
            })

    return render_template(
        "details.html",
        record=record,
        query_seq=query_seq,
        subject_seq=subject_seq,
        query_length=len(query_seq),
        subject_length=len(subject_seq),
        blast_type=blast_type,
        alignment_blocks=alignment_blocks,
        hit=hit,
    )


if __name__ == "__main__":
    app.run(debug=True)
