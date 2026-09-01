from flask import Flask, render_template, request, session, redirect, url_for
from Bio import SeqIO
import re
import sqlite3
import math
import json
import os
from datetime import datetime


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

# NDM variant ke DOI/publication links — jis NDM ki jitni links hain,
# sab list mein aayengi (excel se banaya gaya hai)
PUBLICATIONS_FILE = "database/ndm_publications.json"
NDM_PUBLICATIONS = {}
if os.path.exists(PUBLICATIONS_FILE):
    with open(PUBLICATIONS_FILE, "r", encoding="utf-8") as f:
        NDM_PUBLICATIONS = json.load(f)

# Protein-unique chars — ye DNA mein kabhi nahi aate
PROTEIN_ONLY_CHARS = set("RDEQHILKMFPWYV")

# ── SEARCH HISTORY (persisted locally, never sent online) ───────────────────
HISTORY_FILE = "database/search_history.json"
HISTORY_LIMIT = 100


def load_history():
    """Read persisted search history from the local JSON file."""
    if not os.path.exists(HISTORY_FILE):
        return []
    try:
        with open(HISTORY_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
            return data if isinstance(data, list) else []
    except Exception:
        return []


def save_history(records):
    try:
        with open(HISTORY_FILE, "w", encoding="utf-8") as f:
            json.dump(records, f, indent=2)
    except Exception:
        pass


def add_history_entry(entry):
    """Prepend a new search entry, keeping only the newest HISTORY_LIMIT."""
    records = load_history()
    records.insert(0, entry)
    save_history(records[:HISTORY_LIMIT])


def clear_history():
    save_history([])


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
    """Simple/global method: compare both sequences from offset 0 (no gaps)."""
    if not query or not target:
        return 0.0
    matches = sum(1 for a, b in zip(query, target) if a == b)
    length  = max(len(query), len(target))
    return round((matches / length) * 100, 2)


def smith_waterman_align(query, subject, match=2, mismatch=-1, gap_open=-2, gap_extend=-1):
    """Gapped local alignment (Smith-Waterman, affine gap penalty).
    Used to refine the region found by the BLAST seed-and-extend step.

    Returns (q_ali, s_ali, match_line, score, identity_pct, aln_len).
    """
    n, m = len(query), len(subject)
    if n == 0 or m == 0:
        return ("", "", "", 0, 0.0, 0)

    # Affine-gap dynamic programming matrices
    H = [[0] * (m + 1) for _ in range(n + 1)]  # best score ending at (i,j)
    E = [[0] * (m + 1) for _ in range(n + 1)]  # gap in query (moving down)
    F = [[0] * (m + 1) for _ in range(n + 1)]  # gap in subject (moving right)
    T = [[0] * (m + 1) for _ in range(n + 1)]  # 0=stop,1=diag,2=E,3=F

    max_score = 0
    max_i = max_j = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            E[i][j] = max(H[i-1][j] + gap_open, E[i-1][j] + gap_extend)
            F[i][j] = max(H[i][j-1] + gap_open, F[i][j-1] + gap_extend)
            diag = H[i-1][j-1] + (match if query[i-1] == subject[j-1] else mismatch)
            best = max(0, diag, E[i][j], F[i][j])
            H[i][j] = best
            if best == diag:
                T[i][j] = 1
            elif best == E[i][j]:
                T[i][j] = 2
            elif best == F[i][j]:
                T[i][j] = 3
            else:
                T[i][j] = 0
            if best > max_score:
                max_score, max_i, max_j = best, i, j

    # Traceback
    i, j = max_i, max_j
    q_ali, s_ali = [], []
    while i > 0 and j > 0 and T[i][j] != 0:
        if T[i][j] == 1:
            q_ali.append(query[i-1]); s_ali.append(subject[j-1]); i -= 1; j -= 1
        elif T[i][j] == 2:
            q_ali.append(query[i-1]); s_ali.append("-"); i -= 1
        else:
            q_ali.append("-"); s_ali.append(subject[j-1]); j -= 1

    q_ali.reverse(); s_ali.reverse()
    q_ali_str = "".join(q_ali)
    s_ali_str = "".join(s_ali)

    aln_len  = len(q_ali_str)
    matches  = sum(1 for a, b in zip(q_ali_str, s_ali_str) if a == b and a != "-")
    identity = round((matches / aln_len) * 100, 2) if aln_len else 0.0

    match_line = ""
    for a, b in zip(q_ali_str, s_ali_str):
        if a == "-" or b == "-":
            match_line += "-"
        elif a == b:
            match_line += "|"
        else:
            match_line += "*"

    return q_ali_str, s_ali_str, match_line, max_score, identity, aln_len


def build_word_index(seq, w):
    """BLAST step 1 (seeding): map every word (k-mer) of the query to its positions."""
    index = {}
    for i in range(0, len(seq) - w + 1):
        index.setdefault(seq[i:i+w], []).append(i)
    return index


def blast_align(query, subject, word_size=11, match=2, mismatch=-1,
                gap_open=-2, gap_extend=-1, x_drop=5):
    """BLAST-style heuristic local alignment:
      1. Seed: exact word (k-mer) matches between query and subject.
      2. Extend: ungapped extension of each seed into a maximal HSP
         (High-scoring Segment Pair) with an X-drop threshold.
      3. Refine: gapped Smith-Waterman within the best HSP window.

    Returns (q_ali, s_ali, match_line, score, identity_pct, aln_len).
    """
    n, m = len(query), len(subject)
    if n == 0 or m == 0:
        return ("", "", "", 0, 0.0, 0)

    # For very short inputs fall back to full gapped alignment.
    if n < word_size or m < word_size:
        return smith_waterman_align(query, subject, match, mismatch, gap_open, gap_extend)

    q_index = build_word_index(query, word_size)
    best_hsp = None  # (score, q_start, q_end, s_start, s_end)

    # Step 1 + 2: seed and ungapped extend
    for sj in range(0, m - word_size + 1):
        word = subject[sj:sj+word_size]
        for qi in q_index.get(word, []):
            si = sj
            # score of the seed itself
            seed = word_size * match

            # --- extend right ---
            score = seed
            best_score = seed
            best_off = 0
            off = 0
            while qi + word_size + off < n and si + word_size + off < m:
                score += match if query[qi+word_size+off] == subject[si+word_size+off] else mismatch
                off += 1
                if score > best_score:
                    best_score = score
                    best_off = off
                if score < best_score - x_drop:
                    break
            right_end_q = qi + word_size + best_off
            right_end_s = si + word_size + best_off
            right_score = best_score

            # --- extend left ---
            score = right_score
            best_score_left = right_score
            best_off_left = 0
            off = 0
            while qi - 1 - off >= 0 and si - 1 - off >= 0:
                score += match if query[qi-1-off] == subject[si-1-off] else mismatch
                off += 1
                if score > best_score_left:
                    best_score_left = score
                    best_off_left = off
                if score < best_score_left - x_drop:
                    break
            q_start = qi - best_off_left
            s_start = si - best_off_left
            hsp = (best_score_left, q_start, right_end_q, s_start, right_end_s)
            if best_hsp is None or hsp[0] > best_hsp[0]:
                best_hsp = hsp

    if best_hsp is None or best_hsp[0] <= 0:
        # No seed found: fall back to gapped local alignment across whole seqs.
        return smith_waterman_align(query, subject, match, mismatch, gap_open, gap_extend)

    _, q_start, q_end, s_start, s_end = best_hsp

    # Slight padding around the HSP before gapped refinement.
    pad = 20
    q_win_start = max(0, q_start - pad)
    q_win_end   = min(n, q_end + pad)
    s_win_start = max(0, s_start - pad)
    s_win_end   = min(m, s_end + pad)
    q_win = query[q_win_start:q_win_end]
    s_win = subject[s_win_start:s_win_end]

    qa_win, sa_win, ml_win, score, identity, aln_len = smith_waterman_align(
        q_win, s_win, match, mismatch, gap_open, gap_extend
    )

    if aln_len == 0:
        return ("", "", "", 0, 0.0, 0)

    return qa_win, sa_win, ml_win, score, identity, aln_len


def calculate_evalue(score, query_length, db_size=97, match=2):
    """Score-based expected-value estimate for the BLAST algorithm."""
    if score <= 0:
        return 1.0
    try:
        lambda_ = 2 * math.log(max(match, 2)) / 1.0  # ~per-residue scaling
        evalue = (db_size * (query_length / 100.0)) * math.exp(-lambda_ * score)
        return round(min(max(evalue, 1e-300), 10.0), 6)
    except Exception:
        return 10.0


def search_database(query_sequence, database_file, blast_type="n"):
    """Run the BLAST algorithm against a local database and return top hits."""
    word_size = 11 if blast_type == "n" else 3
    database = load_database(database_file)
    results  = []
    for entry in database:
        qa, sa, ml, score, identity, aln_len = blast_align(
            query_sequence, entry["sequence"], word_size=word_size
        )
        evalue = calculate_evalue(score, len(query_sequence))
        results.append({
            "name":        entry["id"],
            "description": entry["description"],
            "similarity":  identity,
            "evalue":      evalue,
            "sequence":    entry["sequence"],
        })
    results.sort(key=lambda x: x["similarity"], reverse=True)
    return results[:5]


def build_alignment(query_seq, subject_seq, blast_type="n"):
    """Build alignment display blocks using the BLAST algorithm."""
    word_size = 11 if blast_type == "n" else 3
    qa, sa, ml, score = blast_align(query_seq, subject_seq, word_size=word_size)[:4]
    max_len = len(qa)

    chunk_size = 60
    blocks = []
    for i in range(0, max_len, chunk_size):
        q_chunk = qa[i:i+chunk_size]
        m_chunk = ml[i:i+chunk_size]
        s_chunk = sa[i:i+chunk_size]
        if q_chunk.strip("-") or s_chunk.strip("-"):
            blocks.append({
                "query":   q_chunk,
                "match":   m_chunk,
                "subject": s_chunk,
                "start":   i + 1,
            })
    return blocks, ml


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
    return render_template("BlastN.html", results=[], query_sequence="")

@app.route("/blastp")
def blastp():
    # Fresh page load — purana error clear karo
    session["blast_error"]    = ""
    session["blastp_results"] = []
    return render_template("BlastP.html", results=[], query_sequence="")


@app.route("/runblastn", methods=["POST"])
def runblastn():
    raw = request.form["sequence"].strip()

    session["last_blast"]       = "n"
    session["blastn_results"]   = []
    session["blast_error"]      = ""

    # Protein sequence BlastN mein dali?
    if is_protein_sequence(raw):
        session["blastn_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ You have entered the wrong sequence. A protein sequence was detected, but BLASTn accepts only DNA nucleotides (A, T, G, C). Please use BLASTp for protein sequences."
        return redirect(url_for("blastn_results"))

    sequence = re.sub(r"[^ATGCN]", "", raw.upper())

    if len(sequence) < 20:
        session["blastn_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ Insufficient Sequence Length: Please enter a nucleotide sequence containing at least 20 bp to proceed with the analysis."
        return redirect(url_for("blastn_results"))

    session["blastn_query"]   = sequence
    session["query_sequence"] = sequence

    results = search_database(sequence, DATABASE_NUCLEOTIDE, "n")
    session["blastn_results"] = results

    top = results[0] if results else None
    add_history_entry({
        "time":       datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "blast_type": "BLASTn (DNA)",
        "method":     "BLAST algorithm",
        "method_key": "n",
        "query":      sequence,
        "query_len":  len(sequence),
        "top_hit":    top["name"] if top else "No match",
        "identity":   top["similarity"] if top else None,
        "evalue":     top["evalue"] if top else None,
    })

    print("BlastN top hit:", results[0]["name"] if results else "No results")
    return redirect(url_for("blastn_results", seq=sequence))


@app.route("/runblastp", methods=["POST"])
def runblastp():
    raw = request.form["sequence"].strip()

    session["last_blast"]       = "p"
    session["blastp_results"]   = []
    session["blast_error"]      = ""

    # DNA sequence BlastP mein dali?
    if is_dna_sequence(raw):
        session["blastp_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ You have entered the wrong sequence. A DNA/nucleotide sequence was detected, but BLASTp accepts only protein sequences (amino acids). Please use BLASTn for DNA sequences."
        return redirect(url_for("blastp_results"))

    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYVX]", "", raw.upper())

    if len(sequence) < 20:
        session["blastp_query"]   = ""
        session["query_sequence"] = ""
        session["blast_error"]    = "⚠️ Sequence bahut chhoti hai! Minimum 20 amino acids required."
        return redirect(url_for("blastp_results"))

    session["blastp_query"]   = sequence
    session["query_sequence"] = sequence

    results = search_database(sequence, DATABASE_PROTEIN, "p")
    session["blastp_results"] = results

    top = results[0] if results else None
    add_history_entry({
        "time":       datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "blast_type": "BLASTp (Protein)",
        "method":     "BLAST algorithm",
        "method_key": "p",
        "query":      sequence,
        "query_len":  len(sequence),
        "top_hit":    top["name"] if top else "No match",
        "identity":   top["similarity"] if top else None,
        "evalue":     top["evalue"] if top else None,
    })

    print("BlastP top hit:", results[0]["name"] if results else "No results")
    return redirect(url_for("blastp_results", seq=sequence))


@app.route("/blastp_results")
def blastp_results():
    seq_param = request.args.get("seq", "")
    if seq_param:
        sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYVX]", "", seq_param.upper())
        results = search_database(sequence, DATABASE_PROTEIN, "p")
        session["last_blast"]     = "p"
        session["blastp_query"]   = sequence
        session["query_sequence"] = sequence
        session["blastp_results"] = results
    else:
        sequence = session.get("blastp_query", "")
        results  = session.get("blastp_results", [])
    return render_template("BlastP.html", results=results, query_sequence=sequence)

@app.route("/blastn_results")
def blastn_results():
    seq_param = request.args.get("seq", "")
    if seq_param:
        sequence = re.sub(r"[^ATGCN]", "", seq_param.upper())
        results = search_database(sequence, DATABASE_NUCLEOTIDE, "n")
        session["last_blast"]     = "n"
        session["blastn_query"]   = sequence
        session["query_sequence"] = sequence
        session["blastn_results"] = results
    else:
        sequence = session.get("blastn_query", "")
        results  = session.get("blastn_results", [])
    return render_template("BlastN.html", results=results, query_sequence=sequence)


@app.route("/details/<name>")
def details(name):
    blast_type = request.args.get("type") or session.get("last_blast", "p")
    query_seq  = request.args.get("seq", "")

    if not query_seq:
        if blast_type == "n":
            query_seq = session.get("blastn_query", "")
        else:
            query_seq = session.get("blastp_query", "")

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

    hit = None
    if query_seq and subject_seq:
        qa, sa, ml, score, identity, aln_len = blast_align(
            query_seq, subject_seq, word_size=(11 if blast_type == "n" else 3)
        )
        evalue = calculate_evalue(score, len(query_seq))
        hit = {"similarity": identity, "evalue": evalue, "aln_len": aln_len}

    alignment_blocks, match_line = build_alignment(query_seq, subject_seq, blast_type)

    alphafold_url = NDM_ALPHAFOLD_LINKS.get(name, "")
    pdb_url       = NDM_PDB_LINKS.get(name, "")
    publication_links = NDM_PUBLICATIONS.get(name, [])

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
        publication_links=publication_links,
    )


@app.route("/publications/<name>")
def publications(name):
    links = NDM_PUBLICATIONS.get(name, [])
    return render_template(
        "publications.html",
        variant_id=name,
        links=links,
    )


@app.route("/history")
def history():
    records = load_history()
    return render_template("history.html", records=records, count=len(records))


@app.route("/clear_history", methods=["POST"])
def clear_history_route():
    clear_history()
    return redirect(url_for("history"))


if __name__ == "__main__":
    app.run(debug=True)
