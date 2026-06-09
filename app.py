from flask import Flask, render_template, request
from Bio import SeqIO
import re
import sqlite3

app = Flask(__name__)

DATABASE_NUCLEOTIDE = "database/NUCLEOTIDE_ndmfinal.txt"
DATABASE_PROTEIN = "database/PROTEIN_ndmfinal.txt"


# FASTA DATABASE LOAD KARNE KE LIYE
def load_database(file_path):

    records = []

    for record in SeqIO.parse(file_path, "fasta"):

        print("ID= =",record.id)

        records.append({
            "id": record.id,
            "description": record.description,
            "sequence": str(record.seq).upper()
        })

    return records


# SIMILARITY CALCULATE KARNE KE LIYE
def calculate_similarity(query, target):

    matches = 0

    for i in range(min(len(query), len(target))):

        if query[i] == target[i]:
            matches += 1

    length = max(len(query), len(target))

    if length == 0:
        return 0

    return round((matches / length) * 100, 2)


# DATABASE SEARCH
def search_database(query_sequence, database_file):

    database = load_database(database_file)

    results = []

    for entry in database:

        similarity = calculate_similarity(
            query_sequence,
            entry["sequence"]
        )
        
        results.append({
            "name": entry["id"],
            "description": entry["description"],
            "similarity": similarity,
            "sequence": entry["sequence"],
        })

    results.sort(
        key=lambda x: x["similarity"],
        reverse=True
    )

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

    results = search_database(
        sequence,
        DATABASE_NUCLEOTIDE
    )
    print("RESULTS =", results)
    print("FIRST NAME =", results[0]["name"])
    return render_template(
        "BlastN.html",
        results=results
    )


# BLAST P
@app.route("/runblastp", methods=["POST"])
def runblastp():

    sequence = request.form["sequence"]

    sequence = sequence.upper()

    sequence = re.sub(r"[^ARNDCQEGHILKMFPSTWYV]", "", sequence)

    results = search_database(
        sequence,
        DATABASE_PROTEIN
    )
    print("RESULTS =", results)
    return render_template(
        "BlastP.html",
        results=results
    )

@app.route("/details/<name>")
def details(name):

    conn = sqlite3.connect("database/ndm.db")
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(ndm_variants)")
    print(cursor.fetchall())
    cursor.execute("SELECT variant_id FROM ndm_variants")
    print(cursor.fetchall()[:5])

    print("NAME REPR=", repr(name))

    cursor.execute(
        "SELECT * FROM ndm_variants WHERE variant_id=?",
        (name,)
    )

    row = cursor.fetchone() 
    
    print("NAME =", name)
    print("ROW =", row)

    conn.close()

    if row:

        record = {
            "id": row[0],
            "protein_name": row[1],
            "dna_name": row[2],
            "description": row[3],
            "country": row[4],
            "host_bacteria": row[5],
            "publication": row[6],
            "sequence": row[7]
        }

        return render_template(
            "details.html",
            record=record
        )

    return "Record Not Found"
if __name__ == "__main__":
    app.run(debug=True)