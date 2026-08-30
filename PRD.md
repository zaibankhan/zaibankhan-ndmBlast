# PRD — NDM-BLAST: Bioinformatic Sequence Analysis Web Tool

**Product Requirements Document (PRD)**

| Field | Value |
|-------|-------|
| Product Name | NDM-BLAST |
| Version | 1.0 |
| Status | Draft |
| Author | Mohd Zaiban Khan |
| Institution | M.Sc. Bioinformatics, Jamia Millia Islamia, New Delhi, India |
| Last Updated | 2026-08-30 |

---

## 1. Overview / Problem Statement

Antimicrobial resistance (AMR) is a major global health threat. The **New Delhi
Metallo-β-lactamase (NDM)** family of enzymes is one of the most clinically
significant carbapenemases, rendering the bacteria that produce them resistant
to almost all β-lactam antibiotics.

There is a need for a simple, local, and fast tool that lets researchers and
students identify which NDM variant an unknown **DNA (nucleotide)** or
**protein** sequence belongs to, without relying on remote servers.

**NDM-BLAST** is a lightweight, locally-running web application that accepts a
user-provided sequence, compares it against a curated local database of NDM
variants, and returns the best-matching variant along with similarity scores,
e-values, and detailed variant information.

---

## 2. Goals (Non-Functional & Product Goals)

1. Provide fast, accurate local identification of NDM variants.
2. Offer a clean, intuitive web interface usable without specialized training.
3. Support both **nucleotide (BLASTN)** and **protein (BLASTP)** searches.
4. Require no internet connection after installation (fully local).
5. Display transparent, explainable scoring (similarity %, e-value, alignment).
6. Serve as an educational and research tool for bioinformatics coursework.

### Out of Scope (v1)
- Multiple Sequence Alignment (MSA)
- Download of results (PDF/Excel)
- Large-scale BLAST databases (NCBI remote)
- User registration / authentication
- Cloud / multi-user deployment

---

## 3. Target Users / Personas

| Persona | Needs |
|---------|-------|
| **Bioinformatics Student** | A practical tool to study NDM variant identification; easy-to-read results. |
| **Microbiology Researcher** | Quick screening of unknown NDM sequences from lab data. |
| **Educator** | A classroom demo for sequence homology and BLAST concepts. |

---

## 4. Functional Requirements

### FR-1: Home page
- Display welcome screen and navigation to BLASTN, BLASTP, and Information pages.

### FR-2: BLASTN (Nucleotide Search)
- Accept a DNA sequence (characters **A, T, G, C, N**).
- Validate input:
  - Reject protein-only characters with a clear error message.
  - Enforce a **minimum length of 50 nucleotides**.
- Compare against the local NDM nucleotide database (`NDM_genes.fasta`).
- Return top matching variants sorted by similarity.

### FR-3: BLASTP (Protein Search)
- Accept a protein sequence (standard amino acid characters).
- Validate input:
  - Reject DNA-only sequences with a clear error message.
  - Enforce a **minimum length of 20 amino acids**.
- Compare against the local NDM protein database (`NDM_proteins.fasta`).
- Return top matching variants sorted by similarity.

### FR-4: Scoring
- **Similarity score**: percentage of matching residues aligned against the query.
- **E-value**: statistical significance score computed from similarity, query
  length, and database size.

### FR-5: Results Page
- For each hit, show:
  - Variant name (link to details page)
  - Similarity (%)
  - E-value
  - Sequence description
  - Sequence
- Provide "Back to BLAST results" navigation.

### FR-6: Variant Details Page
- Display metadata for the matched variant:
  - Variant ID
  - Protein / DNA names
  - Country of discovery
  - Host bacteria
  - Detection technique
  - Publication reference(s)
- Show pairwise **alignment** (query vs subject with match line).
- Provide links to:
  - **AlphaFold** structure (when available)
  - **RCSB PDB** structure (when available)
  - **Publication database** (NDM-specific DOI/links)
  - **NCBI** (protein or nucleotide lookup)

### FR-7: Publications Page
- List NDM variant–specific publication links loaded from
  `database/ndm_publications.json`.

### FR-8: Share / Copy Link
- Provide a "Share" button that uses the Web Share API or copies the current
  page URL to the clipboard.

---

## 5. Data Model

### 5.1 Databases
| Database | Format | Purpose |
|----------|--------|---------|
| `NDM_genes.fasta` | FASTA | Nucleotide sequences of NDM variants |
| `NDM_proteins.fasta` | FASTA | Protein sequences of NDM variants |
| `ndm_new.db` | SQLite | Structured metadata for variants |
| `ndm_publications.json` | JSON | Publication / DOI links per variant |

### 5.2 SQLite Tables (from `ndm_new.db`)
- `ndm_nucleotide` — variant_id, description, sequence
- `ndm_protein` — variant_id, description, sequence
- `ndm_variants` — variant metadata (protein name, DNA name, country, host
  bacteria, publication, detection technique, etc.)

### 5.3 Key Constants
- AlphaFold structure links mapped per NDM variant (e.g., NDM-1 → AF-C7C422-F1).
- PDB structure links where available (NDM-1, NDM-4, ...).
- `PROTEIN_ONLY_CHARS` — set used to distinguish protein vs. DNA input.

---

## 6. Technical Architecture

| Layer | Technology | Notes |
|-------|------------|-------|
| Backend | **Python 3** | Business logic, input validation, scoring |
| Web Framework | **Flask** | Routing, templating, session management |
| Frontend | **HTML + CSS + JS** | Server-rendered Jinja2 templates |
| Sequence Parsing | **Biopython** | FASTA parsing (`Bio.SeqIO`) |
| Database | **SQLite** | Variant metadata storage |
| Search Logic | Custom (in `app.py`) | Similarity + e-value computation over local FASTA |

### Pipeline
1. User submits a sequence.
2. Backend detects if the input is DNA or protein.
3. Backend validates length and character set.
4. Sequence is compared to the appropriate local database.
5. Top hits are scored (similarity %) and e-values are computed.
6. Results are rendered; user can drill into variant details.

---

## 7. User Flows

### 7.1 Happy Path — Identify an unknown sequence
1. Home page → click **BLASTN** or **BLASTP**.
2. Paste a sequence into the text area.
3. Submit.
4. View ranked results.
5. Click a variant → see alignment + metadata + external links.

### 7.2 Error Paths
- Input too short → message shown ("Minimum 50 nucleotides" / "20 amino acids").
- Wrong sequence type (protein pasted into BLASTN, or vice-versa) → corrective
  message guiding the user to the correct tool.

---

## 8. Non-Functional Requirements

| Category | Requirement |
|----------|-------------|
| Performance | Local search returns results in < ~1 second for the current database size. |
| Usability | Clean layout, clear error messages in simple language (English / Hinglish), prominent navigation. |
| Portability | Runs on Windows/macOS/Linux with Python ≥ 3.8 and the listed dependencies. |
| Reliability | No external network dependency at runtime; all data bundled locally. |
| Maintainability | FASTA/JSON/SQLite data files are separate from code so updates are easy. |
| Security | Session-based state plus a fixed app secret; no user data stored on disk. |

---

## 9. Dependencies

```text
Flask
Biopython
```

---

## 10. Milestones / Roadmap

| Phase | Scope | Status |
|-------|-------|--------|
| M1 (MVP) | BLASTN + BLASTP + local FASTA search | Done |
| M2 | Variant details page, alignment, metadata | Done |
| M3 | AlphaFold / PDB / NCBI / Publication links | Done |
| M4 | Polished 3D UI, share button, README | Done |
| M5 (Future) | Downstream improvements (see section 11) | Planned |

---

## 11. Future Improvements (Backlog)

- Multiple Sequence Alignment (MSA)
- Download results as PDF / Excel
- Advanced search filters (country, host, year)
- Cloud deployment (e.g., Render / Heroku / Docker)
- Integration with public BLAST / NCBI databases
- User authentication and saved search history

---

## 12. Success Metrics

- Correct identification of NDM variant for known test sequences (accuracy).
- Search completion time (target < 1 s).
- User task-completion for the "identify a sequence" flow.
- Feedback from coursework evaluation / demonstrations.

---

## 13. Appendix

- **Developer**: Mohd Zaiban Khan
- **Contact**: zaibankhan12345@gmail.com
- **Repository**: https://github.com/zaibankhan/zaibankhan-ndmBlast
