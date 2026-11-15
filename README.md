# 🧬 NDM-BLAST

A local web-based bioinformatics application developed using **Python, Flask, HTML, CSS, SQLite, and Biopython** for identifying **NDM (New Delhi Metallo-β-lactamase)** variants using nucleotide and protein sequence analysis.

---

## 📌 About the Project

NDM-BLAST is a bioinformatics web application that allows users to search unknown DNA or protein sequences against a local NDM variant database.

The project provides an easy-to-use interface for sequence analysis and displays detailed information about matched NDM variants.

This project was developed as part of my **M.Sc. Bioinformatics** coursework at **Jamia Millia Islamia, New Delhi**.

---

## ✨ Features

- 🔬 BLASTN (Nucleotide Sequence Search)
- 🧪 BLASTP (Protein Sequence Search)
- 🗂 Local NDM Variant Database
- 📊 Bit Score, E-value, Percentage Identity & Query Coverage (NCBI BLAST methodology)
- 📝 Variant Details Page
- 🌐 User-friendly Web Interface
- ⚡ Fast Local Sequence Search
- 💾 Save Results & History (CSV download or email attachment)

---

## 🛠 Technologies Used

| Technology | Purpose |
|------------|----------|
| Python | Backend Programming |
| Flask | Web Framework |
| HTML | Web Pages |
| CSS | User Interface |
| SQLite | Local Database |
| Biopython | Sequence Analysis |
| SQL | Data Management |

---

## 📂 Project Structure

```text
ndmBlast/
│
├── database/
├── templates/
├── static/
├── app.py
├── insert_all.py
├── requirements.txt
└── README.md
```

---

## 🚀 How to Run

### Clone Repository

```bash
git clone https://github.com/zaibankhan12345/zaibankhan-ndmBlast.git
```

### Move into Project Folder

```bash
cd zaibankhan-ndmBlast
```

### Install Required Packages

```bash
pip install -r requirements.txt
```

### Run the Application

```bash
python app.py
```

### Open Browser

```text
http://127.0.0.1:5000
```

---

## 💾 Save / Export

Every results page and the History page has a **Save** button with two options:

- **Download as CSV** — downloads the data straight to your device.
- **Send to Email** — enter an address and the data is emailed as a CSV attachment.

Results exports include: NDM Variant, Bit Score, E-value, Identity (%), Query Coverage (%), Raw Score, Description, Sequence.
History exports include: Time, Type, Query, Query Length, Top Hit, Identity, E-value.

### Enabling email delivery

Email sending is optional and needs a mail account. Copy `.env.example` to
`.env` and fill in your details (the `.env` file is git-ignored):

```text
NDM_MAIL_HOST=smtp.gmail.com
NDM_MAIL_PORT=587
NDM_MAIL_USE_TLS=1
NDM_MAIL_USER=your.address@gmail.com
NDM_MAIL_PASSWORD=your-16-character-app-password
NDM_MAIL_FROM=your.address@gmail.com
```

> **Gmail:** use an *App Password* (Google Account → Security → App passwords),
> not your normal password. Until an account is configured, the "Send to Email"
> option shows a message explaining what needs to be set.

---

## 📷 Project Screenshots

### Home Page

![Home Page](Users\HP\OneDrive\Desktop\project\homepage.png)

### BLASTN Result

![BLASTN Result](Users\HP\OneDrive\Desktop\project\blastn-result.png)

### BLASTP Result

![BLASTP Result](Users\HP\OneDrive\Desktop\project\blastp-result.png)

---

## 🎯 Future Improvements

- Multiple Sequence Alignment
- Download Results (PDF/Excel)
- Better UI Design
- Advanced Search Filters
- Cloud Deployment
- Public Database Integration
- User Authentication

---

## 👨‍💻 Developer

**Mohd Zaiban Khan**

M.Sc. Bioinformatics  
Jamia Millia Islamia  
New Delhi, India

---

## 📧 Contact

Email: zaibankhan12345@gmail.com

---

## ⭐ If you found this project useful, please give it a Star!
