# Bioinformatics AI/ML Research 🚀

This repository contains my work from my **AI/ML Bioinformatics Internship**. 

- Applied **Machine Learning & Deep Learning** techniques to **bioinformatics research problems**.  
- Developed **CNN models** for **biomedical image classification** on molecular datasets.  
- Explored **AI applications in gene sequencing** and **molecular data analysis**. 

---
## 🧪 Research Use Case

This project demonstrates how **AI/ML can accelerate bioinformatics research**:

1. **Protein Sequence Alignment**  
   - Compare evolutionary similarity between species.  
   - Identify conserved protein regions relevant for **drug target discovery**.  

2. **Machine Learning Models**  
   - **Random Forest & SVM**: Useful for **gene expression analysis** and **biomarker identification**.  
   - **CNN**: Prototype for **biomedical image classification**, which can be extended to protein structure images or histopathology slides.  

📌 These workflows can support **personalized medicine, molecular diagnostics, and computational drug discovery** research.  


## 🔬 Contents

### 1. Sequence Alignment (`sequence_alignment/`)
Protein sequence alignment experiments using **Biopython**.  
- Global and local pairwise alignments  
- Substitution matrices: **BLOSUM62**, **PAM250**  
- Similarity score and percentage identity  

| Script         | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `decorin.py`   | Basic alignment with **BLOSUM62**; calculates matches and percentage similarity |
| `laminin.py`   | Global alignment with **BLOSUM62** and **PAM250** for comparison            |
| `gastrin.py`   | Comprehensive: global + local alignments, tested on both matrices           |

---

### 2. Molecular Models (`molecular_models/`)
Machine Learning and Deep Learning models applied to synthetic/public datasets to **demonstrate workflows**.  

- `cnn_molecular_classification.py` → Convolutional Neural Network (demo on CIFAR-10, adapted to **biomedical image classification**)  
- `ml_models.py` → Multiple ML algorithms on synthetic molecular datasets:
  - **Random Forest**
  - **Support Vector Machine (SVM)**
  - **Logistic Regression**

---

## 📂 Repository Structure

```
bioinformatics-ml-research/
│
├── data/                          # (FASTA files, synthetic/public datasets)
│   └── sample_sequences/
│       ├── seq1.fasta
│       ├── seq2.fasta
│       ├── seq3.fasta
│       └── seq4.fasta
│
├── molecular_models/              # ML + CNN prototypes
│   ├── cnn_molecular_classification.py
│   └── ml_models.py
│
├── sequence_alignment/            # Sequence alignment prototypes
│   ├── decorin.py
│   ├── laminin.py
│   └── gastrin.py
│
├── Introduction to BioInformatics.pdf   # Internship reference/notes
├── README.md                      # Documentation
├── requirements.txt                # Dependencies
```

---

## ⚠️ Note
Due to confidentiality, **real datasets are not shared**.  
The provided code demonstrates the **workflow and methodology** using public or synthetic datasets.

---

## 🚀 Example Results
- **CNN**: ~92% accuracy (demo on CIFAR-10 dataset)  
- **Random Forest**: ~85% accuracy (synthetic molecular features)  
- **SVM**: ~83% accuracy (synthetic gene expression dataset)  

---

## 👨‍💻 Author
**Farhan Naushad**  
🔗 [LinkedIn](https://www.linkedin.com/in/farhannaushad01)  
🔗 [GitHub](https://github.com/farhannaushad08)
