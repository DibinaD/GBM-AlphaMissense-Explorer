Glioblastoma Multiforme (GBM) is one of the most aggressive brain cancers with a complex mutational landscape. Identifying functionally important mutations helps guide research and therapeutic development. This project integrates:
TCGA GBM mutation data (PanCancer Atlas)
AlphaMissense AI scores for all possible missense variants (DeepMind)
ClinVar clinical significance annotations
COSMIC somatic mutation catalog and presents the results in a user-friendly R Shiny interface.
-----------------------------------------------------------------------------------------------
Key Features
-------------
AI-driven annotation
Leverage AlphaMissense pathogenicity scores (0–1) and class labels (benign/pathogenic).

Database validation filters
ClinVar-confirmed pathogenic variants
COSMIC-listed cancer mutations

Interactive exploration
Multi-gene selection
Score sliders and checkbox filters
Real-time table and plot updates

Downloadable results
Export filtered mutation sets as CSV

Clean, focused data views
Essential columns (gene, position, protein change, AI score, clinical tags)

GBM-AlphaMissense-Explorer/
├── app.R                            # R Shiny app
├── data_cleaning_and_annotation.py  # Python script for AlphaMissense merge
├── data/                            # Example directory (ignored by Git)
│   ├── AlphaMissense_hg19.tsv       # Raw AI prediction (large file)
│   └── data_mutations.txt           # TCGA MAF (large file)
├── annotated_with_alphamissense.csv # Processed annotation
├── .gitignore                       # large data and temp files
└── README.md                        # This file

Data Download Instructions
--------------------------
Note: Raw and processed data files are large (>50 MB) and not included in this repo.

AlphaMissense predictions (hg19)
Download from Zenodo: https://zenodo.org/records/8360242/files/AlphaMissense_hg19.tsv.gz

TCGA GBM MAF file
Obtain from GDC Data Portal or cBioPortal (PanCancer Atlas GBM cohort)

Processed annotation CSV
Run the Python script data_cleaning_and_annotation.py to merge AlphaMissense with TCGA data and produce annotated_with_alphamissense.csv.

Place the raw files in data/ and the processed CSV at the repo root before running the Shiny app.

Usage
-----------
Select Gene(s): Choose one or more genes to focus on.
Score Filter: Adjust the pathogenicity score slider to narrow AI predictions.
Pathogenic Only: Show only variants classified as pathogenic by AlphaMissense.
ClinVar Only: Show clinically validated mutations.
COSMIC Only: Show variants listed in the COSMIC database.
Download: Export your current filter results as CSV.
Explore tables and interactive plots under different tabs to gain insights.




