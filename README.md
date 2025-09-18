# Chip-Seq Analysis Workflow: Kmt5a Controls Hepatic Metabolic Pathways
This repository contains a bioinformatics workflow for a replication of a ChIP-seq analysis focused on the role of Kmt5a in controlling hepatic metabolic pathways. 

`Original paper:` Nikolaou KC, Moulos P, Harokopos V, Chalepakis G, Talianidis I. Kmt5a Controls Hepatic Metabolic Pathways by Facilitating RNA Pol II Release from Promoter-Proximal Regions. Cell Rep. 2017 Jul 25;20(4):909-922. doi: 10.1016/j.celrep.2017.07.003. PMID: 28746875; PMCID: PMC5542855.

`Data:` https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97338

## Repository Structure
The workflow is a step-by-step guide covering the entire process from raw sequencing reads to the final biological interpretation.

`basic_ChIP_seq_analysis.sh`
- Quality Control (FastQC): Initial quality assessment of raw sequencing reads to identify potential issues such as low quality scores, adapter contamination, or other biases.

- Read Trimming and Filtering (cutadapt, minion): Removal of sequencing adapters and low-quality bases from the raw reads to ensure clean data for alignment.

`alignment_input_ctl.sh`
- Alignment (Bowtie2): Mapping the high-quality reads to a reference genome to determine their genomic location.

`peaks.sh`
- Peak Calling (MACS): Identification of genomic regions where the Kmt5a protein is bound. This tool statistically enriches for signal over a control sample to define "peaks."

`final_steps.sh`
- Visualization (IGV): Visualization of the aligned reads and identified peaks on a genome browser to manually inspect the top three peaks for each ChIP sample. This provides a visual confirmation of the peak calling results.

- Motif Analysis (MEME): Discovery of DNA sequence motifs that are enriched within the identified peaks. This helps to determine the consensus binding site for Kmt5a or co-factors.

- Downstream Analysis (SAMtools, etc.): Additional processing and manipulation of the alignment files (BAM/SAM) using SAMtools for tasks such as sorting, indexing, and filtering.

## Prerequisites

- A Linux-based environment

- Software packages: FastQC, cutadapt, Bowtie2, MACS, SAMtools, IGV, and MEME

- Reference genome for alignment (mm9 - NCBI37)

## Analysis Summary
This analysis aims to show that Kmt5a controls hepatic metabolic pathways by facilitating RNA Pol II release from promoter-proximal regions. The ChIP-seq data will provide evidence of Kmt5a binding at specific genomic locations. The downstream analysis will connect these binding events to relevant genes and pathways, suggesting a mechanism by which Kmt5a influences gene expression by affecting the pausing and release of RNA Polymerase II.
