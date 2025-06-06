#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script pour fusionner et filtrer les rapports (exon et target) selon le cancer associé à chaque sample.
Si le type de cancer n'est pas précisé dans le fichier de métadonnées, toutes les lignes du rapport sont conservées.
Utilisation :
    python3 filterCancerGenesDNA.py --input_dir <répertoire_reports>
                                    --sample_cancer_file <patient_metadata.csv>
                                    --cancer_gene_file <cancer_gene_dict.tsv>
                                    --output_dir <répertoire_sortie>
"""

import os
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Filtrer et fusionner les rapports exon et target selon le cancer associé à chaque sample"
    )
    parser.add_argument("--input_dir", required=True,
                        help="Répertoire contenant les fichiers exon et target (ex: *.exon_cov_report.bed.tsv et *.target_bed_read_cov_report.bed)")
    parser.add_argument("--sample_cancer_file", required=True,
                        help="Fichier CSV avec le mapping Patient -> Cancer (colonnes 'Patient' et 'Cancer')")
    parser.add_argument("--cancer_gene_file", required=True,
                        help="Fichier TSV avec le mapping Cancer -> Gènes (colonnes 'Cancer' et 'Gene')")
    parser.add_argument("--output_dir", required=True,
                        help="Répertoire de sortie pour les fichiers filtrés")
    args = parser.parse_args()

    input_dir = args.input_dir
    sample_cancer_file = args.sample_cancer_file
    cancer_gene_file = args.cancer_gene_file
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    # 1. Lecture du mapping sample -> cancer
    sample_to_cancer = {}
    try:
        with open(sample_cancer_file, newline='', encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                sample = row['Patient'].strip()
                cancer = row['Cancer'].strip()
                sample_to_cancer[sample] = cancer
    except Exception as e:
        print(f"[Error] Lecture du fichier {sample_cancer_file} impossible: {e}")
        return

    # 2. Lecture du mapping cancer -> gènes
    cancer_to_genes = {}
    try:
        with open(cancer_gene_file, newline='', encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                cancer = row['Cancer'].strip()
                gene = row['Gene'].strip()
                if cancer not in cancer_to_genes:
                    cancer_to_genes[cancer] = set()
                cancer_to_genes[cancer].add(gene)
    except Exception as e:
        print(f"[Error] Lecture du fichier {cancer_gene_file} impossible: {e}")
        return

    # 3. Parcours des fichiers dans input_dir et collecte des lignes par sample
    filtered_data = {}  # sample -> liste de lignes
    exon_names = {}     # sample -> ensemble des valeurs "name" issues des fichiers exon

    for filename in os.listdir(input_dir):
        file_path = os.path.join(input_dir, filename)
        if not os.path.isfile(file_path):
            continue

        # Identification du sample : on suppose que le nom du fichier commence par le nom du sample
        sample_found = None
        for sample in sample_to_cancer.keys():
            if filename.startswith(sample):
                sample_found = sample
                break
        if sample_found is None:
            print(f"[Warning] Aucun sample correspondant trouvé pour le fichier {filename}. Ignoré.")
            continue

        # Récupération du cancer associé et vérification du mapping
        cancer = sample_to_cancer[sample_found]
        if not cancer:
            # Si le type de cancer n'est pas spécifié, on ne filtre pas sur les gènes.
            genes_cancer = None
        elif cancer not in cancer_to_genes:
            print(f"[Warning] Le cancer '{cancer}' pour le sample '{sample_found}' n'est pas présent dans le mapping. Fichier {filename} ignoré.")
            continue
        else:
            genes_cancer = cancer_to_genes[cancer]

        # Détermine le type de fichier (exon ou target)
        is_exon = "exon_cov_report" in filename
        is_target = "target_bed_read_cov_report" in filename

        try:
            with open(file_path, newline='', encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter='\t')
                if 'name' not in reader.fieldnames:
                    print(f"[Warning] La colonne 'name' est absente dans {filename}. Ignoré.")
                    continue
                for row in reader:
                    name_field = row.get("name", "").strip()
                    if not name_field:
                        continue
                    # Extraction du gène (la partie avant le premier underscore)
                    gene = name_field.split('_')[0]
                    # Si genes_cancer est défini, on conserve la ligne uniquement si le gène est dans la liste
                    if genes_cancer is not None:
                        if gene not in genes_cancer:
                            continue
                    # Pour les fichiers exon, on ajoute toujours et on mémorise la valeur "name"
                    if is_exon:
                        filtered_data.setdefault(sample_found, []).append(row)
                        exon_names.setdefault(sample_found, set()).add(name_field)
                    elif is_target:
                        # Pour target, on vérifie que la valeur "name" n'est pas uniquement composée d'éléments contenant "Exon"
                        if '+' in name_field:
                            parts = [part.strip() for part in name_field.split('+')]
                            if all("Exon" in part for part in parts):
                                continue
                        else:
                            if "Exon" in name_field:
                                continue
                        filtered_data.setdefault(sample_found, []).append(row)
                    else:
                        continue
        except Exception as e:
            print(f"[Error] Lecture du fichier {filename} impossible: {e}")

    # 4. Déduplication par sample basée sur la colonne "name"
    for sample, rows in filtered_data.items():
        unique_rows = {}
        for row in rows:
            key = row['name']
            if key not in unique_rows:
                unique_rows[key] = row
        filtered_data[sample] = list(unique_rows.values())

    # 5. Écriture des fichiers filtrés pour chaque sample (colonnes: #chrom, start, end, name)
    for sample, rows in filtered_data.items():
        if not rows:
            continue
        output_file = os.path.join(output_dir, f"{sample}.filtered.tsv")
        fieldnames = ['#chrom', 'start', 'end', 'name']
        try:
            with open(output_file, 'w', newline='', encoding="utf-8") as out_f:
                writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                for row in rows:
                    writer.writerow({k: row[k] for k in fieldnames if k in row})
            print(f"Fichier filtré écrit pour le sample '{sample}': {output_file}")
        except Exception as e:
            print(f"[Error] Impossible d'écrire le fichier pour le sample '{sample}': {e}")

if __name__ == '__main__':
    main()

