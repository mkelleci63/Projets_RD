#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script pour filtrer les rapports RNA (exon_cov_report) selon le cancer associé à chaque sample.
Si le type de cancer n'est pas spécifié dans le fichier de métadonnées, toutes les lignes du rapport sont conservées.
Utilisation :
    python3 filterCancerGenesRNA.py --input_dir <répertoire_reports_RNA>
                                    --sample_cancer_file <patient_metadata.csv>
                                    --cancer_gene_file <cancer_gene_dict.tsv>
                                    --output_dir <répertoire_sortie>
"""

import os
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Filtrer les rapports RNA selon le cancer associé au sample"
    )
    parser.add_argument("--input_dir", required=True, help="Répertoire contenant les rapports RNA")
    parser.add_argument("--sample_cancer_file", required=True, help="Fichier CSV avec le mapping Patient -> Cancer")
    parser.add_argument("--cancer_gene_file", required=True, help="Fichier TSV avec le mapping Cancer -> Gènes")
    parser.add_argument("--output_dir", required=True, help="Répertoire de sortie pour les fichiers filtrés")
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)

    # Lecture mapping sample -> cancer
    sample_to_cancer = {}
    try:
        with open(args.sample_cancer_file, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                sample = row['Patient'].strip()
                cancer = row['Cancer'].strip()
                sample_to_cancer[sample] = cancer
    except Exception as e:
        print(f"[Error] Impossible de lire {args.sample_cancer_file}: {e}")
        return

    # Lecture mapping cancer -> gènes
    cancer_to_genes = {}
    try:
        with open(args.cancer_gene_file, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                cancer = row['Cancer'].strip()
                gene = row['Gene'].strip()
                if cancer not in cancer_to_genes:
                    cancer_to_genes[cancer] = set()
                cancer_to_genes[cancer].add(gene)
    except Exception as e:
        print(f"[Error] Impossible de lire {args.cancer_gene_file}: {e}")
        return

    # Parcours des fichiers dans le répertoire d'input
    for filename in os.listdir(args.input_dir):
        file_path = os.path.join(args.input_dir, filename)
        if not os.path.isfile(file_path):
            continue
        # Extraction du nom original du sample depuis le nom de fichier
        original_sample = filename.split('.')[0]
        # Pour la recherche dans le mapping, on remplace 'R-' par 'D-' afin d'avoir la correspondance avec le metadata
        lookup_sample = original_sample.replace('R-', 'D-')
        
        if lookup_sample not in sample_to_cancer:
            print(f"[Warning] Sample '{lookup_sample}' non trouvé dans le mapping. Fichier {filename} ignoré.")
            continue
        
        cancer = sample_to_cancer[lookup_sample]
        if not cancer:
            # Si le type de cancer n'est pas spécifié, on ne filtre pas sur les gènes.
            genes_cancer = None
        elif cancer not in cancer_to_genes:
            print(f"[Warning] Cancer '{cancer}' pour le sample '{lookup_sample}' non trouvé dans le mapping. Fichier {filename} ignoré.")
            continue
        else:
            genes_cancer = cancer_to_genes[cancer]
        
        filtered_rows = []
        try:
            with open(file_path, newline='', encoding='utf-8') as f:
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
                    if genes_cancer is not None:
                        if gene not in genes_cancer:
                            continue
                    filtered_rows.append(row)
        except Exception as e:
            print(f"[Error] Impossible de lire {filename}: {e}")
            continue

        if not filtered_rows:
            print(f"[Info] Aucun résultat filtré pour le sample '{lookup_sample}'.")
            continue

        output_file = os.path.join(args.output_dir, f"{original_sample}.filtered.tsv")
        fieldnames = filtered_rows[0].keys()
        try:
            with open(output_file, 'w', newline='', encoding='utf-8') as out_f:
                writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                for row in filtered_rows:
                    writer.writerow(row)
            print(f"Fichier filtré écrit pour le sample '{original_sample}': {output_file}")
        except Exception as e:
            print(f"[Error] Impossible d'écrire le fichier pour le sample '{original_sample}': {e}")

if __name__ == '__main__':
    main()

