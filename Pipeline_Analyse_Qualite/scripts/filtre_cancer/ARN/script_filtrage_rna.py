#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script pour filtrer les rapports RNA (exon_cov_report) selon le cancer associé à chaque sample,
à partir d'un JSON d'alias de gènes incluant les synonymes et leurs cancers.
La colonne 'name' est remplacée par le gene_id canonique pour toutes les lignes.
Utilisation :
    python3 filterCancerGenesRNA.py --input_dir <reports_RNA_dir> \
                                    --sample_cancer_file <patient_metadata.csv> \
                                    --alias_cancer_file <alias_cancers.json> \
                                    --output_dir <output_dir>
"""

import os
import csv
import json
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Filtrer les rapports RNA selon le cancer associé à chaque sample"
    )
    parser.add_argument("--input_dir", required=True,
                        help="Répertoire contenant les rapports RNA (.exon_cov_report.bed.tsv)")
    parser.add_argument("--sample_cancer_file", required=True,
                        help="Fichier CSV avec le mapping Patient -> Cancer (colonnes 'Patient' et 'Cancer')")
    parser.add_argument("--alias_cancer_file", required=True,
                        help="Fichier JSON mapping alias de gène -> {gene_id, cancers}")
    parser.add_argument("--output_dir", required=True,
                        help="Répertoire de sortie pour les fichiers filtrés")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Lecture mapping sample -> cancer
    sample_to_cancer = {}
    with open(args.sample_cancer_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            sample = row.get('Patient', '').strip()
            cancer = row.get('Cancer', '').strip()
            sample_to_cancer[sample] = cancer

    # 2. Lecture du JSON alias -> {gene_id, cancers}
    with open(args.alias_cancer_file, encoding='utf-8') as f:
        alias_map = json.load(f)
        # alias_map[alias] = {'gene_id': ..., 'cancers': [...]}  

    # 3. Parcours des fichiers dans input_dir
    for filename in os.listdir(args.input_dir):
        file_path = os.path.join(args.input_dir, filename)
        if not os.path.isfile(file_path):
            continue

        # Extraction du sample depuis le nom de fichier
        original_sample = filename.split('.')[0]
        lookup_sample = original_sample.replace('R-', 'D-')
        cancer = sample_to_cancer.get(lookup_sample)
        has_cancer = bool(cancer)

        filtered_rows = []
        try:
            with open(file_path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                if 'name' not in reader.fieldnames:
                    continue
                for row in reader:
                    name_field = row.get('name', '').strip()
                    if not name_field:
                        continue

                    # Extraction de l'alias (avant le premier underscore)
                    alias, sep, rest = name_field.partition('_')

                    # Vérification du cancer
                    if has_cancer:
                        entry = alias_map.get(alias)
                        if not entry or cancer not in entry.get('cancers', []):
                            continue
                    # Substitution du champ 'name' par le gene_id canonique
                    entry = alias_map.get(alias)
                    canonical = entry['gene_id'] if entry else alias
                    row['name'] = canonical + (sep + rest if sep else '')

                    # Conserver la ligne
                    filtered_rows.append(row)
        except Exception as e:
            print(f"[Error] Impossible de lire {filename}: {e}")
            continue

        if not filtered_rows:
            print(f"[Info] Aucun résultat filtré pour le sample '{original_sample}'.")
            continue

        # 4. Écriture du fichier filtré
        output_file = os.path.join(args.output_dir, f"{original_sample}.filtered.tsv")
        fieldnames = ['#chrom', 'start', 'end', 'name']
        with open(output_file, 'w', newline='', encoding='utf-8') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for r in filtered_rows:
                writer.writerow({k: r[k] for k in fieldnames if k in r})

        print(f"Écrit: {output_file} ({len(filtered_rows)} lignes)")

if __name__ == '__main__':
    main()

