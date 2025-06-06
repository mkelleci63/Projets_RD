#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script pour fusionner et filtrer les rapports (exon et target) selon le cancer associé à chaque sample,
à partir d'un JSON d'alias de gènes incluant les synonymes et leurs cancers.
La colonne 'name' est remplacée par le gene_id canonique pour toutes les lignes.
Utilisation :
    python3 filterCancerGenesDNA.py --input_dir <reports_dir> \
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
        description="Filtrer et fusionner les rapports exon et target selon le cancer associé à chaque sample"
    )
    parser.add_argument("--input_dir", required=True,
                        help="Répertoire contenant les fichiers exon et target (.exon_cov_report.bed.tsv et .target_bed_read_cov_report.bed)")
    parser.add_argument("--sample_cancer_file", required=True,
                        help="Fichier CSV avec le mapping Patient -> Cancer (colonnes 'Patient' et 'Cancer')")
    parser.add_argument("--alias_cancer_file", required=True,
                        help="Fichier JSON mapping alias de gène -> {gene_id, cancers}")
    parser.add_argument("--output_dir", required=True,
                        help="Répertoire de sortie pour les fichiers filtrés")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Lecture du mapping sample -> cancer
    sample_to_cancer = {}
    with open(args.sample_cancer_file, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            sample = row.get('Patient', '').strip()
            cancer = row.get('Cancer', '').strip()
            sample_to_cancer[sample] = cancer

    # 2. Lecture du JSON alias -> {gene_id, cancers}
    with open(args.alias_cancer_file, encoding="utf-8") as f:
        alias_map = json.load(f)

    filtered_data = {}  # sample -> list of rows

    # 3. Parcours des fichiers de rapports
    for filename in os.listdir(args.input_dir):
        file_path = os.path.join(args.input_dir, filename)
        if not os.path.isfile(file_path):
            continue

        # Détermination du sample par préfixe de nom de fichier
        sample = None
        for s in sample_to_cancer:
            if filename.startswith(s):
                sample = s
                break
        if sample is None:
            continue

        cancer = sample_to_cancer.get(sample)
        has_cancer = bool(cancer)

        is_exon = 'exon_cov_report' in filename
        is_target = 'target_bed_read_cov_report' in filename

        with open(file_path, newline='', encoding="utf-8") as f:
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
                # Reconstruire name avec suffixe d'origine
                row['name'] = canonical + (sep + rest if sep else '')

                # Filtrage spécifique aux fichiers 'target'
                if is_target:
                    parts = [p.strip() for p in name_field.split('+')]
                    if all('Exon' in p for p in parts):
                        continue
                    if len(parts) == 1 and 'Exon' in parts[0]:
                        continue

                # Conserver la ligne
                filtered_data.setdefault(sample, []).append(row)

    # 4. Déduplication par sample sur la colonne 'name'
    for sample, rows in filtered_data.items():
        unique = {}
        for row in rows:
            key = row['name']
            if key not in unique:
                unique[key] = row
        filtered_data[sample] = list(unique.values())

    # 5. Écriture des fichiers filtrés
    for sample, rows in filtered_data.items():
        if not rows:
            continue
        out_file = os.path.join(args.output_dir, f"{sample}.filtered.tsv")
        fieldnames = ['#chrom','start','end','name']
        with open(out_file, 'w', newline='', encoding="utf-8") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r[k] for k in fieldnames if k in r})

        print(f"Écrit: {out_file} ({len(rows)} lignes)")

if __name__ == '__main__':
    main()

