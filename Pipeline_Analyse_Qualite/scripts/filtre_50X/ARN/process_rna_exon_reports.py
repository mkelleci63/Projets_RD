#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script pour traiter les rapports RNA (exon et exon_read) issus du dossier RnaDragenCaller.

Utilisation :
    python3 process_rna_exon_reports.py <chemins_cov_file> <chemins_read_cov_file> <dossier_sortie>

Ce script procède comme suit :
  1. Charge les noms d’exons depuis les fichiers exon_read_cov_report.bed.
  2. Parcourt le fichier listant les chemins des exon_cov_report.bed et, pour chaque fichier,
     filtre les lignes dont la valeur 'pct_above_50' est inférieure à 100 et ajoute le nom d’exon correspondant.
  3. Écrit un fichier de sortie pour chaque rapport traité dans le répertoire spécifié.
Note : Pour RNA, il n’y a pas de fichiers target à copier.
"""

import os
import csv
import argparse
import sys

def load_exon_names(read_cov_paths_file):
    """
    Charge les valeurs de la colonne 'name' des fichiers exon_read_cov_report.bed
    et les stocke dans un dictionnaire {(chrom, start, end): name}.
    """
    exon_name_dict = {}
    
    try:
        with open(read_cov_paths_file, 'r') as f:
            file_paths = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier des chemins exon_read_cov_report.bed: {e}")
        sys.exit(1)

    for file_path in file_paths:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                required_cols = ['#chrom', 'start', 'end', 'name']
                if not reader.fieldnames or not all(col in reader.fieldnames for col in required_cols):
                    print(f"Erreur : Colonnes manquantes dans {file_path}. Fichier ignoré.")
                    continue

                for row in reader:
                    key = (row['#chrom'], row['start'], row['end'])
                    exon_name_dict[key] = row['name']
        except FileNotFoundError:
            print(f"Erreur : fichier {file_path} introuvable.")

    return exon_name_dict

def process_sample_file(input_file, exon_name_dict, output_dir):
    """
    Traite un fichier exon_cov_report.bed :
      - Filtre les lignes où 'pct_above_50' < 100
      - Ajoute le nom d'exon depuis exon_read_cov_report.bed (ou 'N/A' si non trouvé)
      - Écrit le résultat dans le dossier de sortie avec le suffixe .tsv
    """
    try:
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            required_cols = ['#chrom', 'start', 'end', 'pct_above_50']
            if not reader.fieldnames or not all(col in reader.fieldnames for col in required_cols):
                print(f"Attention : Colonnes manquantes dans {input_file}. Fichier ignoré.")
                return
            data = []
            for row in reader:
                try:
                    pct_val = float(row['pct_above_50'])
                except ValueError:
                    print(f"Attention : Valeur invalide pour 'pct_above_50' dans {input_file} : {row['pct_above_50']}. Ligne ignorée.")
                    continue
                if pct_val < 100:
                    key = (row['#chrom'], row['start'], row['end'])
                    row['name'] = exon_name_dict.get(key, 'N/A')
                    # On garde la valeur d'origine pour le tri
                    row['_pct_val'] = pct_val
                    data.append(row)
    except FileNotFoundError:
        print(f"Erreur : le fichier {input_file} n'existe pas.")
        return
    
    data_sorted = sorted(data, key=lambda x: x['_pct_val'])
    
    base_name = os.path.basename(input_file)
    output_file = os.path.join(output_dir, base_name + '.tsv')
    
    try:
        with open(output_file, 'w', newline='') as f_out:
            fieldnames = ['#chrom', 'start', 'end', 'name', 'pct_above_50']
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in data_sorted:
                writer.writerow({k: row[k] for k in fieldnames})
        print(f"Fichier traité : {input_file} -> {output_file}")
    except Exception as e:
        print(f"[Error] Impossible d'écrire le fichier pour {input_file}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Traitement des rapports RNA (exon_cov_report et exon_read_cov_report)"
    )
    parser.add_argument("cov_paths_file", help="Fichier contenant les chemins des exon_cov_report.bed")
    parser.add_argument("read_cov_paths_file", help="Fichier contenant les chemins des exon_read_cov_report.bed")
    parser.add_argument("output_dir", help="Dossier de sortie pour les fichiers traités")
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Charge les noms d'exons depuis les fichiers exon_read_cov_report.bed
    exon_name_dict = load_exon_names(args.read_cov_paths_file)

    try:
        with open(args.cov_paths_file, 'r') as f:
            paths = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier {args.cov_paths_file}: {e}")
        sys.exit(1)
    
    for file_path in paths:
        process_sample_file(file_path, exon_name_dict, args.output_dir)

if __name__ == "__main__":
    main()

