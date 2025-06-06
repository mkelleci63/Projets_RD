#!/usr/bin/env python3
import os
import sys
import glob
import shutil
import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Copie les fichiers MetricsOutput.tsv depuis l'arborescence et consolide les métriques."
    )
    parser.add_argument("--parent_directory", type=str, required=True,
                        help="Dossier parent à explorer pour trouver les fichiers MetricsOutput.tsv")
    parser.add_argument("--output_file", type=str, required=True,
                        help="Chemin complet du fichier de sortie MetricsSummary_Global.tsv")
    parser.add_argument("--data_dir", type=str, required=True,
                        help="Dossier de travail pour copier les fichiers MetricsOutput.tsv")
    return parser.parse_args()

def copy_metrics_files(parent_path, data_dir):
    """
    Parcours le dossier parent pour trouver, dans les sous-dossiers (nom à 4 caractères commençant par '2')
    et leurs sous-dossiers contenant 'TSO500', les dossiers 'Results' et copie les fichiers MetricsOutput.tsv
    dans data_dir.
    """
    os.makedirs(data_dir, exist_ok=True)
    print(f"Création (si nécessaire) du dossier de destination : {data_dir}")
    
    # Itération sur les sous-dossiers du dossier parent (par exemple des années)
    for entry in os.listdir(parent_path):
        entry_path = os.path.join(parent_path, entry)
        if os.path.isdir(entry_path) and len(entry) == 4 and entry.startswith("2"):
            print(f"Exploration de l'année : {entry}")
            # Dans le dossier d'année, on cible uniquement les dossiers contenant "TSO500"
            for sub_entry in os.listdir(entry_path):
                sub_entry_path = os.path.join(entry_path, sub_entry)
                if os.path.isdir(sub_entry_path) and "TSO500" in sub_entry:
                    print(f"  Dossier TSO500 trouvé : {sub_entry}")
                    # Recherche récursive des dossiers "Results" dans ce dossier
                    for root, dirs, files in os.walk(sub_entry_path):
                        if os.path.basename(root) == "Results":
                            metrics_file = os.path.join(root, "MetricsOutput.tsv")
                            if os.path.isfile(metrics_file):
                                # Création d'un nom de fichier unique basé sur le chemin relatif d'origine
                                relative_path = os.path.relpath(root, parent_path)
                                safe_name = relative_path.replace(os.sep, "_")
                                destination_file = os.path.join(data_dir, f"{safe_name}_MetricsOutput.tsv")
                                shutil.copy2(metrics_file, destination_file)
                                print(f"    Copié : {metrics_file} -> {destination_file}")
                            else:
                                print(f"    Aucun fichier MetricsOutput.tsv dans {root}")
    return data_dir

def safe_float(value):
    if value == "NA":
        return np.nan
    try:
        return float(value)
    except ValueError:
        return np.nan

def process_section(section_data):
    """
    Traite une section d'un fichier MetricsOutput.tsv et retourne un dictionnaire
    des métriques avec leurs valeurs, échantillons et seuils.
    """
    metrics = {}
    header = section_data[0].split("\t")
    samples = header[3:]
    for row in section_data[1:]:
        row_data = row.split("\t")
        if len(row_data) < 4:
            continue
        metric_name = row_data[0]
        lsl = safe_float(row_data[1])
        usl = safe_float(row_data[2])
        values = [safe_float(x) for x in row_data[3:]]
        metrics[metric_name] = {
            "Values": values,
            "Samples": samples,
            "LSL": lsl,
            "USL": usl
        }
    return metrics

def consolidate_metrics_from_files(data_dir, output_file_path):
    """
    Parcours les fichiers MetricsOutput.tsv présents dans data_dir, consolide les métriques
    et génère un fichier TSV de sortie.
    """
    file_paths = glob.glob(os.path.join(data_dir, "*MetricsOutput.tsv"))
    if not file_paths:
        print("Aucun fichier MetricsOutput.tsv trouvé dans", data_dir)
        sys.exit(1)
    
    global_metrics = {}
    
    for file_path in file_paths:
        print(f"Traitement du fichier : {file_path}")
        with open(file_path, "r") as file:
            lines = file.readlines()
        
        sections = {}
        current_section = None
        for line in lines:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                current_section = line.strip("[]")
                sections[current_section] = []
            elif current_section and line:
                sections[current_section].append(line)
        
        for section, section_data in sections.items():
            if section_data:
                metrics = process_section(section_data)
                for metric, data in metrics.items():
                    if metric not in global_metrics:
                        global_metrics[metric] = {"Values": [], "Samples": [], "LSL": data["LSL"], "USL": data["USL"]}
                    global_metrics[metric]["Values"].extend(data["Values"])
                    global_metrics[metric]["Samples"].extend(data["Samples"])
    
    output_data = []
    for metric, data in global_metrics.items():
        values = data["Values"]
        samples = data["Samples"]
        lsl = data["LSL"]
        usl = data["USL"]
        valid_values = [v for v in values if not np.isnan(v)]
        if valid_values:
            mean = np.mean(valid_values)
            std_dev = np.std(valid_values)
            p5 = np.percentile(valid_values, 5)
            p95 = np.percentile(valid_values, 95)
            seuil_inf = lsl if not np.isnan(lsl) else p5
            seuil_sup = usl if not np.isnan(usl) else p95
            used_samples = sorted([samples[i] for i in range(len(values)) if not np.isnan(values[i])])
            output_data.append([
                metric, mean, std_dev, seuil_inf, seuil_sup, ", ".join(used_samples)
            ])
    
    output_df = pd.DataFrame(output_data, columns=[
        "Metric (UOM)", "Moyenne", "Ecart-type", "Seuil Inf", "Seuil Sup", "Echantillons utilisés"
    ])
    
    dir_output = os.path.dirname(output_file_path)
    if dir_output:
        os.makedirs(dir_output, exist_ok=True)
    
    output_df.to_csv(output_file_path, sep="\t", index=False)
    print("Fichier de sortie généré :", output_file_path)

def main():
    args = parse_arguments()
    
    print("=== Phase de copie des fichiers MetricsOutput.tsv ===")
    copied_dir = copy_metrics_files(args.parent_directory, args.data_dir)
    
    print("\n=== Phase de consolidation des métriques ===")
    consolidate_metrics_from_files(copied_dir, args.output_file)

if __name__ == "__main__":
    main()

