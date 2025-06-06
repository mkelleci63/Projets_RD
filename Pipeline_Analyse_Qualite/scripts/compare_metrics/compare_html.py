#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import json

# Définition et parsing des arguments
parser = argparse.ArgumentParser(
    description="Générer un rapport HTML à partir de MetricsOutput et MetricsSummary_Global."
)
parser.add_argument("--global_file", required=True, help="Chemin vers le fichier MetricsSummary_Global.tsv")
parser.add_argument("--input_file", required=True, help="Chemin vers le fichier MetricsOutput.tsv")
parser.add_argument("--output_file", required=True, help="Chemin vers le fichier HTML généré")
args = parser.parse_args()

# Charge le fichier des seuils globaux
seuils_df = pd.read_csv(args.global_file, sep="\t")

# Charge le fichier MetricsOutput
with open(args.input_file, "r") as file:
    lines = file.readlines()

# --- Extraction des données et seuils locaux ---
def process_section(section_data):
    metrics = {}
    header = section_data[0].split("\t")
    # Colonnes : 0 = nom de la métrique, 1 = LSL Guideline, 2 = USL Guideline, 3+ = valeurs par échantillon
    samples = header[3:]
    for row in section_data[1:]:
        row_data = row.split("\t")
        metric_name = row_data[0]
        try:
            lsl = float(row_data[1]) if row_data[1] not in ["NA", ""] else None
        except:
            lsl = None
        try:
            usl = float(row_data[2]) if row_data[2] not in ["NA", ""] else None
        except:
            usl = None
        values = [float(x) if x != "NA" else np.nan for x in row_data[3:]]
        metrics[metric_name] = {"LSL": lsl, "USL": usl, "Values": dict(zip(samples, values))}
    return metrics

# Extraction des sections
sections = {}
current_section = None
for line in lines:
    line = line.strip()
    if line.startswith("["):
        current_section = line.strip("[]")
        sections[current_section] = []
    elif current_section and line:
        sections[current_section].append(line)

# Sections pertinentes
relevant_sections = [
    "Run QC Metrics",
    "DNA Library QC Metrics",
    "DNA Library QC Metrics for Small Variant Calling and TMB",
    "DNA Library QC Metrics for MSI",
    "DNA Library QC Metrics for CNV",
    "RNA Library QC Metrics",
    "DNA Expanded Metrics",
    "RNA Expanded Metrics"
]

# Récupération de toutes les métriques et tous les échantillons
all_metrics = set(seuils_df["Metric (UOM)"].values)
all_samples = set()
for section in relevant_sections:
    if section in sections:
        metrics = process_section(sections[section])
        for metric, data in metrics.items():
            all_samples.update(data["Values"].keys())
all_samples = sorted(all_samples)

# Couleurs par section
section_colors = {
    "Run QC Metrics": "#FFCCCC",
    "DNA Library QC Metrics": "#CCE5FF",
    "DNA Library QC Metrics for Small Variant Calling and TMB": "#CCFFCC",
    "DNA Library QC Metrics for MSI": "#FFFFCC",
    "DNA Library QC Metrics for CNV": "#FFCCFF",
    "RNA Library QC Metrics": "#FFE5CC",
    "DNA Expanded Metrics": "#CCFFFF",
    "RNA Expanded Metrics": "#E5CCFF"
}

# Dictionnaire d'explications pour les métriques
metric_explanations = {
    "PCT_Q30_R1 (%)": "Percentage of Read 1 reads with quality score ≥ 30",
    "PCT_PF_READS (%)": "Total percentage of reads passing filter",
    "PCT_Q30_R2 (%)": "Percentage of Read 2 reads with quality score ≥ 30",
    "CONTAMINATION_SCORE (NA)": "The contamination score is based on VAF distribution of SNPs",
    "PCT_EXON_50X (%)": "Percent exon bases with 50x fragment coverage",
    "MEDIAN_INSERT_SIZE (bp)": "The median fragment length in the sample",
    "MEDIAN_EXON_COVERAGE (Count)": "Median exon fragment coverage across all exon bases",
    "USABLE_MSI_SITES (Count)": "The number of MSI sites usable for MSI calling",
    "MEDIAN_BIN_COUNT_CNV_TARGET (Count)": "The median raw bin count per CNV target",
    "GENE_SCALED_MAD (Count)": "GENE_SCALED_MAD (Count)",
    "TOTAL_ON_TARGET_READS (Count)": "The total number of reads that map to the target regions",
    "MEDIAN_CV_GENE_500X (NA)": ("The median CV for all genes with median coverage > 500x. "
                                  "Higher CV indicates possible issues lors de la préparation de la librairie."),
    "MEDIAN_INSERT_SIZE (Count)": "The median fragment length in the sample",
    "PCT_TARGET_250X (%)": "Percentage of target bases with 250X fragment coverage.",
    "MEDIAN_TARGET_COVERAGE (Count)": "Median depth across all target loci",
    "PCT_SOFT_CLIPPED_BASES (%)": "Percentage of bases soft-clipped during alignement",
    "PCT_READ_ENRICHMENT (%)": "Percentage of reads overlapping les régions cibles",
    "MEDIAN_TARGET_HRD_COVERAGE (Count)": "Median target coverage for HRD analysis",
    "PCT_ALIGNED_READS (%)": "Proportion of aligned reads that pass QC",
    "PCT_CHIMERIC_READS (%)": "Percentage of reads alignées en deux parties",
    "PCT_EXON_100X (%)": "Percentage of exon bases with 100X coverage",
    "MEAN_TARGET_COVERAGE (Count)": "Mean depth across all target loci",
    "PCT_TARGET_0.4X_MEAN (%)": "Percentage of target bases with coverage ≥ 0.4x the mean",
    "TOTAL_PF_READS (Count)": "Total number of reads passing filter",
    "MEAN_FAMILY_SIZE (Count)": "Mean size of UMI families",
    "PCT_TARGET_50X (%)": "Percentage of target bases with 50X coverage",
    "PCT_TARGET_100X (%)": "Percentage of target bases with 100X coverage",
    "PCT_USABLE_UMI_READS (%)": "Percentage of reads with valid UMI",
    "PCT_Q30_BASES (%)": "Average percentage of bases ≥ Q30",
    "PCT_CONTAMINATION_EST (%)": "Estimated contamination percentage",
    "ALLELE_DOSAGE_RATIO (NA)": "Allele Dosage Ratio estimate",
    "PCT_ON_TARGET_READS (%)": "Percentage of reads on target",
    "SCALED_MEDIAN_GENE_COVERAGE (Count)": "Scaled median coverage par gène",
    "RNA_PCT_Q30_BASES (%)": "Average percentage of RNA bases ≥ Q30",
    "GENE_ABOVE_MEDIAN_CUTOFF (Count)": "Number of genes above median coverage cutoff",
    "GENE_MEDIAN_COVERAGE (Count)": "Median coverage depth per gene"
}

# --- Fonction d'évaluation générale (pour toutes les métriques sauf MEAN_FAMILY_SIZE) ---
def evaluate_metric(value, local_lsl, local_usl, global_lsl, global_usl):
    """
    Renvoie un tuple (bg_color, tooltip, add_star)
    Si la valeur est inférieure ou supérieure aux seuils attendus, le tooltip inclut le seuil exact.
    """
    # Priorité aux seuils LSL/USL d'ILLUMINA s'ils sont définis
    if local_lsl is not None or local_usl is not None:
        if local_lsl is not None and local_usl is not None:
            if value < local_lsl:
                return ("red", f"Valeur non conforme (inférieure au seuil de {local_lsl:.2f})", False)
            elif value > local_usl:
                return ("red", f"Valeur non conforme (supérieure au seuil de {local_usl:.2f})", False)
            else:
                return ("white", "Valeur conforme", False)
        elif local_lsl is None and local_usl is not None:
            if value > local_usl:
                return ("red", f"Valeur non conforme (supérieure au seuil de {local_usl:.2f})", False)
            else:
                return ("white", "Valeur conforme", False)
        elif local_usl is None and local_lsl is not None:
            if value < local_lsl:
                return ("red", f"Valeur non conforme (inférieure au seuil de {local_lsl:.2f})", False)
            else:
                return ("white", "Valeur conforme", False)
    else:
        if global_lsl is not None and value < global_lsl:
            return ("red", f"Valeur non conforme (inférieure au seuil de {global_lsl:.2f})", False)
        elif global_usl is not None and value > global_usl:
            return ("white", f"Valeur au-dessus du seuil moyen supérieur du groupe ({global_usl:.2f})", True)
        else:
            return ("white", "Valeur conforme", False)

# Génére la ligne d'en-tête dans une balise <thead>
header_row = "<tr><th>Metric</th>" + "".join([f"<th>{sample}</th>" for sample in all_samples]) + "</tr>"

# Génére le corps du tableau dans une balise <tbody>
body_rows = ""
for section in relevant_sections:
    if section in sections:
        metrics = process_section(sections[section])
        for metric in all_metrics:
            if metric in metrics:
                local_lsl = metrics[metric].get("LSL")
                local_usl = metrics[metric].get("USL")
                global_vals = seuils_df.loc[seuils_df["Metric (UOM)"] == metric, ["Seuil Inf", "Seuil Sup"]].values
                global_lsl = global_vals[0][0] if len(global_vals) > 0 and not pd.isna(global_vals[0][0]) else None
                global_usl = global_vals[0][1] if len(global_vals) > 0 and not pd.isna(global_vals[0][1]) else None

                # Pour la métrique MEAN_FAMILY_SIZE (Count), seuils à 2 et 5
                if metric.upper() == "MEAN_FAMILY_SIZE (COUNT)":
                    seuil_text = "<small>LSL: 2.00, USL: 5.00</small>"
                else:
                    # Formate les seuils avec 2 décimales
                    if local_lsl is not None or local_usl is not None:
                        lsl_str = f"{local_lsl:.2f}" if local_lsl is not None else "N/A"
                        usl_str = f"{local_usl:.2f}" if local_usl is not None else "N/A"
                        seuil_text = f"<small>LSL: {lsl_str}, USL: {usl_str}</small>"
                    else:
                        inf_str = f"{global_lsl:.2f}" if global_lsl is not None else "N/A"
                        sup_str = f"{global_usl:.2f}" if global_usl is not None else "N/A"
                        seuil_text = f"<small>Global: Inf: {inf_str}, Sup: {sup_str}</small>"

                values_dict = metrics[metric]["Values"]
                color_section = section_colors.get(section, "#FFFFFF")
                tooltip_base = metric_explanations.get(metric, "")
                row_html = f"<tr><td class='metric-name' title='{tooltip_base}' style='background-color:{color_section}'>{metric}{seuil_text}</td>"
                for sample in all_samples:
                    value = values_dict.get(sample, np.nan)
                    if not np.isnan(value):
                        # Cas particulier pour MEAN_FAMILY_SIZE (Count)
                        if metric.upper() == "MEAN_FAMILY_SIZE (COUNT)":
                            if value < 2:
                                cell_color = "red"
                                cell_tooltip = "Valeur non conforme (inférieure au seuil de 2.00)"
                                add_star = False
                            elif 2 <= value <= 5:
                                cell_color = "white"
                                cell_tooltip = "Valeur conforme"
                                add_star = False
                            elif 5 < value <= 10:
                                cell_color = "orange"
                                cell_tooltip = "Valeur borderline (entre 5.00 et 10.00)"
                                add_star = False
                            elif value > 10:
                                cell_color = "red"
                                cell_tooltip = "Valeur non conforme (supérieure au seuil de 10.00)"
                                add_star = False
                        else:
                            cell_color, cell_tooltip, add_star = evaluate_metric(value, local_lsl, local_usl, global_lsl, global_usl)
                        
                        if add_star:
                            row_html += f"<td title='{cell_tooltip}' style='position: relative;'>{value}<span class='green-circle'></span></td>"
                        else:
                            row_html += f"<td title='{cell_tooltip}' style='background-color:{cell_color};'>{value}</td>"
                    else:
                        row_html += f"<td title='Pas de donnée'>NA</td>"
                row_html += "</tr>"
                body_rows += row_html

# Construction du HTML final avec header cloné en haut et container d'information fixé en bas
html_content = f"""
<html>
<head>
    <meta charset="UTF-8">
    <title>Metrics Comparison</title>
    <style>
        html, body {{
            margin: 0;
            padding: 0;
            height: 100%;
        }}
        body {{
            font-family: Arial, sans-serif;
        }}
        .legend {{
            padding: 10px;
            border-bottom: 2px solid black;
        }}
        .legend h3 {{
            margin: 0;
        }}
        .legend-item {{
            display: inline-block;
            margin-right: 15px;
        }}
        .legend-color {{
            width: 20px;
            height: 20px;
            display: inline-block;
            margin-right: 5px;
            border: 1px solid black;
        }}
        .table-container {{
            overflow-y: auto;
            /* Le padding-bottom sera ajusté via JS pour réserver l'espace du container du bas */
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
        }}
        th, td {{
            border: 1px solid black;
            text-align: center;
            padding: 8px;
        }}
        /* La ligne d'en-tête d'origine reste inchangée */
        thead th {{
            background-color: #ADD8E6;
        }}
        td.metric-name {{
            background-color: #ADD8E6;
            font-weight: bold;
            position: relative;
        }}
        td.metric-name small {{
            display: block;
            font-size: 15px;
            font-weight: normal;
        }}
        /* Style pour le rond vert */
        .green-circle {{
            position: absolute;
            top: 2px;
            right: 2px;
            display: inline-block;
            width: 16px;
            height: 16px;
            background-color: green;
            border-radius: 50%;
        }}
        /* Container fixe en bas réduit */
        #bottom-info {{
            position: fixed;
            bottom: 0;
            left: 0;
            right: 0;
            z-index: 50;
            padding: 5px 10px;
            border-top: 1px solid black;
            background-color: #f9f9f9;
            font-size: 0.9em;
        }}
        #bottom-info ul {{
            list-style-type: none;
            padding-left: 0;
            margin: 0;
        }}
        #bottom-info li {{
            margin-bottom: 3px;
        }}
    </style>
</head>
<body>
    <div class="legend">
        <h3>Légende des sections :</h3>
"""
# Génére la légende des sections
for section, color in section_colors.items():
    html_content += f"""
        <div class="legend-item">
            <span class="legend-color" style="background-color:{color}"></span>{section}
        </div>
    """
html_content += f"""
    </div>
    <div class="table-container">
        <table>
            <thead>
                {header_row}
            </thead>
            <tbody>
                {body_rows}
            </tbody>
        </table>
    </div>
    
    <!-- Container fixe pour la légende et l'en-tête cloné en haut -->
    <div id="fixed-header" style="position: fixed; top: 0; left: 0; right: 0; z-index: 50; background-color: white; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
        <!-- Le contenu sera cloné via JavaScript -->
    </div>
    
    <!-- Container fixe en bas -->
    <div id="bottom-info">
        <p><strong>Indicateurs :</strong></p>
        <ul>
            <li><span class="green-circle" style="position: static;"></span> Valeur au-dessus du seuil moyen supérieur du groupe.</li>
            <li style="color: red;">Rouge : Indique une valeur non conforme (inférieure ou supérieure au seuil défini).</li>
        </ul>
        <p><strong>Groupe :</strong> Ensemble des fichiers MetricsOutput des runs passés au TSO500 mis ensemble pour calculer des seuils pour les métriques sans seuils fixés par Illumina.</p>
    </div>
    
    <script>
        function updateTableContainerPadding() {{
            var bottomInfo = document.getElementById('bottom-info');
            var bottomHeight = bottomInfo.offsetHeight;
            // Définir le padding-bottom du conteneur de la table pour que le bas du contenu ne soit pas recouvert
            document.querySelector('.table-container').style.paddingBottom = (bottomHeight + 20) + "px";
        }}

        window.addEventListener('load', function() {{
            // Cloner la légende et l'en-tête du tableau
            var legendOriginal = document.querySelector('.legend');
            var legendClone = legendOriginal.cloneNode(true);
            var table = document.querySelector('.table-container table');
            var theadOriginal = table.querySelector('thead');
            var theadClone = theadOriginal.cloneNode(true);
            
            // Créer un conteneur pour le tableau en-tête cloné
            var tableCloneContainer = document.createElement('div');
            tableCloneContainer.style.overflow = "hidden";
            var clonedTable = document.createElement('table');
            clonedTable.style.borderCollapse = "collapse";
            clonedTable.style.width = table.offsetWidth + "px";
            clonedTable.appendChild(theadClone);
            tableCloneContainer.appendChild(clonedTable);

            // Synchroniser la largeur des cellules
            var originalThs = theadOriginal.querySelectorAll("th");
            var clonedThs = theadClone.querySelectorAll("th");
            for (var i = 0; i < originalThs.length; i++) {{
                var width = window.getComputedStyle(originalThs[i]).width;
                clonedThs[i].style.width = width;
            }}

            var fixedHeader = document.getElementById('fixed-header');
            fixedHeader.appendChild(legendClone);
            fixedHeader.appendChild(tableCloneContainer);

            updateTableContainerPadding();
        }});

        window.addEventListener('resize', updateTableContainerPadding);
    </script>
</body>
</html>
"""

# Sauvegarde le HTML généré
with open(args.output_file, "w") as f:
    f.write(html_content)

print(f"Fichier HTML généré : {args.output_file}")

