#!/usr/bin/env python3
import os
import sys
import pandas as pd

def extraire_valeurs_umi_metrics(csv_file):
    """
    Extrait les valeurs pertinentes du fichier .umi_metrics.csv en prenant
    uniquement la première occurrence de chaque index d'intérêt.
    Index recherchés :
      - Number of reads
      - Families discarded by min-support-reads
      - Consensus pairs emitted
    """
    number_of_reads = None
    families_discarded = None
    consensus_pairs = None

    with open(csv_file, 'r', encoding='utf-8') as f:
        for line in f:
            champs = line.strip().split(',')
            if len(champs) < 4:
                continue
            index_name = champs[2].strip()
            value_str = champs[3].strip()
            try:
                value = int(value_str)
            except ValueError:
                continue
            if index_name == "Number of reads" and number_of_reads is None:
                number_of_reads = value
            elif index_name == "Families discarded by min-support-reads" and families_discarded is None:
                families_discarded = value
            elif index_name == "Consensus pairs emitted" and consensus_pairs is None:
                consensus_pairs = value
            if number_of_reads is not None and families_discarded is not None and consensus_pairs is not None:
                break
    return {
        "number_of_reads": number_of_reads,
        "families_discarded_by_min_support_reads": families_discarded,
        "consensus_pairs_emitted": consensus_pairs
    }

def extraire_valeurs_mapping_metrics(csv_file):
    """
    Extrait les valeurs pertinentes du fichier .mapping_metrics.csv en prenant
    uniquement la première occurrence de chaque index.
    Index recherchés :
      - Total input reads
      - Number of duplicate marked reads
      - Number of unique reads (excl. duplicate marked reads)
    """
    total_input_reads = None
    duplicate_marked_reads = None
    unique_reads = None

    with open(csv_file, 'r', encoding='utf-8') as f:
        for line in f:
            champs = line.strip().split(',')
            if len(champs) < 5:
                continue
            index_name = champs[2].strip()
            value_str = champs[3].strip()
            try:
                value = int(value_str)
            except ValueError:
                continue
            if index_name == "Total input reads" and total_input_reads is None:
                total_input_reads = value
            elif index_name == "Number of duplicate marked reads" and duplicate_marked_reads is None:
                duplicate_marked_reads = value
            elif index_name == "Number of unique reads (excl. duplicate marked reads)" and unique_reads is None:
                unique_reads = value
            if total_input_reads is not None and duplicate_marked_reads is not None and unique_reads is not None:
                break
    return {
        "total_input_reads": total_input_reads,
        "duplicate_marked_reads": duplicate_marked_reads,
        "unique_reads": unique_reads
    }

def trouver_echantillons_adn(racine):
    """
    Recherche dans la racine les dossiers d’échantillons contenant un fichier .umi_metrics.csv.
    En cas de présence d’un dossier 'ADN', on cherche dans ce dossier, sinon directement dans 'racine'.
    Retourne une liste de tuples (nom_echantillon, chemin_du_dossier).
    """
    dossiers_echantillons = []
    chemin_adn = os.path.join(racine, 'ADN')
    
    if os.path.isdir(chemin_adn):
        # On est dans l'organisation ADN/ARN
        dossiers = [os.path.join(chemin_adn, nom) for nom in sorted(os.listdir(chemin_adn)) if os.path.isdir(os.path.join(chemin_adn, nom))]
    else:
        # Recherche directe dans la racine
        dossiers = [os.path.join(racine, nom) for nom in sorted(os.listdir(racine)) if os.path.isdir(os.path.join(racine, nom))]
    
    for dossier in dossiers:
        # On vérifie la présence d'un fichier .umi_metrics.csv
        for f in os.listdir(dossier):
            if f.endswith('.umi_metrics.csv'):
                nom = os.path.basename(dossier)
                dossiers_echantillons.append((nom, dossier))
                break
    return dossiers_echantillons

def trouver_echantillons_arn(racine):
    """
    Recherche dans la racine les dossiers d’échantillons qui :
      - Contiennent un fichier .mapping_metrics.csv
      - Ne contiennent PAS de fichier .umi_metrics.csv (donc pas ADN)
    En cas de présence d’un dossier 'ARN', on cherche dans ce dossier, sinon directement dans 'racine'.
    Retourne une liste de tuples (nom_echantillon, chemin_du_dossier).
    """
    dossiers_echantillons = []
    chemin_arn = os.path.join(racine, 'ARN')
    
    if os.path.isdir(chemin_arn):
        dossiers = [os.path.join(chemin_arn, nom) for nom in sorted(os.listdir(chemin_arn)) if os.path.isdir(os.path.join(chemin_arn, nom))]
    else:
        dossiers = [os.path.join(racine, nom) for nom in sorted(os.listdir(racine)) if os.path.isdir(os.path.join(racine, nom))]
        
    for dossier in dossiers:
        fichiers = os.listdir(dossier)
        # Si un fichier .umi_metrics.csv est présent, ce dossier sera déjà traité en ADN
        if any(f.endswith('.umi_metrics.csv') for f in fichiers):
            continue
        # Si on trouve un fichier .mapping_metrics.csv alors c'est un échantillon ARN
        for f in fichiers:
            if f.endswith('.mapping_metrics.csv'):
                nom = os.path.basename(dossier)
                dossiers_echantillons.append((nom, dossier))
                break
    return dossiers_echantillons

def generer_tableau_adn(echantillons_adn):
    """
    Génère un DataFrame pour les échantillons ADN.
    Lignes générées :
      - Reads bruts ( paired-end ) : Number of reads
      - Paires reads brutes : Number of reads / 2
      - Reads consensus : Consensus pairs emitted
      - Reads uniques filtrés ( family_size < 2 ) : Families discarded by min-support-reads
      - Molécules consensus : (Consensus pairs emitted + Families discarded by min-support-reads)
      - Taux duplicats % : (1 - (Molécules consensus / Paires reads brutes)) * 100 (arrondi à 2 décimales)
    Colonnes : noms des échantillons.
    """
    data_adn = {
        "Reads bruts ( paired-end )": [],
        "Paires reads brutes": [],
        "Reads consensus": [],
        "Reads uniques filtrés ( family_size < 2 )": [],
        "Molécules consensus": [],
        "Taux duplicats %": []
    }
    echantillons = []

    for nom_echantillon, chemin_echantillon in echantillons_adn:
        csv_file = None
        for f in os.listdir(chemin_echantillon):
            if f.endswith(".umi_metrics.csv"):
                csv_file = os.path.join(chemin_echantillon, f)
                break
        if not csv_file:
            continue

        valeurs = extraire_valeurs_umi_metrics(csv_file)
        number_of_reads = valeurs["number_of_reads"] or 0
        families_discarded = valeurs["families_discarded_by_min_support_reads"] or 0
        consensus_pairs = valeurs["consensus_pairs_emitted"] or 0

        reads_bruts = number_of_reads
        paires_reads_brutes = reads_bruts / 2.0
        reads_consensus = consensus_pairs
        reads_uniques_filtres = families_discarded
        molecules_consensus = reads_consensus + reads_uniques_filtres

        taux_duplicats = 0
        if paires_reads_brutes > 0:
            taux_duplicats = (1 - (molecules_consensus / paires_reads_brutes)) * 100
        taux_duplicats = round(taux_duplicats, 2)

        data_adn["Reads bruts ( paired-end )"].append(reads_bruts)
        data_adn["Paires reads brutes"].append(paires_reads_brutes)
        data_adn["Reads consensus"].append(reads_consensus)
        data_adn["Reads uniques filtrés ( family_size < 2 )"].append(reads_uniques_filtres)
        data_adn["Molécules consensus"].append(molecules_consensus)
        data_adn["Taux duplicats %"].append(taux_duplicats)
        echantillons.append(nom_echantillon)

    if not echantillons:
        return pd.DataFrame()

    df_adn = pd.DataFrame(data_adn, columns=data_adn.keys()).transpose()
    df_adn.columns = echantillons
    return df_adn

def generer_tableau_arn(echantillons_arn):
    """
    Génère un DataFrame pour les échantillons ARN.
    Lignes générées :
      - Reads bruts ( paired-end ) : Total input reads
      - Reads marqués comme uniques : Number of unique reads (excl. duplicate marked reads)
      - Reads marqués comme duplicats : Number of duplicate marked reads
      - Taux duplicats % : (Number of duplicate marked reads / Total input reads) * 100
    Colonnes : noms des échantillons.
    """
    data_arn = {
        "Reads bruts ( paired-end )": [],
        "Reads marqués comme uniques": [],
        "Reads marqués comme duplicats": [],
        "Taux duplicats %": []
    }
    echantillons = []

    for nom_echantillon, chemin_echantillon in echantillons_arn:
        csv_file = None
        for f in os.listdir(chemin_echantillon):
            if f.endswith(".mapping_metrics.csv"):
                csv_file = os.path.join(chemin_echantillon, f)
                break
        if not csv_file:
            continue

        valeurs = extraire_valeurs_mapping_metrics(csv_file)
        total_input_reads = valeurs["total_input_reads"] or 0
        duplicate_marked_reads = valeurs["duplicate_marked_reads"] or 0
        unique_reads = valeurs["unique_reads"] or 0

        reads_bruts = total_input_reads
        reads_duplicats = duplicate_marked_reads
        reads_uniques = unique_reads

        taux_duplicats = 0
        if reads_bruts > 0:
            taux_duplicats = (reads_duplicats / reads_bruts) * 100
        taux_duplicats = round(taux_duplicats, 2)

        data_arn["Reads bruts ( paired-end )"].append(reads_bruts)
        data_arn["Reads marqués comme uniques"].append(reads_uniques)
        data_arn["Reads marqués comme duplicats"].append(reads_duplicats)
        data_arn["Taux duplicats %"].append(taux_duplicats)
        echantillons.append(nom_echantillon)

    if not echantillons:
        return pd.DataFrame()

    df_arn = pd.DataFrame(data_arn, columns=data_arn.keys()).transpose()
    df_arn.columns = echantillons
    return df_arn

def main(input_directory):
    """
    Génère un fichier Excel dans QC_Analyse/Taux_duplicats à partir des données.
    Les dossiers contenant un fichier .umi_metrics.csv seront traités dans la feuille ADN
    et ceux sans seront évalués pour la présence d'un fichier .mapping_metrics.csv dans la feuille ARN.
    """
    qc_analyse = os.path.join(input_directory, 'QC_Analyse')
    taux_duplicats_dir = os.path.join(qc_analyse, 'Taux_duplicats')
    os.makedirs(taux_duplicats_dir, exist_ok=True)
    output_excel = os.path.join(taux_duplicats_dir, 'Resultats_duplicats.xlsx')

    # Récupère les échantillons en respectant l'exclusivité
    echantillons_adn = trouver_echantillons_adn(input_directory)
    echantillons_arn = trouver_echantillons_arn(input_directory)

    df_adn = generer_tableau_adn(echantillons_adn)
    df_arn = generer_tableau_arn(echantillons_arn)

    with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
        workbook = writer.book

        # Formats de base
        fmt_base = {
            'font_name': 'Liberation Sans',
            'font_size': 10,
            'border': 1,
            'align': 'center',
            'valign': 'vcenter'
        }
        fmt_empty = workbook.add_format(dict(fmt_base))
        fmt_header = workbook.add_format(dict(fmt_base, bold=True, bg_color='#BFD3C1'))
        fmt_index_normal = workbook.add_format(dict(fmt_base, bg_color='#D6E2E9', text_wrap=True))
        fmt_taux = workbook.add_format(dict(fmt_base, bold=True, bg_color='#F4E3B2', text_wrap=True))
        fmt_cells = workbook.add_format(fmt_base)

        row_height = 42.5
        col_width = 26.4

        ######################
        # Feuille ADN
        if not df_adn.empty:
            df_adn.to_excel(writer, sheet_name='ADN', header=True, index=True)
            ws_adn = writer.sheets['ADN']
            ws_adn.set_column(0, 0, col_width)
            for col in range(1, len(df_adn.columns) + 1):
                ws_adn.set_column(col, col, col_width)
            for row in range(len(df_adn.index) + 1):
                ws_adn.set_row(row, row_height)
                for col in range(len(df_adn.columns) + 1):
                    if row == 0 and col == 0:
                        ws_adn.write(row, col, "", fmt_empty)
                    elif row == 0:
                        ws_adn.write(row, col, df_adn.columns[col - 1], fmt_header)
                    elif col == 0:
                        if df_adn.index[row - 1] == "Taux duplicats %":
                            ws_adn.write(row, col, df_adn.index[row - 1], fmt_taux)
                        else:
                            ws_adn.write(row, col, df_adn.index[row - 1], fmt_index_normal)
                    else:
                        if df_adn.index[row - 1] == "Taux duplicats %":
                            ws_adn.write(row, col, df_adn.iloc[row - 1, col - 1], fmt_taux)
                        else:
                            ws_adn.write(row, col, df_adn.iloc[row - 1, col - 1], fmt_cells)

        ######################
        # Feuille ARN
        if not df_arn.empty:
            df_arn.to_excel(writer, sheet_name='ARN', header=True, index=True)
            ws_arn = writer.sheets['ARN']
            ws_arn.set_column(0, 0, col_width)
            for col in range(1, len(df_arn.columns) + 1):
                ws_arn.set_column(col, col, col_width)
            for row in range(len(df_arn.index) + 1):
                ws_arn.set_row(row, row_height)
                for col in range(len(df_arn.columns) + 1):
                    if row == 0 and col == 0:
                        ws_arn.write(row, col, "", fmt_empty)
                    elif row == 0:
                        ws_arn.write(row, col, df_arn.columns[col - 1], fmt_header)
                    elif col == 0:
                        if df_arn.index[row - 1] == "Taux duplicats %":
                            ws_arn.write(row, col, df_arn.index[row - 1], fmt_taux)
                        else:
                            ws_arn.write(row, col, df_arn.index[row - 1], fmt_index_normal)
                    else:
                        if df_arn.index[row - 1] == "Taux duplicats %":
                            ws_arn.write(row, col, df_arn.iloc[row - 1, col - 1], fmt_taux)
                        else:
                            ws_arn.write(row, col, df_arn.iloc[row - 1, col - 1], fmt_cells)

    print(f"Fichier Excel généré : {output_excel}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : python duplicats.py <chemin_vers_dossier_contenant_echantillons>")
        sys.exit(1)
    input_directory = sys.argv[1]
    main(input_directory)

