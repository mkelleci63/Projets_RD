#!/usr/bin/env python3
"""
variants_infos_final.py:
  À partir des 3 MAFs bruts dans un dossier donné,
  extrait les variants :
    1) présents dans FFPE mais absents de Leuco, puis
    2) parmi ceux-ci, présents également dans ctDNA/cfDNA.

  Produit un TSV avec un en-tête en deux lignes :
    - 1re ligne : colonnes d'origine vides, puis noms d'échantillons (FFPE, ctDNA)
      au-dessus de chacune de leurs deux colonnes (Variant_Frequencies, Total_Depth)
    - 2e ligne : noms des colonnes d'origine (sans Variant_Classification), puis
      Variant_Frequencies et Total_Depth pour chaque échantillon

Usage :
    python variants_infos_final.py -i /chemin/vers/dossier_input

Le fichier de sortie 'ffpe_ctdna_matrix.tsv' sera placé dans
/chemin/vers/dossier_input/variants_informatifs/.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--input_dir', required=True,
                   help="Dossier contenant les fichiers MAF")
    return p.parse_args()

def read_maf(path: Path) -> pd.DataFrame:
    cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'HGVSc', 'HGVSp', 'Consequence',
        'Variant_Frequencies', 'Total_Depth'
    ]
    df = pd.read_csv(
        path, sep='\t', comment='#', dtype=str, low_memory=False,
        encoding='latin-1', on_bad_lines='warn'
    )
    missing = set(cols) - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} : colonnes manquantes {sorted(missing)}")
    return df[cols].copy()

def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"Erreur : {input_dir} n’est pas un dossier valide.", file=sys.stderr)
        sys.exit(1)

    # repérer par mot-clé (insensible à la casse), inclut 'cfdna' pour cfDNA
    ffpe_f = leuco_f = ctDNA_f = None
    for path in input_dir.glob("*.maf"):
        name = path.name.lower()
        if 'ffpe' in name:
            ffpe_f = path
        elif 'leuco' in name:
            leuco_f = path
        elif 'ctdna' in name or 'cfdna' in name:
            ctDNA_f = path
    if not all([ffpe_f, leuco_f, ctDNA_f]):
        print("Erreur : fichiers FFPE, Leuco ou ctDNA/cfDNA manquants.", file=sys.stderr)
        sys.exit(1)

    # lecture
    df_ffpe  = read_maf(ffpe_f)
    df_leuco = read_maf(leuco_f)
    df_ctDNA = read_maf(ctDNA_f)

    # colonnes d'origine (sans Variant_Classification) et clés d'identité
    orig_cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'HGVSc', 'HGVSp', 'Consequence'
    ]
    key_cols = ['HGVSc']
    extra_cols = ['Variant_Frequencies','Total_Depth']

    # filtrages
    keys_ffpe  = set(map(tuple, df_ffpe[key_cols].values))
    keys_leuco = set(map(tuple, df_leuco[key_cols].values))
    keys_ctDNA = set(map(tuple, df_ctDNA[key_cols].values))

    ffpe_only = keys_ffpe - keys_leuco
    final_keys = ffpe_only & keys_ctDNA

    # sous-ensemble FFPE
    mask_ffpe = df_ffpe[key_cols].apply(tuple, axis=1).isin(final_keys)
    df_sel = df_ffpe[mask_ffpe].reset_index(drop=True)

    # récupérer fréquences/profondeurs ctDNA pour ces variants
    df_ct = df_ctDNA[key_cols + extra_cols].copy()
    df_ct = df_ct.rename(columns={
        'Variant_Frequencies': 'ctDNA_Variant_Frequencies',
        'Total_Depth':        'ctDNA_Total_Depth'
    })
    # fusionner sur clés
    df_final = df_sel.merge(df_ct, on=key_cols, how='left')

    # préparation output
    out_dir = input_dir / 'variants_informatifs'
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / 'ffpe_ctdna_matrix.tsv'

    # construction des en-têtes
    samples = ['FFPE','ctDNA']
    first_header = [''] * len(orig_cols)
    for s in samples:
        first_header += [s, s]
    second_header = orig_cols.copy()
    for _ in samples:
        second_header += ['Variant_Frequencies','Total_Depth']

    # écriture
    with open(out_file, 'w', newline='') as fh:
        fh.write('\t'.join(first_header) + '\n')
        fh.write('\t'.join(second_header) + '\n')
        df_final.to_csv(fh, sep='\t', header=False, index=False)

    print(f"Écrit {len(df_final)} variants dans {out_file}")

if __name__ == '__main__':
    main()
