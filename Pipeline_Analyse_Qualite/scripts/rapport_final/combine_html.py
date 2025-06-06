#!/usr/bin/env python3
import os
import json
import re
import sys

def make_options(files):
    """
    Génère une chaîne HTML des <option> pour chaque échantillon extrait de `files`.
    """
    lines = []
    for f in sorted(files):
        sample = os.path.splitext(f)[0].split('_')[0]
        key = sample
        lines.append(f'            <option value="{key}">{sample}</option>')
    return "\n".join(lines)


def generate_index(directory='.'):    
    """
    Produit un Rapport_QC.html avec quatre onglets :
      - Informations qualité run : MetricsComparison.html
      - Variants : fichiers CombinedVariantOutput*.html
      - Low coverage 50X : fichiers *_coverage.html
      - Détails duplicats : Resultats_duplicats.html
    Les fichiers .html doivent être dans le même dossier que l'index.
    """
    src = directory

    # 1) Fichiers fixes
    metrics_file   = 'MetricsComparison.html'
    duplicats_file = 'Resultats_duplicats.html'

    # 2) Collecte de tous les .html présents
    all_files = [f for f in os.listdir(src) if f.endswith('.html')]

    # 3) Filtrage pour chaque onglet
    variant_files = [f for f in all_files if 'CombinedVariantOutput' in f]
    lowcov_files  = [f for f in all_files if f.endswith('_coverage.html')]

    # 4) Construction du dict pages
    pages = {
        'MetricsComparison': metrics_file,
        'Duplicats': duplicats_file
    }
    # ajouter variants (clé = sample)
    for f in variant_files:
        sample = os.path.splitext(f)[0].split('_')[0]
        pages[sample] = f
    # ajouter low coverage (clé = sample)
    for f in lowcov_files:
        sample = os.path.splitext(f)[0].split('_')[0]
        pages[sample] = f

    # 5) Génération HTML
    html = f"""<!DOCTYPE html>
<html lang=\"fr\">
<head>
  <meta charset=\"utf-8\">
  <title>Rapport QC</title>
  <style>
    body {{ margin:0; padding:0; font-family:Arial,Helvetica,sans-serif }}
    #contentFrame {{ width:100%; height:calc(100vh - 80px); border:none }}
    #container {{ position:fixed; bottom:0; width:100%; background:#f0f0f0;
                   padding:15px; text-align:center; box-shadow:0 -2px 5px rgba(0,0,0,.2) }}
    button, select {{ font-size:16px; font-weight:bold; padding:10px 20px; margin:5px;
                     border:0; border-radius:5px; background:#d3d3d3; cursor:pointer }}
    button:hover, select:hover {{ background:#c0c0c0 }}
    select {{ appearance:button; text-align-last:center; min-width:150px }}
  </style>
</head>
<body>
  <!-- Volet d'affichage -->
  <iframe id=\"contentFrame\" src=\"{metrics_file}\"></iframe>

  <!-- Barre de navigation -->
  <div id=\"container\">
    <!-- Bouton Infos QC -->
    <button data-page=\"MetricsComparison\">Informations qualité run</button>

    <!-- Menu Variants -->
    <select id=\"variantSelect\">
      <option value=\"\">Variants</option>
{make_options(variant_files)}
    </select>

    <!-- Menu Low coverage -->
    <select id=\"lowCoverageSelect\">
      <option value=\"\">Low coverage 50X</option>
{make_options(lowcov_files)}
    </select>

    <!-- Bouton Duplicats -->
    <button data-page=\"Duplicats\">Détails duplicats</button>
  </div>

  <script>
    // Mapping clé->fichier
    const pages = {json.dumps(pages, ensure_ascii=False)};
    const frame = document.getElementById('contentFrame');

    function showPage(key) {{
      const url = pages[key];
      if (url) frame.src = url;
    }}

    // Evenements
    document.querySelectorAll('button[data-page]').forEach(btn =>
      btn.addEventListener('click', () => showPage(btn.dataset.page))
    );

    document.getElementById('variantSelect').addEventListener('change', e => {{
      if (e.target.value) showPage(e.target.value);
      e.target.selectedIndex = 0;
    }});

    document.getElementById('lowCoverageSelect').addEventListener('change', e => {{
      if (e.target.value) showPage(e.target.value);
      e.target.selectedIndex = 0;
    }});
  </script>
</body>
</html>
"""

    # 6) Écriture du fichier final
    output_file = os.path.join(directory, 'Rapport_QC.html')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
    print('✅ Fichier index généré :', output_file)


if __name__ == '__main__':
    base_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    generate_index(base_dir)

