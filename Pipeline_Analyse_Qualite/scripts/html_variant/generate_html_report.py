#!/usr/bin/env python3
import argparse
import os
import csv

# Ensemble des sections spécifiques à l'ARN (en minuscules)
ARN_SECTIONS = {"splice variants", "fusions"}

def html_escape(s: str) -> str:
    """ Protège les caractères spéciaux pour éviter qu'ils ne soient interprétés comme des balises HTML. """
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace('"', "&quot;"))

def parse_variant_file(content):
    """
    Découpe le fichier TSV en sections, en utilisant les lignes d'en-tête
    au format [NomSection]. Chaque section est un dictionnaire avec 'name' et 'lines'.
    Les lignes contenant certains textes indésirables sont ignorées.
    """
    unwanted_phrases = [
        "TUMOR FRACTION AND PLOIDY ARE BETA FEATURES. REFER TO USER GUIDE FOR DETAILS.",
        "ABSOLUTE COPY NUMBER AND LOSS OF HETEROZYGOSITY ARE BETA FEATURES. REFER TO USER GUIDE FOR DETAILS.",
        "ABSOLUTE COPY NUMBER IS A BETA FEATURE. REFER TO USER GUIDE FOR DETAILS.",
        "TUMOR FRACTION, PLOIDY, ABSOLUTE COPY NUMBER AND LOSS OF HETEROZYGOSITY ARE BETA FEATURES. REFER TO USER GUIDE FOR DETAILS."
    ]
    lines = content.splitlines()
    sections, current = [], None
    for line in lines:
        t = line.strip()
        if not t or any(phrase in t for phrase in unwanted_phrases):
            continue
        if t.startswith('[') and t.endswith(']'):
            if current and current["name"].lower() != "default":
                sections.append(current)
            current = {"name": t[1:-1].strip(), "lines": []}
        else:
            if current is None:
                current = {"name": "Default", "lines": []}
            current["lines"].append(line)
    if current and current["name"].lower() != "default":
        sections.append(current)
    return sections

def parse_analysis_details(sections):
    """ Recherche la section [Analysis Details] et construit un dictionnaire pour les paires clé/valeur. """
    details = {}
    for sec in sections:
        if sec["name"].lower() == "analysis details":
            for line in sec["lines"]:
                parts = [p.strip() for p in line.split('\t') if p.strip()]
                i = 0
                while i < len(parts) - 1:
                    key = parts[i]
                    val = parts[i+1]
                    details[key] = val
                    i += 2
    return details

def update_copy_number_lines(lines, tumor_percentage):
    """
    Met à jour (ou crée) la colonne 'Absolute Copy Number' en se basant sur la formule :
        n = ((200 * Y) - 2 * (100 - X)) / X
    où Y est la valeur de la colonne 'Fold Change' et X = '%_tumeur'.
    """
    try:
        X = float(tumor_percentage)
    except ValueError:
        print(f"Erreur : '%_tumeur' non numérique : {tumor_percentage}")
        return lines

    if not lines:
        return lines

    header = lines[0].rstrip("\n").split("\t")
    try:
        idx_fold = header.index("Fold Change")
    except ValueError:
        print("Erreur : La colonne 'Fold Change' est absente dans le header => impossible de calculer.")
        return lines

    # S'assurer qu'il y a la colonne "Absolute Copy Number"
    try:
        idx_abs = header.index("Absolute Copy Number")
    except ValueError:
        header.append("Absolute Copy Number")
        idx_abs = len(header) - 1

    new_lines = ["\t".join(header)]
    for line in lines[1:]:
        line = line.rstrip("\n")
        if not line.strip():
            new_lines.append(line)
            continue
        cells = line.split("\t")
        if len(cells) < len(header):
            cells += [""] * (len(header) - len(cells))
        try:
            fold_val = float(cells[idx_fold])
            n = ((200 * fold_val) - 2 * (100 - X)) / X
            new_val = f"{n:.2f}"
        except ValueError:
            new_val = cells[idx_abs] if idx_abs < len(cells) else ""
        cells[idx_abs] = new_val
        new_lines.append("\t".join(cells))
    return new_lines

def remove_empty_columns_from_table(lines):
    """ Supprime les colonnes vides du tableau (TSV). Retourne un header et les données filtrées. """
    if not lines:
        return [], []
    header = [cell.strip() for cell in lines[0].split("\t")]
    data = []
    for line in lines[1:]:
        cells = [cell.strip() for cell in line.split("\t")]
        if len(cells) < len(header):
            cells.extend([""] * (len(header) - len(cells)))
        data.append(cells)
    keep = []
    for i, h in enumerate(header):
        if not (h == "" and all(row[i] == "" for row in data)):
            keep.append(i)
    new_header = [header[i] for i in keep]
    new_data = [[row[i] for i in keep] for row in data]
    return new_header, new_data

def generate_table_from_lines(lines):
    """ Génère un tableau HTML simple à partir de lignes TSV. """
    header, data = remove_empty_columns_from_table(lines)
    html = "<table style='max-width:100%;'>\n<thead><tr>"
    for col in header:
        html += f"<th>{html_escape(col)}</th>"
    html += "</tr></thead>\n<tbody>\n"
    for row in data:
        html += "<tr>"
        for cell in row:
            html += f"<td>{html_escape(cell)}</td>"
        html += "</tr>\n"
    html += "</tbody></table>\n"
    return html

def generate_small_variants_table_html(lines):
    """
    Génère le tableau HTML pour la section "Small Variants",
    avec tri possible sur "Allele Frequency" et "Depth".
    """
    header, data = remove_empty_columns_from_table(lines)
    header_html = "<tr>\n"
    for i, cell in enumerate(header):
        cl = cell.lower()
        cell_escaped = html_escape(cell)
        if cl == "allele frequency":
            header_html += (
                f'<th data-col-index="{i}">{cell_escaped} '
                f'<span class="sort-icon" onclick="sortTable(\'allele frequency\', \'asc\')">&#9650;</span> '
                f'<span class="sort-icon" onclick="sortTable(\'allele frequency\', \'desc\')">&#9660;</span>'
                "</th>\n"
            )
        elif cl == "depth":
            header_html += (
                f'<th data-col-index="{i}">{cell_escaped} '
                f'<span class="sort-icon" onclick="sortTable(\'depth\', \'asc\')">&#9650;</span> '
                f'<span class="sort-icon" onclick="sortTable(\'depth\', \'desc\')">&#9660;</span>'
                "</th>\n"
            )
        else:
            header_html += f'<th data-col-index="{i}">{cell_escaped}</th>\n'
    header_html += "</tr>\n"
    
    body_html = ""
    for row in data:
        body_html += "<tr>\n"
        for cell in row:
            body_html += f"<td>{html_escape(cell)}</td>\n"
        body_html += "</tr>\n"
    table_html = f"""
<table id="smallVariantsTable">
  <thead>{header_html}</thead>
  <tbody>{body_html}</tbody>
</table>
"""
    search_bar = (
        '<input type="text" id="smallVariantsSearchInput" placeholder="Rechercher un gène..." '
        'onkeyup="filterTableRows(\'smallVariantsSearchInput\', \'smallVariantsTable\', 0)" '
        'style="width:100%; padding:8px; margin-bottom:10px; border:1px solid #ccc; border-radius:4px;">'
    )
    return f'<div>{search_bar}<div class="scroll-wrapper" style="overflow-x:auto; overflow-y:auto; max-height:400px; width:100%;">{table_html}</div></div>'

def generate_table_copy_number_variants(lines, gene_cytoband):
    """
    Génère le tableau HTML pour la section [Copy Number Variants],
    en affichant 7 colonnes dans l'ordre :
      Gene, Chromosome/scaffold, Karyotype, Strand, Fold Change, Copy Number Variant, Absolute Copy Number.
    
    Améliorations :
      - Convertit "1" en "+" et "-1" en "-" pour la colonne Strand
      - Ajout d'icônes de tri sur Chromosome/scaffold, Karyotype, Fold Change, Absolute Copy Number
      - Le tri sur Chromosome/scaffold prend en compte X, Y après les autosomes
      - Le tri sur Karyotype prend d’abord p < q, puis un tri numérique sur ce qui suit
    """
    header = lines[0].rstrip("\n").split("\t")
    data = [line.rstrip("\n").split("\t") for line in lines[1:]]

    # On identifie les index de chaque colonne
    try:
        idx_gene = header.index("Gene")
    except ValueError:
        idx_gene = None
    try:
        idx_fold = header.index("Fold Change")
    except ValueError:
        idx_fold = None
    try:
        idx_cnv = header.index("Copy Number Variant")
    except ValueError:
        idx_cnv = None
    try:
        idx_abs = header.index("Absolute Copy Number")
    except ValueError:
        idx_abs = None
    # Strand
    try:
        idx_strand = header.index("Strand")
    except ValueError:
        idx_strand = None

    # On construit le nouveau header en y ajoutant les icônes de tri
    # 0 => Gene, 1 => Chromosome/scaffold, 2 => Karyotype, 3 => Strand, 4 => Fold Change, 5 => Copy Number Variant, 6 => Absolute Copy Number
    new_header = []
    col_names = [
        "Gene",
        "Chromosome/scaffold",
        "Karyotype",
        "Strand",
        "Fold Change",
        "Copy Number Variant",
        "Absolute Copy Number"
    ]
    # On va créer un dictionnaire pour associer le nom de la colonne à son index dans new_header
    col_index_map = {}
    for i, cn in enumerate(col_names):
        col_index_map[cn] = i
    # On ajoute des icônes de tri sur Chromosome/scaffold, Karyotype, Fold Change, et Absolute Copy Number
    def col_with_sort_icons(col_name):
        # On autorise le tri sur 4 colonnes
        if col_name.lower() in ["chromosome/scaffold", "fold change", "absolute copy number", "karyotype"]:
            return (
                f'<th data-col-name="{col_name}">{html_escape(col_name)} '
                f'<span class="sort-icon" onclick="sortCnvTable(\'{col_name}\', \'asc\')">&#9650;</span> '
                f'<span class="sort-icon" onclick="sortCnvTable(\'{col_name}\', \'desc\')">&#9660;</span>'
                '</th>'
            )
        else:
            return f'<th data-col-name="{col_name}">{html_escape(col_name)}</th>'
    # On construit la ligne d'entête
    header_html = "<tr>\n"
    for col_name in col_names:
        header_html += col_with_sort_icons(col_name) + "\n"
    header_html += "</tr>\n"

    # On parcourt les lignes pour construire new_data
    new_data = []
    for row in data:
        # Récupération du gène
        gene = row[idx_gene] if idx_gene is not None and idx_gene < len(row) else ""
        # Récupération des infos du gène depuis le fichier gene_cytoband
        gene_details = gene_cytoband.get(gene, {"Chromosome/scaffold": "", "Karyotype": "", "Strand": ""})

        # Conversion du brin (1 => +, -1 => -)
        actual_strand = gene_details["Strand"]
        if actual_strand == "1":
            actual_strand = "+"
        elif actual_strand == "-1":
            actual_strand = "-"
        # Si la ligne TSV a un champ Strand, on l'écrase ou on l'ignore
        if idx_strand is not None and idx_strand < len(row):
            # On prend la valeur la plus pertinente
            # (GeneDetails en priorité, sinon la valeur du TSV)
            row_strand = row[idx_strand].strip()
            if row_strand in ["1", "+"]:
                row_strand = "+"
            elif row_strand in ["-1", "-"]:
                row_strand = "-"
            # On fusionne
            if actual_strand.strip():
                strand = actual_strand
            else:
                strand = row_strand
        else:
            strand = actual_strand

        fold_val = row[idx_fold] if idx_fold is not None and idx_fold < len(row) else ""
        cnv_val  = row[idx_cnv] if idx_cnv is not None and idx_cnv < len(row) else ""
        abs_val  = row[idx_abs] if idx_abs is not None and idx_abs < len(row) else ""

        # Nouvelle ligne
        new_row = [
            gene,
            gene_details["Chromosome/scaffold"],  # Col 1
            gene_details["Karyotype"],            # Col 2
            strand,                               # Col 3
            fold_val,                             # Col 4
            cnv_val,                              # Col 5
            abs_val                               # Col 6
        ]
        new_data.append(new_row)

    # Construction du tableau HTML
    html_table = f"<table id='cnvTable'>\n<thead>\n{header_html}</thead>\n<tbody>\n"
    for row in new_data:
        html_table += "<tr>\n"
        for cell in row:
            html_table += f"<td>{html_escape(cell)}</td>\n"
        html_table += "</tr>\n"
    html_table += "</tbody>\n</table>\n"

    # Barre de recherche
    search_bar = (
        '<input type="text" id="cnvSearchInput" placeholder="Rechercher un gène..." '
        'onkeyup="filterTableRows(\'cnvSearchInput\', \'cnvTable\', 0)" '
        'style="width:100%; padding:8px; margin-bottom:10px; border:1px solid #ccc; border-radius:4px;">'
    )

    # On englobe le tout
    return f'<div>{search_bar}<div class="scroll-wrapper" style="overflow-x:auto; overflow-y:auto; max-height:400px; width:100%;">{html_table}</div></div>'

def generate_table_with_search_for_loh(lines):
    """ Génère le tableau HTML pour la section [Loss of Heterozygosity], avec barre de recherche. """
    header, data = remove_empty_columns_from_table(lines)
    html_table = "<table id='lohTable'>\n<thead>\n<tr>\n"
    for col in header:
        html_table += f"<th>{html_escape(col)}</th>\n"
    html_table += "</tr>\n</thead>\n<tbody>\n"
    for row in data:
        html_table += "<tr>\n"
        for cell in row:
            html_table += f"<td>{html_escape(cell)}</td>\n"
        html_table += "</tr>\n"
    html_table += "</tbody>\n</table>\n"
    search_bar = (
        '<input type="text" id="lohSearchInput" placeholder="Rechercher un gène..." '
        'onkeyup="filterTableRows(\'lohSearchInput\', \'lohTable\', 0)" '
        'style="width:100%; padding:8px; margin-bottom:10px; border:1px solid #ccc; border-radius:4px;">'
    )
    return f'<div>{search_bar}<div class="scroll-wrapper" style="overflow-x:auto; overflow-y:auto; max-height:400px; width:100%;">{html_table}</div></div>'

def generate_section_html(section, gene_cytoband=None):
    """
    Génère le HTML pour une section sous forme d'accordéon.
    Ajout du paramètre gene_cytoband afin de personnaliser "Copy Number Variants".
    """
    sec_name = section["name"]
    lower_name = sec_name.lower()
    if lower_name == "small variants":
        content = generate_small_variants_table_html(section["lines"])
    elif lower_name == "copy number variants" and gene_cytoband is not None:
        content = generate_table_copy_number_variants(section["lines"], gene_cytoband)
    elif lower_name == "loss of heterozygosity":
        content = generate_table_with_search_for_loh(section["lines"])
    else:
        if section["lines"] and "\t" in section["lines"][0]:
            content = generate_table_from_lines(section["lines"])
        else:
            content = "<pre>" + "\n".join(html_escape(l) for l in section["lines"]) + "</pre>"
    html = f"""
<div class="accordion-section">
  <div class="accordion-header" onclick="toggleSection(this)">{html_escape(sec_name)}</div>
  <div class="accordion-content">{content}</div>
</div>
"""
    return html

def generate_html_two_blocks(sections, sample_name, dna_id, rna_id, gene_cytoband):
    """ Construit un document HTML avec deux blocs : ADN et ARN + tri custom pour CNV. """
    dna_sections = []
    rna_sections = []
    for sec in sections:
        lower_name = sec["name"].lower()
        if lower_name in {"analysis details", "sequencing run details"}:
            continue
        if lower_name in ARN_SECTIONS:
            rna_sections.append(sec)
        else:
            dna_sections.append(sec)
    
    dna_block = ""
    if dna_id.upper() != "NA" and dna_sections:
        dna_block += '<h2>ADN</h2>\n'
        for sec in dna_sections:
            dna_block += generate_section_html(sec, gene_cytoband)
    rna_block = ""
    if rna_id.upper() != "NA" and rna_sections:
        rna_block += '<h2>ARN</h2>\n'
        for sec in rna_sections:
            rna_block += generate_section_html(sec, gene_cytoband)
    
    if dna_id.upper() != "NA" and rna_id.upper() != "NA":
        suffix = " (ADN + ARN)"
    elif dna_id.upper() != "NA":
        suffix = " (ADN)"
    elif rna_id.upper() != "NA":
        suffix = " (ARN)"
    else:
        suffix = ""
    
    final_title = f"Échantillon : {sample_name}{suffix}"
    
    html_template = f"""<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8">
  <title>{html_escape(final_title)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; background-color: #f7f7f7; margin: 0; padding: 20px; }}
    .central-container {{ max-width: 90%; margin: auto; }}
    h1 {{ text-align: center; }}
    h2 {{ margin-top: 30px; color: #333; }}
    .accordion-section {{ border: 1px solid #B0C4DE; margin-bottom: 10px; border-radius: 4px; background-color: #fff; }}
    .accordion-header {{ background-color: #B0C4DE; padding: 10px; cursor: pointer; font-weight: bold; }}
    .accordion-content {{ padding: 10px; border-top: 1px solid #B0C4DE; display: none; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th, td {{ padding: 6px; border: 1px solid #aaa; text-align: left; }}
    .scroll-wrapper {{ overflow-x: auto; overflow-y: auto; max-height: 400px; width: 100%; border: 1px solid #ddd; }}
    input[type="text"] {{
      width: 100%;
      padding: 8px;
      margin-bottom: 10px;
      border: 1px solid #ccc;
      border-radius: 4px;
      box-sizing: border-box;
    }}
    .sort-icon {{
      margin-left: 4px;
      cursor: pointer;
    }}
  </style>
</head>
<body>
  <div class="central-container">
    <h1>{html_escape(final_title)}</h1>
    {dna_block}
    {rna_block}
  </div>
  <script>
    function toggleSection(header) {{
      var content = header.nextElementSibling;
      if (content.style.display === "none" || content.style.display === "") {{
        content.style.display = "block";
      }} else {{
        content.style.display = "none";
      }}
    }}

    // Tri pour la table smallVariantsTable
    function getColumnIndex(columnName, tableId) {{
      var headers = document.querySelectorAll("#" + tableId + " thead th");
      for (var i = 0; i < headers.length; i++) {{
        if (headers[i].textContent.trim().toLowerCase().includes(columnName)) {{
          return parseInt(headers[i].getAttribute("data-col-index")) || i;
        }}
      }}
      return -1;
    }}
    function sortTable(columnName, direction) {{
      var table = document.getElementById("smallVariantsTable");
      var tbody = table.tBodies[0];
      var rows = Array.from(tbody.rows);
      var colIndex = getColumnIndex(columnName, "smallVariantsTable");
      if (colIndex === -1) return;
      rows.sort(function(a, b) {{
        var aText = a.cells[colIndex].textContent.trim();
        var bText = b.cells[colIndex].textContent.trim();
        var aVal = parseFloat(aText);
        var bVal = parseFloat(bText);
        // Si ce n'est pas numérique, on compare en string
        if (isNaN(aVal) || isNaN(bVal)) {{
          aVal = aText.toLowerCase();
          bVal = bText.toLowerCase();
        }}
        if (direction === 'asc') return (aVal > bVal) ? 1 : -1;
        else return (aVal < bVal) ? 1 : -1;
      }});
      rows.forEach(function(row) {{
        tbody.appendChild(row);
      }});
    }}
    function filterTableRows(inputId, tableId, colIndex) {{
      var input = document.getElementById(inputId);
      var filter = input.value.toUpperCase();
      var table = document.getElementById(tableId);
      var tbody = table.tBodies[0];
      var rows = tbody.getElementsByTagName("tr");
      for (var i = 0; i < rows.length; i++) {{
        var cell = rows[i].getElementsByTagName("td")[colIndex];
        if (cell) {{
          var txtValue = cell.textContent || cell.innerText;
          rows[i].style.display = txtValue.toUpperCase().indexOf(filter) > -1 ? "" : "none";
        }}
      }}
    }}

    // Tri pour la table cnvTable
    // Les colonnes sur lesquelles on veut trier spécialement : Chromosome/scaffold, Karyotype, Fold Change, Absolute Copy Number
    function parseChromosome(value) {{
      // Enlève "chr" s'il existe
      value = value.replace(/^chr/i, '');
      // Convertit X => 23, Y => 24, M => 25, etc.
      if (value.toUpperCase() === 'X') return 23;
      if (value.toUpperCase() === 'Y') return 24;
      // Sinon on essaie de parser un nombre
      var num = parseInt(value);
      if (!isNaN(num)) return num;
      // Sinon => on le place en dernier
      return 99999;
    }}
    function parseKaryotype(value) {{
      // ex : p15.33, q22.2
      // on sépare le bras (p/q) et la partie numérique
      // bras p => 0, bras q => 1
      // num => ex 15.33 => 15.33
      value = value.trim();
      if (!value) return [999, 999]; // s'il n'y a rien
      var bras = value.charAt(0).toLowerCase() === 'p' ? 0 : 1;  // p => 0, q => 1
      var remain = value.substring(1); // ex : 15.33
      // on parse remain => float
      var num = parseFloat(remain) || 0;
      return [bras, num];
    }}
    function sortCnvTable(colName, direction) {{
      var table = document.getElementById("cnvTable");
      if (!table) return;
      var tbody = table.tBodies[0];
      var rows = Array.from(tbody.rows);

      // On récupère l'index de colonne
      var headers = table.querySelectorAll("thead th");
      var colIndex = -1;
      for (var i = 0; i < headers.length; i++) {{
        if (headers[i].getAttribute("data-col-name") && headers[i].getAttribute("data-col-name").toLowerCase() === colName.toLowerCase()) {{
          colIndex = i;
          break;
        }}
      }}
      if (colIndex === -1) return;

      rows.sort(function(a, b) {{
        var aText = a.cells[colIndex].textContent.trim();
        var bText = b.cells[colIndex].textContent.trim();

        if (colName.toLowerCase() === 'fold change' || colName.toLowerCase() === 'absolute copy number') {{
          // Tri numérique
          var aVal = parseFloat(aText);
          var bVal = parseFloat(bText);
          if (isNaN(aVal)) aVal = -9999999;
          if (isNaN(bVal)) bVal = -9999999;
          if (direction === 'asc') return (aVal > bVal) ? 1 : -1;
          else return (aVal < bVal) ? 1 : -1;
        }} else if (colName.toLowerCase() === 'chromosome/scaffold') {{
          // Tri personnalisé
          var aVal = parseChromosome(aText);
          var bVal = parseChromosome(bText);
          if (aVal === bVal) {{
            // Si tie, on compare la karyotype (colIndex=2)
            var aKaryo = a.cells[2].textContent.trim();
            var bKaryo = b.cells[2].textContent.trim();
            var A = parseKaryotype(aKaryo);
            var B = parseKaryotype(bKaryo);
            // Compare bras puis numéro
            if (A[0] === B[0]) {{
              return (direction === 'asc') ? (A[1] - B[1]) : (B[1] - A[1]);
            }} else {{
              return (direction === 'asc') ? (A[0] - B[0]) : (B[0] - A[0]);
            }}
          }}
          if (direction === 'asc') return (aVal > bVal) ? 1 : -1;
          else return (aVal < bVal) ? 1 : -1;
        }} else if (colName.toLowerCase() === 'karyotype') {{
          // Tri perso p < q, puis num
          var A = parseKaryotype(aText);
          var B = parseKaryotype(bText);
          if (A[0] === B[0]) {{
            if (direction === 'asc') return A[1] - B[1];
            else return B[1] - A[1];
          }} else {{
            if (direction === 'asc') return A[0] - B[0];
            else return B[0] - A[0];
          }}
        }}
        // Sinon tri alpha
        var aVal = aText.toLowerCase();
        var bVal = bText.toLowerCase();
        if (direction === 'asc') return (aVal > bVal) ? 1 : -1;
        else return (aVal < bVal) ? 1 : -1;
      }});

      rows.forEach(function(row) {{
        tbody.appendChild(row);
      }});
    }}
  </script>
</body>
</html>
"""
    return html_template

def load_metadata(metadata_file):
    """ Charge le fichier CSV de métadonnées et retourne un dictionnaire associant le Patient à la valeur de %_tumeur. """
    metadata = {}
    try:
        with open(metadata_file, newline='', encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                key = row["Patient"].strip()
                metadata[key] = row["%_tumeur"].strip()
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier de métadonnées : {e}")
    return metadata

def load_gene_cytoband(cytoband_file):
    """
    Charge le fichier gene_cytoband.txt et retourne un dictionnaire
    où la clé est le nom du gène et la valeur est un dictionnaire contenant :
      - 'Chromosome/scaffold'
      - 'Karyotype'
      - 'Strand'
    """
    gene_info = {}
    try:
        with open(cytoband_file, newline='', encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gene = row["Gene"].strip()
                gene_info[gene] = {
                    "Chromosome/scaffold": row["Chromosome/scaffold"].strip(),
                    "Karyotype": row["Karyotype"].strip(),
                    "Strand": row["Strand"].strip()
                }
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier gene_cytoband.txt : {e}")
    return gene_info

def extract_sample_name(tsv_filepath):
    """
    Extrait le nom de l'échantillon à partir du nom du fichier TSV.
    On suppose que le fichier est nommé {sample}_CombinedVariantOutput.tsv.
    """
    base = os.path.basename(tsv_filepath)
    parts = base.split("_")
    return parts[0] if parts else base

def main():
    parser = argparse.ArgumentParser(description="Génère un rapport HTML à partir d'un fichier TSV agrégé.")
    parser.add_argument("--input_file", required=True, help="Chemin du fichier TSV d'entrée (CombinedVariantOutput.tsv)")
    parser.add_argument("--output_file", required=True, help="Chemin du fichier HTML de sortie")
    parser.add_argument("--metadata_file", required=True, help="Chemin du fichier CSV de métadonnées (patient_metadata.csv)")
    parser.add_argument("--gene_cytoband_file", required=True, help="Chemin du fichier gene_cytoband.txt")
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Erreur : Le fichier '{args.input_file}' n'existe pas.")
        exit(1)
    if not os.path.exists(args.metadata_file):
        print(f"Erreur : Le fichier de métadonnées '{args.metadata_file}' n'existe pas.")
        exit(1)
    if not os.path.exists(args.gene_cytoband_file):
        print(f"Erreur : Le fichier gene_cytoband '{args.gene_cytoband_file}' n'existe pas.")
        exit(1)
    
    gene_cytoband = load_gene_cytoband(args.gene_cytoband_file)
    
    with open(args.input_file, "r", encoding="utf-8") as fin:
        content = fin.read()
    
    sections = parse_variant_file(content)
    analysis_details = parse_analysis_details(sections)
    
    # Récupération des informations d'analyse
    dna_id = analysis_details.get("DNA Sample ID", "NA")
    rna_id = analysis_details.get("RNA Sample ID", "NA")
    sample_name = analysis_details.get("Pair ID", extract_sample_name(args.input_file))
    
    # Mise à jour de la section Copy Number Variants si une valeur de %_tumeur est disponible
    metadata = load_metadata(args.metadata_file)
    tumor_percentage = metadata.get(sample_name, None)
    if tumor_percentage is not None:
        for sec in sections:
            if sec["name"].lower() == "copy number variants":
                sec["lines"] = update_copy_number_lines(sec["lines"], tumor_percentage)
    
    # Construction du rapport HTML en deux blocs (ADN et ARN) + tri custom CNV
    html_result = generate_html_two_blocks(sections, sample_name, dna_id, rna_id, gene_cytoband)
    
    with open(args.output_file, "w", encoding="utf-8") as fout:
        fout.write(html_result)
    
    print(f"Rapport HTML généré : {args.output_file}")

if __name__ == "__main__":
    main()

