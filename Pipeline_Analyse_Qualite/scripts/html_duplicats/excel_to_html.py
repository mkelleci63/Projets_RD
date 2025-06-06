import argparse
import io
import re

from bs4 import BeautifulSoup
from xlsx2html import xlsx2html
import openpyxl

def remove_colgroups(html):
    """
    Supprime toutes les balises <colgroup> et <col> du fragment HTML.
    """
    soup = BeautifulSoup(html, "html.parser")
    for colgroup in soup.find_all("colgroup"):
        colgroup.decompose()
    for col in soup.find_all("col"):
        col.decompose()
    return str(soup)

def is_cell_empty(cell):
    """
    Détermine si une cellule est vide en vérifiant son contenu après avoir supprimé
    les balises et entités HTML.
    """
    content_html = cell.decode_contents(formatter="html").strip()
    # Supprime toutes les balises
    text_no_tags = re.sub(r'<[^>]*>', '', content_html)
    # Supprime les entités HTML (ex. &nbsp;)
    text_no_entities = re.sub(r'&[^;]+;', '', text_no_tags)
    return not text_no_entities.strip()

def fix_scientific_notation(html):
    """
    Convertit les nombres en notation scientifique (ex. 1.23e+05) en notation classique.
    """
    soup = BeautifulSoup(html, "html.parser")
    pattern = re.compile(r'^-?\d(\.\d+)?[eE][+-]?\d+$')
    for cell in soup.find_all(['td', 'th']):
        text = cell.get_text(strip=True)
        if pattern.match(text):
            try:
                value = float(text)
                # On formate en notation classique
                new_text = format(value, 'f')
                cell.string.replace_with(new_text)
            except ValueError:
                pass
    return str(soup)

def remove_empty_rows(html):
    """
    Supprime les lignes <tr> qui sont entièrement vides.
    """
    soup = BeautifulSoup(html, "html.parser")
    for tr in soup.find_all("tr"):
        cells = tr.find_all(['td', 'th'])
        if all(is_cell_empty(cell) for cell in cells):
            tr.decompose()
    return str(soup)

def convert_excel_to_html(input_file, output_file):
    # Charge le classeur pour obtenir la liste des feuilles
    wb = openpyxl.load_workbook(input_file, read_only=True)
    
    html_parts = []
    # En-tête du document HTML
    html_parts.append("<html>")
    html_parts.append("<head>")
    html_parts.append("<meta charset='UTF-8'>")
    html_parts.append("<title>Détails duplicats</title>")
    html_parts.append("<style>")
    html_parts.append("  table { border-collapse: collapse; margin-bottom: 20px; }")
    html_parts.append("  table, th, td { border: 1px solid #000; }")
    # On force l'affichage sur une seule ligne ET on augmente le padding
    html_parts.append("  td, th {")
    html_parts.append("      font-size: 18px !important;")
    html_parts.append("      white-space: nowrap !important;")
    html_parts.append("      padding: 10px 15px !important;")  # plus d'espace
    # Optionnel : impose une largeur minimale
    # html_parts.append("      min-width: 100px;")  # A décommenter si on veut un min-width
    html_parts.append("  }")
    html_parts.append("  h2 { margin-top: 30px; }")
    html_parts.append("</style>")
    html_parts.append("</head>")
    html_parts.append("<body>")
    html_parts.append("<h1>Détails duplicats</h1>")
    
    # Parcours des feuilles
    for sheet_name in wb.sheetnames:
        html_sheet = xlsx2html(input_file, sheet=sheet_name)
        
        # Si xlsx2html renvoie un StringIO, on le convertit en chaîne
        if isinstance(html_sheet, io.StringIO):
            html_sheet = html_sheet.getvalue()
        
        # 1) Supprime <colgroup> et <col>
        html_sheet = remove_colgroups(html_sheet)
        # 2) Corrige la notation scientifique
        html_sheet = fix_scientific_notation(html_sheet)
        # 3) Supprime les lignes entièrement vides
        html_sheet = remove_empty_rows(html_sheet)
        
        # Ajoute un titre pour la feuille
        html_parts.append(f"<h2>{sheet_name}</h2>")
        html_parts.append(html_sheet)
    
    html_parts.append("</body>")
    html_parts.append("</html>")
    
    # Écriture du fichier final
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(html_parts))
    
    print(f"Conversion terminée : {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convertir un fichier Excel en HTML en conservant les styles, "
                    "mais en supprimant les <colgroup>, forçant l'affichage sur une seule ligne, etc."
    )
    parser.add_argument("input", help="Chemin vers le fichier Excel (.xlsx)")
    parser.add_argument("output", help="Chemin vers le fichier HTML généré")
    args = parser.parse_args()
    
    convert_excel_to_html(args.input, args.output)

