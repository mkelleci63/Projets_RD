#!/usr/bin/env python3
import csv
import sys

def insert_suffix(patient, suffix):
    if '-' in patient:
        prefix, rest = patient.split('-', 1)
        return f"{prefix}{suffix}-{rest}"
    else:
        return f"{patient}{suffix}"

def main(input_file, output_file):
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile)
        header = next(reader)
        rows = list(reader)
    
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for row in rows:
            # Écriture de la ligne d'origine
            writer.writerow(row)
            patient = row[0]
            # Dupliquer la ligne avec différentes variantes du nom du patient
            row_inserted_d = row.copy()
            row_inserted_d[0] = insert_suffix(patient, 'D')
            writer.writerow(row_inserted_d)
            
            row_inserted_r = row.copy()
            row_inserted_r[0] = insert_suffix(patient, 'R')
            writer.writerow(row_inserted_r)
            
            row_appended_d = row.copy()
            row_appended_d[0] = f"{patient}D"
            writer.writerow(row_appended_d)
            
            row_appended_r = row.copy()
            row_appended_r[0] = f"{patient}R"
            writer.writerow(row_appended_r)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: generate_metadata.py <input_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

