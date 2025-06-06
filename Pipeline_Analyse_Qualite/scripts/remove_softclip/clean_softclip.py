#!/usr/bin/env python
import pysam
import argparse
import sys

def remove_soft_clips(aln, header):
    """
    Retire les bases soft-clippées en début et en fin de lecture.
    Le CIGAR est ajusté pour ne conserver que les opérations de matching.
    Les bases soft-clippées (S) sont enlevées de la séquence et des qualités.
    """
    # Si le CIGAR est vide ou qu'il n'y a pas de soft clip, retourner l'alignement inchangé.
    if aln.cigar is None or not any(op == 4 for op, length in aln.cigar):
        return aln

    new_seq = aln.query_sequence
    new_qual = aln.qual
    new_cigar = list(aln.cigar)

    # Suppression du soft clip en début
    if new_cigar[0][0] == 4:
        sc = new_cigar[0][1]
        new_seq = new_seq[sc:]
        new_qual = new_qual[sc:]
        new_cigar = new_cigar[1:]
        # La position de départ (POS) reste inchangée car déjà définie sur le premier aligned base

    # Suppression du soft clip en fin
    if new_cigar and new_cigar[-1][0] == 4:
        sc = new_cigar[-1][1]
        new_seq = new_seq[:-sc]
        new_qual = new_qual[:-sc]
        new_cigar = new_cigar[:-1]

    # Création d'une nouvelle instance d'alignement à partir du dictionnaire et du header
    new_aln = aln.__class__.from_dict(aln.to_dict(), header)
    new_aln.cigar = new_cigar
    new_aln.query_sequence = new_seq
    new_aln.qual = new_qual
    return new_aln

def main():
    parser = argparse.ArgumentParser(description='Nettoyage des soft-clips dans un BAM.')
    parser.add_argument('--input', required=True, help='Fichier BAM en entrée')
    parser.add_argument('--output', required=True, help='Fichier BAM filtré en sortie')
    args = parser.parse_args()

    try:
        in_bam = pysam.AlignmentFile(args.input, "rb")
    except Exception as e:
        sys.stderr.write("Erreur lors de l'ouverture du BAM d'entrée: {}\n".format(e))
        sys.exit(1)

    header = in_bam.header

    try:
        out_bam = pysam.AlignmentFile(args.output, "wb", header=header)
    except Exception as e:
        sys.stderr.write("Erreur lors de la création du BAM de sortie: {}\n".format(e))
        sys.exit(1)

    for aln in in_bam.fetch(until_eof=True):
        try:
            new_aln = remove_soft_clips(aln, header)
            out_bam.write(new_aln)
        except Exception as e:
            sys.stderr.write("Erreur lors du traitement d'un alignement: {}\n".format(e))
    in_bam.close()
    out_bam.close()

if __name__ == '__main__':
    main()

