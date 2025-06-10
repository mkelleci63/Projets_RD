# Projet analyse QC

- Ce dépot contient un pipeline permettant d'effectuer une analyse QC pour des runs passant des analyses ctBLUNG,ctPROST,TSO500/TSO500_HRD
- L'accent est mis sur la détection des régions en low_coverage 50X et sur le calcul des taux de duplicats des échantillons.

# Contenu du répertoire QC_ANALYSE : 
- run_pipeline.sh : Ce script doit être exécuté en premier, il permet de prendre en arguments les informations sur le run ( année,analyse,chemin )
- main.nf : Fichier principal listant tous les process du workflow
- environnement.yml : Fichier contenant toutes les dépendances
- nextflow.config : Fichier indiquant l'emplacement de création de l'environnement conda et quels process utilisent cet environnement
- metadata : Dossier contenant l'ensemble des fichiers de metadonnées nécessaires aux analyses effectuées par le pipeline
- scripts : Dossier contenant l'ensemble des scripts appellés et utilisés dans le pipeline
    
# Lancement du pipeline : 
- run_pipeline.sh <chemin_run> <type_analyse> <chemin_fichier_métadonnées>
- exemple : run_pipeline.sh /mnt/dragen_data/2024/TSO500/241120_NDX550301_RUO_0365_AHLKNYBGXW_RD-2.6.0 TSO500 /mnt/dev/mkelleci/patient_metadata.csv

# Output : 
- Un dossier QC_Analyse sera disponible dans le dossier de sortie du run.
- Dans ce dossier se trouvera un rapport html contenant l'ensemble des analyses QC du run

# Prérequis :
- Nextflow doit être installé dans l'environnement de travail
- Le pipeline doit être executé dans l'environnement de base du terminal et non à l'intérieur d'un environnement pré-existant
- Pour les analyses TSO500/TSO500_HRD, il faut bien faire attention à ce que le fichier /metadata/patient_metadata.csv contienne les bonnes informations :
    - Pour ces deux analyses il faut préciser le nom de tous les échantillons du run ( ADN et/ou ARN )
    - Le pipeline est fonctionnel pour les runs ayant passé une analyse TSO500 avec la version 2.6.0 au minimum
    
# Remarques :
- Le fichier de métadonnées patient_metadata.csv contient par défaut les échantillons du run /mnt/dragen_data/2024/TSO500/241120_NDX550301_RUO_0365_AHLKNYBGXW_RD-2.6.0
- Le fichier de métadonnées doit être au format csv avec la présence du header suivant : Patient,Cancer,%_tumeur
