#!/usr/bin/env bash

# Fonction d'affichage de l'usage
usage() {
    echo "Usage:"
    echo "  Pour TSO500/TSO500_HRD : $0 <chemin_run> <type_analyse> <chemin_sample_cancer_file>"
    echo "  Pour ctBLUNG, ctPROST, LUCATISS : $0 <chemin_run> <type_analyse>"
    exit 1
}

# Fonction pour comparer deux versions (v1 >= v2 ?)
compare_versions() {
    local v1="$1"
    local v2="$2"
    IFS='.' read -r v1_major v1_minor v1_patch v1_extra <<< "$v1"
    IFS='.' read -r v2_major v2_minor v2_patch v2_extra <<< "$v2"
    v1_major=${v1_major:-0}
    v1_minor=${v1_minor:-0}
    v1_patch=${v1_patch:-0}
    v2_major=${v2_major:-0}
    v2_minor=${v2_minor:-0}
    v2_patch=${v2_patch:-0}
    if (( v1_major > v2_major )); then
        return 0
    elif (( v1_major < v2_major )); then
        return 1
    else
        if (( v1_minor > v2_minor )); then
            return 0
        elif (( v1_minor < v2_minor )); then
            return 1
        else
            if (( v1_patch >= v2_patch )); then
                return 0
            else
                return 1
            fi
        fi
    fi
}

# 1. Récupère le chemin absolu du dossier contenant ce script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 2. Vérifie le nombre d'arguments en fonction du type d'analyse
if [ "$#" -lt 2 ]; then
    usage
fi

run_dir="$1"
type_analyse="$2"

if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    if [ "$#" -ne 3 ]; then
        usage
    fi
    input_metadata_file="$3"
else
    input_metadata_file=""
fi

# Création du dossier QC_Analyse et du sous-dossier erreurs dans le dossier run
output_dir="${run_dir}/QC_Analyse"
erreurs_dir="${output_dir}/Execution_pipeline/Erreurs"
mkdir -p "${erreurs_dir}"

# Le fichier d'erreur sera créé dans ce dossier
error_log_file="${erreurs_dir}/erreurs_données.txt"
> "$error_log_file"

#############################################
# Normalisation du fichier de métadonnées (pour TSO500/TSO500_HRD uniquement)
#############################################
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    normalized_file="${SCRIPT_DIR}/metadata_samples.csv"
    {
      echo "Patient,Cancer,%_tumeur"
      line_num=2
      tail -n +2 "$input_metadata_file" | while IFS= read -r line; do
          if [ -z "$line" ]; then
              continue
          fi
          nf=$(echo "$line" | awk -F',' '{print NF}')
          if [ "$nf" -eq 3 ]; then
              echo "$line"
          elif [ "$nf" -eq 2 ]; then
              patient=$(echo "$line" | cut -d',' -f1 | xargs)
              field2=$(echo "$line" | cut -d',' -f2 | xargs)
              if [[ $field2 =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                  echo "$patient,,$field2"
              else
                  echo "$patient,$field2,"
              fi
          elif [ "$nf" -eq 1 ]; then
              patient=$(echo "$line" | xargs)
              echo "$patient,,"
          elif [ "$nf" -gt 3 ]; then
              patient=$(echo "$line" | cut -d',' -f1 | xargs)
              extra_fields=$(echo "$line" | cut -d',' -f4-)
              echo "[Error] Ligne $line_num: Pour l'échantillon '$patient', nombre de champs supérieur à 3 ($nf champs détectés). Champs en trop: $extra_fields" >> "$error_log_file"
              echo "$line"
          else
              echo "[Error] Ligne $line_num: nombre de champs inattendu ($nf)." >> "$error_log_file"
              exit 1
          fi
          line_num=$((line_num+1))
      done
    } > "$normalized_file"
fi

#############################################
# Vérification du format du fichier de métadonnées (pour TSO500/TSO500_HRD uniquement)
#############################################
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    declare -a metadata_format_errors
    declare -a metadata_format_warnings

    if [ ! -f "$input_metadata_file" ]; then
        metadata_format_errors+=("[Error] Le fichier de métadonnées fourni n'existe pas : ${input_metadata_file}")
    else
        header_line=$(head -n 1 "$input_metadata_file" | tr -d '\r\n')
        expected_header="Patient,Cancer,%_tumeur"
        if [ "$header_line" != "$expected_header" ]; then
            metadata_format_errors+=("[Error] Format du fichier de métadonnées incorrect. Attendu : '$expected_header'. Le fichier fourni contient : '$header_line'")
        fi
    fi

    if [ -f "$normalized_file" ]; then
        line_num=2
        while IFS= read -r line; do
            if [ "$line" = "Patient,Cancer,%_tumeur" ]; then
                continue
            fi
            nf=$(echo "$line" | awk -F',' '{print NF}')
            if [ "$nf" -gt 3 ]; then
                patient=$(echo "$line" | cut -d',' -f1 | xargs)
                extra_fields=$(echo "$line" | cut -d',' -f4-)
                metadata_format_errors+=("[Error] Ligne $line_num: Pour l'échantillon '$patient', nombre de champs supérieur à 3 ($nf champs détectés). Champs en trop: $extra_fields")
                line_num=$((line_num + 1))
                continue
            fi
            patient=$(echo "$line" | cut -d',' -f1 | xargs)
            if [ "$nf" -eq 3 ]; then
                cancer=$(echo "$line" | cut -d',' -f2 | xargs)
                tumor=$(echo "$line" | cut -d',' -f3 | xargs)
            elif [ "$nf" -eq 2 ]; then
                second_field=$(echo "$line" | cut -d',' -f2 | xargs)
                if [[ $second_field =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    cancer=""
                    tumor="$second_field"
                else
                    cancer="$second_field"
                    tumor=""
                fi
            elif [ "$nf" -eq 1 ]; then
                cancer=""
                tumor=""
            fi
            if [ -z "$cancer" ]; then
                metadata_format_warnings+=("[Warning] Ligne $line_num: Pour l'échantillon '$patient', le type de cancer n'est pas précisé.")
            fi
            if [ -z "$tumor" ]; then
                metadata_format_warnings+=("[Warning] Ligne $line_num: Pour l'échantillon '$patient', le % de cellules tumorales n'est pas précisé.")
            fi
            line_num=$((line_num + 1))
        done < "$normalized_file"
    fi

    for err in "${metadata_format_errors[@]}"; do
        echo "$err" >> "$error_log_file"
        echo "$err"
    done
    for warn in "${metadata_format_warnings[@]}"; do
        echo "$warn" >> "$error_log_file"
        echo "$warn"
    done
    if [ ${#metadata_format_errors[@]} -gt 0 ]; then
        exit 1
    fi
fi

#############################################
# Pour TSO500/TSO500_HRD : génération des métadonnées étendues
#############################################
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    metadata_dir="${SCRIPT_DIR}/metadata"
    output_metadata_file="${metadata_dir}/metadata_extended.csv"
    python "${SCRIPT_DIR}/scripts/metadata/generate_metadata.py" "${normalized_file}" "${output_metadata_file}"
    sample_cancer_file="${output_metadata_file}"
else
    sample_cancer_file=""
fi

#############################################
# Section 1 : Vérification du type d'analyse
#############################################
declare -a analysis_errors
allowed_analyses=("ctBLUNG" "ctPROST" "LUCATISS" "TSO500" "TSO500_HRD")
type_valid=0
for allowed in "${allowed_analyses[@]}"; do
    if [ "$type_analyse" == "$allowed" ]; then
        type_valid=1
        break
    fi
done
if [ $type_valid -eq 0 ]; then
    analysis_errors+=("[Error] Le type d'analyse '$type_analyse' est inconnu. Les analyses possibles sont : ${allowed_analyses[*]}.")
fi

#############################################
# Section 2 : Vérification du contenu du fichier de métadonnées (pour TSO500/TSO500_HRD)
#############################################
declare -a metadata_errors
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    valid_cancers=("Breast" "CNS" "Colorectal" "Lung" "Melanoma" "Ovarian" "Prostate" "Thyroid" "Uterine and Cervical")
    valid_cancers_str="Breast, CNS, Colorectal, Lung, Melanoma, Ovarian, Prostate, Thyroid, Uterine and Cervical"
    if [ ! -f "$sample_cancer_file" ]; then
        metadata_errors+=("[Error] Le fichier de métadonnées étendu n'existe pas : ${sample_cancer_file}")
    fi
    declare -A sample_errors
    add_error() {
        local sample="$1"
        local error_msg="$2"
        if [ -z "${sample_errors[$sample]}" ]; then
            sample_errors[$sample]="    - $error_msg"
        else
            sample_errors[$sample]="${sample_errors[$sample]}"$'\n'"    - $error_msg"
        fi
    }
    if [ -f "$normalized_file" ]; then
        while IFS= read -r line; do
            if [ "$line" = "Patient,Cancer,%_tumeur" ]; then continue; fi
            nf=$(echo "$line" | awk -F',' '{print NF}')
            if [ "$nf" -gt 3 ]; then
                patient=$(echo "$line" | cut -d',' -f1 | xargs)
                extra_fields=$(echo "$line" | cut -d',' -f4-)
                add_error "$patient" "Nombre de champs supérieur à 3 ($nf champs détectés). Champs en trop: $extra_fields"
                continue
            fi
            patient=$(echo "$line" | cut -d',' -f1 | xargs)
            if [ "$nf" -eq 3 ]; then
                cancer=$(echo "$line" | cut -d',' -f2 | xargs)
                percent=$(echo "$line" | cut -d',' -f3 | xargs)
            elif [ "$nf" -eq 2 ]; then
                second_field=$(echo "$line" | cut -d',' -f2 | xargs)
                if [[ $second_field =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    cancer=""
                    percent="$second_field"
                else
                    cancer="$second_field"
                    percent=""
                fi
            else
                cancer=""
                percent=""
            fi
            if [ -n "$cancer" ]; then
                valid=0
                for allowed in "${valid_cancers[@]}"; do
                    if [ "$cancer" == "$allowed" ]; then valid=1; break; fi
                done
                if [ $valid -eq 0 ]; then
                    add_error "$patient" "Le type de cancer '$cancer' n'est pas reconnu. Les valeurs acceptées sont : $valid_cancers_str."
                fi
            fi
            if [ -n "$percent" ]; then
                if ! [[ $percent =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    add_error "$patient" "Le pourcentage '$percent' n'est pas un nombre valide."
                else
                    if ! [[ $percent =~ ^[0-9]+$ ]]; then
                        add_error "$patient" "Le pourcentage '$percent' n'est pas un entier."
                    fi
                    if (( $(echo "$percent < 0" | bc -l) )) || (( $(echo "$percent > 100" | bc -l) )); then
                        add_error "$patient" "Le pourcentage '$percent' doit être compris entre 0 et 100."
                    fi
                fi
            fi
        done < "$normalized_file"
    fi
    for sample in "${!sample_errors[@]}"; do
        metadata_errors+=("Échantillon : $sample"$'\n'"${sample_errors[$sample]}")
    done
fi

#############################################
# Section 3 : Vérification de la structure du run
#############################################
declare -a structure_errors
declare -a structure_warnings
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    if [ ! -d "$run_dir" ]; then
        structure_errors+=("[Error] Le dossier du run n'existe pas : ${run_dir}")
    fi
    results_dir="${run_dir}/Results"
    if [ ! -d "$results_dir" ]; then
        structure_errors+=("[Error] Dans le run, le dossier 'Results' est absent.")
    fi
    logs_dir="${run_dir}/Logs_Intermediates"
    if [ ! -d "$logs_dir" ]; then
        structure_errors+=("[Error] Dans le run, le dossier 'Logs_Intermediates' est absent.")
    fi

    metrics_file="${results_dir}/MetricsOutput.tsv"
    if [ ! -f "$metrics_file" ]; then
        structure_errors+=("[Error] Pour l'analyse QC TSO500, le fichier 'MetricsOutput.tsv' est absent ou mal nommé dans ${results_dir}.")
    else
        workflow_line=$(grep -m1 "^Workflow Version" "$metrics_file")
        if [ -z "$workflow_line" ]; then
            structure_errors+=("[Error] Impossible de trouver la ligne 'Workflow Version' dans $metrics_file.")
        else
            workflow_version=$(echo "$workflow_line" | awk -F'\t' '{print $2}')
            compare_versions "$workflow_version" "2.6.0"
            if [ $? -eq 1 ]; then
                structure_errors+=("[Error] Attention, la version ($workflow_version) est antérieure à 2.6.0. L'analyse QC n'est pas disponible pour ce run.")
            fi
        fi
    fi

    for sample_dir in "$results_dir"/*; do
        if [ -d "$sample_dir" ]; then
            sample_name=$(basename "$sample_dir")
            if [ ! -f "$sample_dir/${sample_name}_CombinedVariantOutput.tsv" ]; then
                structure_errors+=("[Error] Le fichier CombinedVariantOutput.tsv est absent ou mal nommé dans le dossier '$sample_name' de Results.")
            fi
        fi
    done

elif [[ "$type_analyse" == "ctBLUNG" || "$type_analyse" == "ctPROST" || "$type_analyse" == "LUCATISS" ]]; then
    # Pour ces analyses, on collecte les avertissements sans stopper l'analyse.
    if [[ "$type_analyse" == "ctBLUNG" || "$type_analyse" == "ctPROST" ]]; then
        for sample_path in "$run_dir"/*; do
            if [ -d "$sample_path" ]; then
                sample_name=$(basename "$sample_path")
                if [[ "$sample_name" =~ [Ff][Aa][Ss][Tt][Qq] ]]; then continue; fi
                if [ ! -f "$sample_path/${sample_name}.mapping_metrics.csv" ]; then
                    structure_warnings+=("[Warning] Dans le dossier '$sample_path', le fichier '${sample_name}.mapping_metrics.csv' est manquant ou mal nommé. Vérifiez la présence du fichier nécessaire pour le calcul des taux de duplicats. Si ce n'est pas un dossier échantillon, ne prenez pas en compte ce message.")
                fi
            fi
        done
    elif [[ "$type_analyse" == "LUCATISS" ]]; then
        adn_dir="$run_dir/ADN"
        arn_dir="$run_dir/ARN"
        if [ ! -d "$adn_dir" ]; then
             structure_warnings+=("[Warning] Le dossier 'ADN' est absent dans le run pour l'analyse LUCATISS. Vérifiez le chemin : ${run_dir}/ADN. Si ce n'est pas un dossier échantillon, ne prenez pas en compte ce message.")
        else
             for sample_path in "$adn_dir"/*; do
                  if [ -d "$sample_path" ]; then
                      sample_name=$(basename "$sample_path")
                      if [ ! -f "$sample_path/${sample_name}.umi_metrics.csv" ]; then
                          structure_warnings+=("[Warning] Dans le dossier '$sample_path', le fichier '${sample_name}.umi_metrics.csv' est manquant ou mal nommé. Vérifiez la présence du fichier nécessaire pour le calcul des taux de duplicats. Si ce n'est pas un dossier échantillon, ne prenez pas en compte ce message.")
                      fi
                  fi
             done
        fi
        if [ ! -d "$arn_dir" ]; then
             structure_warnings+=("[Warning] Le dossier 'ARN' est absent dans le run pour l'analyse LUCATISS. Vérifiez le chemin : ${run_dir}/ARN. Si ce n'est pas un dossier échantillon, ne prenez pas en compte ce message.")
        else
             for sample_path in "$arn_dir"/*; do
                  if [ -d "$sample_path" ]; then
                      sample_name=$(basename "$sample_path")
                      if [ ! -f "$sample_path/${sample_name}.mapping_metrics.csv" ]; then
                          structure_warnings+=("[Warning] Dans le dossier '$sample_path', le fichier '${sample_name}.mapping_metrics.csv' est manquant ou mal nommé. Vérifiez la présence du fichier nécessaire pour le calcul des taux de duplicats. Si ce n'est pas un dossier échantillon, ne prenez pas en compte ce message.")
                      fi
                  fi
             done
        fi
    fi
fi

#############################################
# Section 4 : Vérification des noms d'échantillons dans Results (TSO500/TSO500_HRD uniquement)
#############################################
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    detected_samples=()
    if [ -d "$results_dir" ]; then
        for sample_path in "$results_dir"/*; do
            if [ -d "$sample_path" ]; then
                detected_samples+=("$(basename "$sample_path")")
            fi
        done
    fi

    missing_samples=()
    if [ -d "$results_dir" ]; then
        for sample_path in "$results_dir"/*; do
            if [ ! -d "$sample_path" ]; then continue; fi
            sample_name=$(basename "$sample_path")
            if ! grep -q "^${sample_name}" "$input_metadata_file"; then
                missing_samples+=("$sample_name")
            fi
        done
        if [ ${#missing_samples[@]} -ne 0 ]; then
           metadata_errors+=("[Error] Les échantillons suivants sont présents dans le run mais absents ou erronés dans le fichier de métadonnées : ${missing_samples[*]}.
Les échantillons détectés dans le run sont : ${detected_samples[*]}. Veuillez vérifier la saisie exacte des informations.")
        fi
    fi
fi

#############################################
# Affichage et écriture des erreurs et avertissements
#############################################
print_section() {
    local title="$1"
    shift
    local messages=("$@")
    if [ ${#messages[@]} -ne 0 ]; then
        echo "========================================"
        echo "$title"
        echo "========================================"
        for msg in "${messages[@]}"; do
            echo "$msg"
        done
        echo ""
    fi
}

global_messages=""
if [ ${#analysis_errors[@]} -ne 0 ]; then
    global_messages+="========================================\nErreurs sur le type d'analyse:\n========================================\n"
    for err in "${analysis_errors[@]}"; do
        global_messages+="$err\n"
    done
    global_messages+="\n"
fi
if [ ${#metadata_errors[@]} -ne 0 ]; then
    global_messages+="========================================\nErreurs dans le fichier de métadonnées:\n========================================\n"
    for err in "${metadata_errors[@]}"; do
        global_messages+="$err\n----------------------------------------\n"
    done
    global_messages+="\n"
fi
if [ ${#structure_errors[@]} -ne 0 ]; then
    global_messages+="========================================\nErreurs dans la structure du run:\n========================================\n"
    for err in "${structure_errors[@]}"; do
        global_messages+="$err\n"
    done
    global_messages+="\n"
fi
if [ ${#structure_warnings[@]} -ne 0 ]; then
    global_messages+="========================================\nAvertissements sur la structure du run:\n========================================\n"
    for warn in "${structure_warnings[@]}"; do
        global_messages+="$warn\n"
    done
    global_messages+="\n"
fi
if [ -n "$global_messages" ]; then
    echo -e "$global_messages" >> "$error_log_file"
    echo -e "$global_messages"
    # Pour TSO500/TSO500_HRD, une erreur critique arrête l'analyse
    if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
        exit 1
    fi
fi

echo "[Info] Toutes les vérifications ont été passées avec succès."

#############################################
# Suite du script : lancement de Nextflow
#############################################
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    input="${results_dir}/MetricsOutput.tsv"
    echo "[Info] Analyse ${type_analyse} détectée."
    echo "[Info] Fichier d'entrée: ${input}"
else
    input="${run_dir}"
    echo "[Info] Analyse ${type_analyse} détectée."
    echo "[Info] Répertoire d'entrée: ${input}"
fi

ENV_FILE="${SCRIPT_DIR}/environnement.yml"
if [ ! -f "${ENV_FILE}" ]; then
    echo "[Error] Impossible de trouver le fichier environnement.yml à l'emplacement attendu : ${ENV_FILE}" >> "$error_log_file"
    exit 1
fi

echo "[Info] Lancement de Nextflow depuis : ${SCRIPT_DIR}/main_test.nf"
echo "[Info] Environnement conda : ${ENV_FILE}"
if [[ "$type_analyse" == "TSO500" || "$type_analyse" == "TSO500_HRD" ]]; then
    echo "[Info] Fichier de métadonnées (sample_cancer_file) : ${sample_cancer_file}"
fi

nextflow run "${SCRIPT_DIR}/main_test.nf" \
    --input_file "${input}" \
    --analysis_type "${type_analyse}" \
    --sample_cancer_file "${sample_cancer_file}" \
    --env_file "${ENV_FILE}" \
    --with-dag workflow.png

