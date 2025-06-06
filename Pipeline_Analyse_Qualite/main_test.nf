#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//////////////////////////////////////////////////////////////////////////
// Paramètres globaux et chemins vers les scripts/métadonnées
//////////////////////////////////////////////////////////////////////////

params.input_file               = null    // Pour TSO500 : chemin vers MetricsOutput.tsv, sinon chemin du run
params.analysis_type            = null    // Valeurs attendues : TSO500, TSO500_HRD, ctBLUNG, ctPROST, LUCATISS

// Scripts communs
params.consolidate_script       = "${projectDir}/scripts/seuils_ref_metrics/consolidate_metrics.py"
params.compare_script           = "${projectDir}/scripts/compare_metrics/compare_html.py"

// Traitement DNA et ARN
params.exon_report_script       = "${projectDir}/scripts/filtre_50X/ADN/process_exon_reports.py"
params.filter_script            = "${projectDir}/scripts/filtre_cancer/ADN/script_filtrage.py"
params.rna_exon_report_script   = "${projectDir}/scripts/filtre_50X/ARN/process_rna_exon_reports.py"
params.rna_filter_script        = "${projectDir}/scripts/filtre_cancer/ARN/script_filtrage_rna.py"

// Métadonnées
params.sample_cancer_file       = params.sample_cancer_file
params.cancer_gene_file         = "${projectDir}/metadata/cancer_gene_dict.csv"
params.alias_cancer_file        = "${projectDir}/metadata/filtered_alias_lookup.json"


// Rapport HTML DNA
params.generate_html_script     = "${projectDir}/scripts/html_variant/generate_html_report.py"

// Nouveau paramètre pour le fichier gene_cytoband
params.gene_cytoband_file       = "${projectDir}/metadata/gene_cytoband.txt"

// Références pour le coverage
params.fasta_file               = "${projectDir}/metadata/GCF_000001405.25_GRCh37.p13_genomic.fna"
params.gtf_file                 = "${projectDir}/metadata/GRCh37_latest_genomic.gtf"

// Pour le workflow TSO500
params.duplicats_script         = "${projectDir}/scripts/html_duplicats/TSO500/duplicats_TSO500.py"
params.excel_to_html_script     = "${projectDir}/scripts/html_duplicats/excel_to_html.py"

// Pour le workflow secondaire (ctBLUNG, ctPROST, LUCATISS)
params.duplicats_lucatiss_script = "${projectDir}/scripts/html_duplicats/LUCATISS/duplicats_LUCATISS.py"
params.duplicats_ct_script       = "${projectDir}/scripts/html_duplicats/ctBLUNG_ctPROST/duplicats_ctBLUNG_ctPROST.py"

// Pour le script combinant les fichiers HTML
params.combine_html_script       = "${projectDir}/scripts/rapport_final/combine_html.py"

// Fixation d'un seuil à 50X pour supprimer les fichiers samtools depth inutiles
params.coverage_threshold = params.coverage_threshold ?: 50


//////////////////////////////////////////////////////////////////////////
// Workflow principal
//////////////////////////////////////////////////////////////////////////

workflow {

  if( !params.analysis_type ) {
      error "Vous devez spécifier le type d'analyse (--analysis_type). Exemples : TSO500, TSO500_HRD, ctBLUNG, ctPROST, LUCATISS"
  }

  if( params.analysis_type in ['TSO500','TSO500_HRD'] ) {
      println "Exécution du workflow TSO500 (${params.analysis_type})"

      // Dans cette branche, params.input_file = chemin vers MetricsOutput.tsv.
      def chosen_file = params.input_file
      if( !chosen_file ) {
          log.info "Aucun fichier input spécifié. Ouverture de la boîte de dialogue..."
          def cmd = "python ${params.file_chooser_script}".execute()
          cmd.waitFor()
          def inputPath = cmd.in.text.trim()
          if( !inputPath || inputPath == 'null' ) {
              error "Aucun fichier sélectionné : arrêt du pipeline."
          }
          chosen_file = inputPath
          println "Vous avez sélectionné : ${chosen_file}"
      }
      
      // Déduction du chemin du run
      def resultsDir = new File(chosen_file).parent
      def parentDir  = new File(resultsDir).parent
      def run_name = new File(parentDir).name
      println "Nom du run : ${run_name}"

      // Dossiers de sortie
      params.output_dir            = "${parentDir}/QC_Analyse"
      params.combined_dir          = "${params.output_dir}/Rapport_Final"
      params.htm_dupli_dir         = "${params.taux_duplicats_dir}/Html"
      [ params.combined_dir ].each { new File(it).mkdirs() }
      println "Dossiers de sortie créés sous : ${params.output_dir}"

      // Recherche du dossier 'dragen_data'
      def currentDir = new File(chosen_file)
      while( currentDir && currentDir.getName() != "dragen_data" ) {
          currentDir = currentDir.getParentFile()
      }
      if( !currentDir ) {
          error "Le dossier 'dragen_data' n'a pas été trouvé dans l'arborescence de ${chosen_file}"
      }
      def parent_directory = currentDir.toString()
      println "Dossier 'dragen_data' trouvé : ${parent_directory}"

      // Processus intermédiaires
      def metricsSummary = consolidateMetrics(parent_directory)

      // 1) compareMetrics => produit MetricsComparison.html
      def compareOut = compareMetrics(metricsSummary, file(chosen_file), file(params.compare_script))

      def logsDir = new File(parentDir, "Logs_Intermediates")
      def logsChannel = Channel.value(logsDir.toString())
      
      def exonReportsPaths     = listExonReports(logsChannel)
      def exonReadReportsPaths = listExonReadReports(logsChannel)
      def targetReportsPaths   = listTargetReports(logsChannel)
      def dnaExonFiles         = processExonReports(exonReportsPaths, exonReadReportsPaths, targetReportsPaths)
      def mergedDnaExonFiles   = dnaExonFiles.mix().collect()
      
      def exonReportsPathsRna     = listExonReportsRna(logsChannel)
      def exonReadReportsPathsRna = listExonReadReportsRna(logsChannel)
      def rnaExonFiles            = processExonReportsRna(exonReportsPathsRna, exonReadReportsPathsRna)
      def mergedRnaExonFiles      = rnaExonFiles.collect()
      
      def filteredDNA_files = filterCancerGenesDNA(mergedDnaExonFiles)
      def filteredDNATuples = filteredDNA_files.flatten().map { file ->
           def sample = file.getName().tokenize('.')[0]
           tuple(sample, file)
      }
      filteredDNATuples.subscribe { println "Filtered DNA tuple: $it" }
      
      def filteredRNA_files = filterCancerGenesRNA(mergedRnaExonFiles)
      def filteredRNATuples = filteredRNA_files.flatten().map { file ->
           def sample = file.getName().tokenize('.')[0]
           tuple(sample, file)
      }
      filteredRNATuples.subscribe { println "Filtered RNA tuple: $it" }
      
      def combinedVariantList = listCombinedVariantOutputs(resultsDir)
      combinedVariantList.flatMap { it.readLines() }.set { tsvPathsDNA }

      // 2) generateHtmlReport => produit *.html (variants)
      def variantsHtml = generateHtmlReport(tsvPathsDNA)

      // --- Traitement des BAM sans étape de copie explicite ---
      // On récupère les chemins des fichiers BAM/BAI directement depuis les fichiers texte générés.
      def bamListFileDNA = listBamFilesDNA(logsChannel)
      def bamListFileRNA = listBamFilesRNA(logsChannel)
      
      // Conversion des lignes en objets fichier pour DNA
      def dnaBamFiles = bamListFileDNA.flatMap { it.readLines() }
                                      .map { line -> file(line) }
      def dnaBamTuples = dnaBamFiles.map { f ->
           if( f.name.endsWith('.bam') || f.name.endsWith('.bam.bai') ) {
               def sample = f.name.split('\\.')[0]
               tuple(sample, f)
           }
      }.filter { it != null }
      dnaBamTuples.subscribe { println "DNA tuple: $it" }
      
      // Conversion pour RNA
      def rnaBamFiles = bamListFileRNA.flatMap { it.readLines() }
                                      .map { line -> file(line) }
      def rnaBamTuples = rnaBamFiles.map { f ->
           if( f.name.endsWith('.bam') || f.name.endsWith('.bam.bai') ) {
               def sample = f.name.split('\\.')[0]
               tuple(sample, f)
           }
      }.filter { it != null }
      rnaBamTuples.subscribe { println "RNA tuple: $it" }
      
      // Jointure des tuples (BAM, fichier région) pour DNA
      def joinedDNABamRegion = dnaBamTuples.join(filteredDNATuples)
      joinedDNABamRegion.subscribe { println "Join DNA (BAM, region): $it" }
      def reducedDNATuplesInput = joinedDNABamRegion.map { sample, bam, region -> tuple(sample, bam, region) }
      def reducedDNA = reduceBamRegions_DNA(reducedDNATuplesInput)
      def reducedDNATupleFinal = reducedDNA.map { sample, reducedBam, reducedBai -> tuple(sample, reducedBam, reducedBai) }
      reducedDNATupleFinal.subscribe { println "Reduced DNA tuple: $it" }
      
      // Jointure pour RNA
      def joinedRNABamRegion = rnaBamTuples.join(filteredRNATuples)
      joinedRNABamRegion.subscribe { println "Join RNA (BAM, region): $it" }
      def reducedRNATuplesInput = joinedRNABamRegion.map { sample, bam, region -> tuple(sample, bam, region) }
      def reducedRNA = reduceBamRegions_RNA(reducedRNATuplesInput)
      def reducedRNATupleFinal = reducedRNA.map { sample, reducedBam, reducedBai -> tuple(sample, reducedBam, reducedBai) }
      reducedRNATupleFinal.subscribe { println "Reduced RNA tuple: $it" }
      
      // --- Ajout de la jointure pour le calcul de la profondeur ---
      def finalJoinDNATuples = reducedDNATupleFinal.join(filteredDNATuples)
      finalJoinDNATuples.subscribe { println "Final Join (Reduced BAM, filtered TSV) - DNA: $it" }
      
      def finalJoinRNATuples = reducedRNATupleFinal.join(filteredRNATuples)
      finalJoinRNATuples.subscribe { println "Final Join (Reduced BAM, filtered TSV) - RNA: $it" }
      
      // Création des canaux contenant les régions pour le calcul de la profondeur
      def regionDepthDNA = finalJoinDNATuples.flatMap { sample, bam, bai, tsvFile ->
          def regions = []
          tsvFile.eachLine { line, index ->
              if( index == 0 || line.toLowerCase().startsWith("chrom") || line.startsWith("#") )
                  return
              def parts    = line.split('\t')
              def chrom    = parts[0].startsWith("chr") ? parts[0] : "chr" + parts[0]
              def rawStart = parts[1] as Integer
              def rawEnd   = parts[2] as Integer
              def flank    = 50
              // élargir la région de flank bp de chaque côté
              def start    = Math.max(1, rawStart - flank)
              def end      = rawEnd + flank
              def rawName  = parts[3]
              def safeName = rawName.replaceAll('[^A-Za-z0-9._-]+','_').take(50)
              regions.add( tuple(sample, chrom, start.toString(), end.toString(), safeName, bam, bai) )
          }
          return regions
      }
      
      def regionDepthRNA = finalJoinRNATuples.flatMap { sample, bam, bai, tsvFile ->
          def regions = []
          tsvFile.eachLine { line, index ->
              if( index == 0 || line.toLowerCase().startsWith("chrom") || line.startsWith("#") )
                  return
              def parts    = line.split('\t')
              def chrom    = parts[0].startsWith("chr") ? parts[0] : "chr" + parts[0]
              def rawStart = parts[1] as Integer
              def rawEnd   = parts[2] as Integer
              def flank    = 50
              // élargir la région de flank bp de chaque côté
              def start    = Math.max(1, rawStart - flank)
              def end      = rawEnd + flank
              def rawName  = parts[3]
              def safeName = rawName.replaceAll('[^A-Za-z0-9._-]+','_').take(50)
              regions.add( tuple(sample, chrom, start.toString(), end.toString(), safeName, bam, bai) )
          }
          return regions
      }
        
      // Processus de calcul de la profondeur via samtools depth
      def depthDNATuples = calculateDepthDNA(regionDepthDNA)
      def depthRNATuples = calculateDepthRNA(regionDepthRNA)
      
      // Filtrage des fichiers depth : ne garder que ceux qui contiennent au moins une base avec une profondeur < 50
      def filteredDepthDNATuples = depthDNATuples.filter { tuple ->
          def sample = tuple[0]
          def depthFile = tuple[1]
          def keep = false
          depthFile.eachLine { line ->
              def fields = line.tokenize()
              if( fields.size() >= 3 ) {
                  try {
                      if( fields[2].toInteger() < params.coverage_threshold ) {
                          keep = true
                      }
                  } catch(Exception e) { }
              }
          }
          return keep
      }
      filteredDepthDNATuples.subscribe { println "Filtered Depth DNA: $it" }
      
      def filteredDepthRNATuples = depthRNATuples.filter { tuple ->
          def sample = tuple[0]
          def depthFile = tuple[1]
          def keep = false
          depthFile.eachLine { line ->
              def fields = line.tokenize()
              if( fields.size() >= 3 ) {
                  try {
                      if( fields[2].toInteger() < params.coverage_threshold ) {
                          keep = true
                      }
                  } catch(Exception e) { }
              }
          }
          return keep
      }
      filteredDepthRNATuples.subscribe { println "Filtered Depth RNA: $it" }

      // ============================================================
      //  Rapports COVERAGE (HTML) pour ADN et ARN
      // ============================================================

      def groupedDepthDNA = filteredDepthDNATuples
                            .groupTuple(by:0)
                            .map { sample, files ->
                                tuple(sample,
                                      files,
                                      file(params.gtf_file),
                                      file(params.fasta_file),
                                      file("${params.fasta_file}.fai"))
                            }

      def groupedDepthRNA = filteredDepthRNATuples
                            .groupTuple(by:0)
                            .map { sample, files ->
                                tuple(sample,
                                      files,
                                      file(params.gtf_file),
                                      file(params.fasta_file),
                                      file("${params.fasta_file}.fai"))
                            }

      def coverageHtmlDNA = coverage2htmlDNA(groupedDepthDNA)
      def coverageHtmlRNA = coverage2htmlRNA(groupedDepthRNA)

      def chDuplicats = duplicatsTSO500(parentDir)
      chDuplicats.subscribe { println "Fichier Excel duplicats généré : $it" }
      
      // 3) excelToHtml => produit Resultats_duplicats.html
      def duplicatesHtml = excelToHtml(chDuplicats)
      duplicatesHtml.subscribe { println "Fichier HTML généré : $it" }

      // ============== Mélange des sorties HTML ==============
      def allHtml = compareOut
                      .mix(variantsHtml)
                      .mix(duplicatesHtml)
                      .mix(coverageHtmlDNA)
                      .mix(coverageHtmlRNA)
      def combinedResult = generateCombinedHtml(allHtml.collect())
      combinedResult.subscribe { println "Fichier combined.html généré : $it" }

  }
  // Branche secondaire pour ctBLUNG, ctPROST, LUCATISS
  else if( params.analysis_type in ['ctBLUNG','ctPROST','LUCATISS'] ) {
      println "Exécution du workflow secondaire pour ${params.analysis_type}"
      
      def run_dir = params.input_file
      if( !run_dir ) {
          log.info "Aucun dossier run spécifié. Ouverture de la boîte de dialogue..."
          def cmd = "python ${params.file_chooser_script}".execute()
          cmd.waitFor()
          def inputPath = cmd.in.text.trim()
          if( !inputPath || inputPath == 'null' ) {
              error "Aucun dossier sélectionné : arrêt du pipeline."
          }
          run_dir = inputPath
          println "Vous avez sélectionné : ${run_dir}"
      }
      println "Run directory pour l'analyse secondaire : ${run_dir}"
      
      params.output_dir         = "${run_dir}/QC_Analyse"
      params.taux_duplicats_dir = "${params.output_dir}/Taux_duplicats"
      params.combined_dir       = "${params.output_dir}/Rapport_Final"
      params.htm_dupli_dir      = "${params.taux_duplicats_dir}/Html"
      [ params.output_dir, params.taux_duplicats_dir, params.htm_dupli_dir ].each { new File(it).mkdirs() }
      
      def dupExcel
      if( params.analysis_type == 'LUCATISS' ) {
          dupExcel = duplicatsLucatiss(run_dir)
      } else {
          dupExcel = duplicatsCT(run_dir)
      }
      dupExcel.subscribe { println "Fichier Excel duplicats secondaire généré : $it" }
      
      def htmlOut = excelToHtml(dupExcel)
      htmlOut.subscribe { println "Fichier HTML duplicats secondaire généré : $it" }
  }
  else {
      error "Type d'analyse inconnu : ${params.analysis_type}"
  }
}

//////////////////////////////////////////////////////////////////////////
// Définition des Processus
//////////////////////////////////////////////////////////////////////////

//
// Process consolidateMetrics pour générer un fichier avec les valeurs seuils de référence pour les métriques
//
process consolidateMetrics {
    tag "Consolidate Metrics"
    stageInMode 'copy'
    stageOutMode 'copy'
    
    input:
        val parent_dir
    output:
        file "MetricsOutput_Reference_Seuils.tsv"
    script:
    """
    python ${params.consolidate_script} --parent_directory "${parent_dir}" --output_file "./MetricsOutput_Reference_Seuils.tsv" --data_dir "temp_MetricsOutput"
    """
}

//
// Process compareMetrics pour générer un fichier html du MetricsOutput du run qui est comparé aux valeurs de réf du fichier généré par le process consolidateMetrics
//
process compareMetrics {
    scratch true
    tag "Compare Metrics (HTML)"
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.combined_dir, mode: 'copy'
    input:
        path global_file
        path input_file
        path script_file
    output:
        file "MetricsComparison.html"
    script:
    """
    python ${script_file} --global_file "${global_file}" --input_file "${input_file}" --output_file "MetricsComparison.html"
    """
}

//
// Process listExonReports pour lister les fichiers exon_cov_report_bed pour les échantillons DNA 
//
process listExonReports {
    tag "List Exon Reports (DNA)"
    input:
        val logs_dir
    output:
        file "chemins_exon_cov_report_bed.txt"
    script:
    """
    exec > chemins_exon_cov_report_bed.txt 2>&1
    cd "${logs_dir}/DnaDragenCaller" || { echo "Le dossier DnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       report_file="${logs_dir}/DnaDragenCaller/\$sample/\${sample}.exon_cov_report.bed"
       if [ -f "\$report_file" ]; then
           realpath "\$report_file"
       fi
    done
    """
}

//
// Process listExonReadReports pour lister les fichiers exon_read_cov_report_bed pour les échantillons DNA  
//
process listExonReadReports {
    tag "List Exon Read Reports (DNA)"
    input:
        val logs_dir
    output:
        file "chemins_exon_read_cov_report_bed.txt"
    script:
    """
    exec > chemins_exon_read_cov_report_bed.txt 2>&1
    cd "${logs_dir}/DnaDragenCaller" || { echo "Le dossier DnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       read_file="${logs_dir}/DnaDragenCaller/\$sample/\${sample}.exon_read_cov_report.bed"
       if [ -f "\$read_file" ]; then
           realpath "\$read_file"
       fi
    done
    """
}

//
// Process listTargetReports pour lister les fichiers target_bed_read_cov_report_bed pour les échantillons DNA  
//
process listTargetReports {
    tag "List Target Bed Read Reports (DNA)"
    input:
        val logs_dir
    output:
        file "chemins_target_bed_read_cov_report_bed.txt"
    script:
    """
    exec > chemins_target_bed_read_cov_report_bed.txt 2>&1
    cd "${logs_dir}/DnaDragenCaller" || { echo "Le dossier DnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       target_file="${logs_dir}/DnaDragenCaller/\$sample/\${sample}.target_bed_read_cov_report.bed"
       if [ -f "\$target_file" ]; then
           realpath "\$target_file"
       fi
    done
    """
}

//
// Process processExonReports pour filtrer les exons qui ont des bases avec moins de 50X de coverage ( DNA )   
//
process processExonReports {
    tag "Process Exon Reports (DNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        path exon_paths_file
        path read_paths_file
        path target_paths_file
    output:
        path "*.tsv"
        path "*.bed"
    script:
    """
    python ${params.exon_report_script} \\
      ${exon_paths_file} \\
      ${read_paths_file} \\
      ${target_paths_file} \\
      .
    """
}

//
// Process filterCancerGenesDNA pour appliquer un second filtre sur les exons et ne garder que ceux lié au cancer de l'échantillon ( DNA ) 
//
process filterCancerGenesDNA {
    tag "Filter Cancer Genes (DNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        file report_files
    output:
        file("filtres_all_DNA/*")
    script:
    """
    mkdir -p tmp_reports
    cp ${report_files.join(' ')} tmp_reports/
    python3 ${params.filter_script} \\
         --input_dir tmp_reports \\
         --sample_cancer_file ${params.sample_cancer_file} \\
         --alias_cancer_file ${params.alias_cancer_file} \\
         --output_dir tmp_reports/filtres_all_DNA
    mv tmp_reports/filtres_all_DNA .
    """
}

//
// Process filterCancerGenesRNA pour appliquer un second filtre sur les exons et ne garder que ceux lié au cancer de l'échantillon ( RNA ) 
//
process filterCancerGenesRNA {
    tag "Filter Cancer Genes (RNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        file report_files
    output:
        file("filtres_all_RNA/*")
    script:
    """
    mkdir -p tmp_reports
    cp ${report_files*.toString().join(' ')} tmp_reports/
    python3 ${params.rna_filter_script} \\
         --input_dir tmp_reports \\
         --sample_cancer_file ${params.sample_cancer_file} \\
         --alias_cancer_file ${params.alias_cancer_file} \\
         --output_dir tmp_reports/filtres_all_RNA
    mv tmp_reports/filtres_all_RNA .
    """
}

//
// Process listExonReportsRna pour lister les fichiers exon_cov_report_bed pour les échantillons RNA  
//
process listExonReportsRna {
    tag "List Exon Reports (RNA)"
    input:
        val logs_dir
    output:
        file "chemins_exon_cov_report_bed_rna.txt"
    script:
    """
    exec > chemins_exon_cov_report_bed_rna.txt 2>&1
    cd "${logs_dir}/RnaDragenCaller" || { echo "Le dossier RnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       report_file="${logs_dir}/RnaDragenCaller/\$sample/\${sample}.exon_cov_report.bed"
       if [ -f "\$report_file" ]; then
           realpath "\$report_file"
       fi
    done
    """
}

//
// Process listExonReadReportsRna pour lister les fichiers exon_read_cov_report_bed pour les échantillons RNA  
//
process listExonReadReportsRna {
    tag "List Exon Read Reports (RNA)"
    input:
        val logs_dir
    output:
        file "chemins_exon_read_cov_report_bed_rna.txt"
    script:
    """
    exec > chemins_exon_read_cov_report_bed_rna.txt 2>&1
    cd "${logs_dir}/RnaDragenCaller" || { echo "Le dossier RnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       read_file="${logs_dir}/RnaDragenCaller/\$sample/\${sample}.exon_read_cov_report.bed"
       if [ -f "\$read_file" ]; then
           realpath "\$read_file"
       fi
    done
    """
}

//
// Process processExonReportsRna pour filtrer les exons qui ont des bases avec moins de 50X de coverage ( RNA )
//
process processExonReportsRna {
    tag "Process Exon Reports (RNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        path exon_paths_file_rna
        path read_paths_file_rna
    output:
        path "*.tsv"
    script:
    """
    python ${params.rna_exon_report_script} \\
       ${exon_paths_file_rna} \\
       ${read_paths_file_rna} \\
       .
    """
}

//
// Process listCombinedVariantOutputs pour lister les fichiers CombinedVariantOutput à utiliser
//
process listCombinedVariantOutputs {
    tag "List CombinedVariantOutput TSV"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        val results_dir
    output:
        file "chemins_combined_variant_outputs.txt"
    script:
    """
    exec > chemins_combined_variant_outputs.txt 2>&1
    cd "${results_dir}" || { echo "Le dossier Results n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       combined_file="${results_dir}/\$sample/\${sample}_CombinedVariantOutput.tsv"
       if [ -f "\$combined_file" ]; then
           realpath "\$combined_file"
       fi
    done
    """
}

//
// Process generateHtmlReport pour générer des fichiers html à partir des CombinedVariantOutput
// Ajout de la prise en compte du fichier gene_cytoband.txt
//
process generateHtmlReport {
    tag { "${tsv_path.toString().tokenize('/')[-1]}" }
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.combined_dir, mode: 'copy'
    input:
        val tsv_path
    output:
        file "*.html"
    script:
    """
    base=\$(basename "${tsv_path}" .tsv)
    python ${params.generate_html_script} \\
      --input_file "${tsv_path}" \\
      --output_file "\${base}_CombinedVariantOutput_DNA.html" \\
      --metadata_file "${params.sample_cancer_file}" \\
      --gene_cytoband_file "${params.gene_cytoband_file}"
    """
}

//
// Process listBamFilesDNA pour lister les fichiers BAM/BAI DNA à copier
//
process listBamFilesDNA {
    tag "List Bam Files (DNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        val logs_dir
    output:
        file "chemins_bam_files_DNA.txt"
    script:
    """
    exec > chemins_bam_files_DNA.txt 2>&1
    cd "${logs_dir}/DnaDragenCaller" || { echo "Le dossier DnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       bam_file="${logs_dir}/DnaDragenCaller/\$sample/\${sample}.bam"
       bai_file="${logs_dir}/DnaDragenCaller/\$sample/\${sample}.bam.bai"
       if [ -f "\$bam_file" ]; then
           realpath "\$bam_file"
       fi
       if [ -f "\$bai_file" ]; then
           realpath "\$bai_file"
       fi
    done
    """
}

//
// Process listBamFilesRNA pour lister les fichiers BAM/BAI RNA à copier
//
process listBamFilesRNA {
    tag "List Bam Files (RNA)"
    stageInMode 'copy'
    stageOutMode 'copy'
    input:
        val logs_dir
    output:
        file "chemins_bam_files_RNA.txt"
    script:
    """
    exec > chemins_bam_files_RNA.txt 2>&1
    cd "${logs_dir}/RnaDragenCaller" || { echo "Le dossier RnaDragenCaller n'existe pas"; exit 1; }
    for sample_dir in */ ; do
       sample=\$(basename "\$sample_dir")
       bam_file="${logs_dir}/RnaDragenCaller/\$sample/\${sample}.bam"
       bai_file="${logs_dir}/RnaDragenCaller/\$sample/\${sample}.bam.bai"
       if [ -f "\$bam_file" ]; then
           realpath "\$bam_file"
       fi
       if [ -f "\$bai_file" ]; then
           realpath "\$bai_file"
       fi
    done
    """
}

//
// Process reduceBamRegions_DNA pour ne cibler que certaines régions des bam DNA et gagner du temps
//
process reduceBamRegions_DNA {
    tag "\$sample"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    input:
        tuple val(sample), file(bam), file(region)
    output:
        tuple val(sample), file("${sample}_reduced.bam"), file("${sample}_reduced.bam.bai")
    script:
    """
    samtools view -b ${bam} -L ${region} > ${sample}_reduced.bam
    samtools index ${sample}_reduced.bam
    """
}

//
// Process reduceBamRegions_RNA pour ne cibler que certaines régions des bam RNA et gagner du temps
//
process reduceBamRegions_RNA {
    tag "\$sample"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    input:
        tuple val(sample), file(bam), file(region)
    output:
        tuple val(sample), file("${sample}_reduced.bam"), file("${sample}_reduced.bam.bai")
    script:
    """
    samtools view -b ${bam} -L ${region} > ${sample}_reduced.bam
    samtools index ${sample}_reduced.bam
    """
}

//
// Process filterSoftClip_Bam_DNA pour enlever les bases soft-clippés des bam DNA
//
process filterSoftClip_Bam_DNA {
    tag "\$sample"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    input:
        tuple val(sample), file(reducedBam), file(reducedBai)
    output:
        tuple val(sample), file("${sample}_filtered.bam"), file("${sample}_filtered.bam.bai")
    script:
    """
    python ${projectDir}/scripts/remove_softclip/clean_softclip.py --input ${reducedBam} --output ${sample}_filtered.bam
    samtools index ${sample}_filtered.bam
    """
}

//
// Process filterSoftClip_Bam_RNA pour enlever les bases soft-clippés des bam RNA
//
process filterSoftClip_Bam_RNA {
    tag "\$sample"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    input:
        tuple val(sample), file(reducedBam), file(reducedBai)
    output:
        tuple val(sample), file("${sample}_filtered.bam"), file("${sample}_filtered.bam.bai")
    script:
    """
    python ${projectDir}/scripts/remove_softclip/clean_softclip.py --input ${reducedBam} --output ${sample}_filtered.bam
    samtools index ${sample}_filtered.bam
    """
}

//
// Process calculateDepthDNA pour avoir les profondeurs de couverture pour échantillons ADN
//
process calculateDepthDNA {
    tag "\${sample}_\${name}"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    
    input:
        tuple val(sample), val(chrom), val(start), val(end), val(name), path(bam), path(bai)
    output:
        tuple val(sample), file("${sample}/${sample}_${name}_depth.txt")
    script:
    """
    mkdir -p ${sample}
    samtools depth -s -r ${chrom}:${start}-${end} -a -q 0 -Q 1 --excl-flags 0x400 ${bam} > ${sample}/${sample}_${name}_depth.txt
    """
}

//
// Process calculateDepthRNA pour avoir les profondeurs de couverture pour échantillons ARN
//
process calculateDepthRNA {
    tag "\${sample}_\${name}"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    
    input:
        tuple val(sample), val(chrom), val(start), val(end), val(name), path(bam), path(bai)
    output:
        tuple val(sample), file("${sample}/${sample}_${name}_depth.txt")
    script:
    """
    mkdir -p ${sample}
    samtools depth -s -r ${chrom}:${start}-${end} -a -q 0 -Q 1 --excl-flags 0x400 ${bam} > ${sample}/${sample}_${name}_depth.txt
    """
}

//
// Process duplicats pour le workflow TSO500
//
process duplicatsTSO500 {
    tag "Duplicats TSO500"
    stageInMode 'copy'
    stageOutMode 'copy'
    conda "$projectDir/environnement.yml"
    input:
        val run_directory
    output:
        file "Resultats_duplicats.xlsx"
    script:
    """
    python ${params.duplicats_script} ${run_directory}
    rsync -av ${run_directory}/QC_Analyse/Taux_duplicats/Resultats_duplicats.xlsx .
    """
}

//
// Process duplicats pour LUCATISS
//
process duplicatsLucatiss {
    tag "Duplicats LUCATISS"
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.taux_duplicats_dir, mode: 'copy'
    conda "$projectDir/environnement.yml"
    input:
        val run_directory
    output:
        file "Resultats_duplicats.xlsx"
    script:
    """
    python ${params.duplicats_lucatiss_script} ${run_directory}
    rsync -av ${run_directory}/QC_Analyse/Taux_duplicats/Resultats_duplicats.xlsx .
    """
}

//
// Process duplicats pour ctBLUNG / ctPROST
//
process duplicatsCT {
    tag "Duplicats ctBLUNG/ctPROST"
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.taux_duplicats_dir, mode: 'copy'
    conda "$projectDir/environnement.yml"
    input:
        val run_directory
    output:
        file "Resultats_duplicats.xlsx"
    script:
    """
    python ${params.duplicats_ct_script} ${run_directory}
    rsync -av ${run_directory}/QC_Analyse/Taux_duplicats/Resultats_duplicats.xlsx .
    """
}

//
// Process Excel to HTML (commun aux deux workflows)
//
process excelToHtml {
    tag "Excel to HTML"
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.htm_dupli_dir, mode: 'copy'
    publishDir { params.combined_dir }, mode: 'copy', when: { params.analysis_type in ['TSO500','TSO500_HRD'] }
    conda "$projectDir/environnement.yml"
    input:
        file excel_file
    output:
        file "Resultats_duplicats.html"
    script:
    """
    python ${params.excel_to_html_script} ${excel_file} Resultats_duplicats.html
    """
}

//
// Process generateCombinedHtml pour générer le rapport complet des analyses TSO500
//
process generateCombinedHtml {
    tag "Combine all HTML"
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.combined_dir, mode: 'copy'
    conda "$projectDir/environnement.yml"

    input:
        path html_files

    output:
        path "*.html"

    script:
    """
    python ${params.combine_html_script} . .
    """
}

//
// Process coverage2htmlDNA : génère un rapport HTML « coverage » par échantillon ADN
//
process coverage2htmlDNA {
    tag { "Coverage DNA ${sample}" }                // nom lisible dans le log
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.combined_dir, mode: 'copy'

    //------------------------------------------------------------------
    // ENTRÉES
    //------------------------------------------------------------------
    /*
       Le tuple contient :
         • sample        : nom de l’échantillon (val)
         • depth_files   : tous les *.depth.txt filtrés (path[])
         • gtf           : fichier d’annotation GTF (path)
         • fasta         : séquence FASTA (path)
         • fai           : index .fai (path optionnel)
    */
    input:
        tuple val(sample), path(depth_files), path(gtf), path(fasta), path(fai)

    //------------------------------------------------------------------
    // SORTIE
    //------------------------------------------------------------------
    output:
        file "${sample}_coverage.html"

    //------------------------------------------------------------------
    // SCRIPT
    //------------------------------------------------------------------
    script:
    """
    # 1) Répertoire temporaire et sous‑dossier échantillon
    tmp=\$(mktemp -d)
    mkdir -p \$tmp/${sample}
    cp ${depth_files.join(' ')} \$tmp/${sample}/

    # 2) Les références sont déjà « stagées » ici ; on retient juste leur nom
    gtf_basename=\$(basename ${gtf})
    fasta_basename=\$(basename ${fasta})

    # 3) Dossier de sortie local
    outdir=\$(pwd)/outdir
    mkdir -p \$outdir

    # 4) Appel du script Python (inchangé)
    python ${projectDir}/scripts/low_coverage/coverage2html_V2.py \\
          \$gtf_basename \\
          \$fasta_basename \\
          \$tmp \\
          \$outdir

    # 5) Renommage/placement du fichier pour Nextflow
    mv \$outdir/${sample}.html ${sample}_coverage.html
    """
}

//
// Process coverage2htmlRNA : génère un rapport HTML « coverage » par échantillon ARN
//
process coverage2htmlRNA {

    tag { "Coverage RNA ${sample}" }
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir params.combined_dir, mode: 'copy'

    //------------------------------------------------------------------
    // ENTRÉES
    //------------------------------------------------------------------
    input:
        tuple val(sample), path(depth_files), path(gtf), path(fasta), path(fai)

    //------------------------------------------------------------------
    // SORTIE
    //------------------------------------------------------------------
    output:
        file "${sample}_coverage.html"

    //------------------------------------------------------------------
    // SCRIPT
    //------------------------------------------------------------------
    script:
    """
    # 1) Répertoire temporaire et sous‑dossier échantillon
    tmp=\$(mktemp -d)
    mkdir -p \$tmp/${sample}
    cp ${depth_files.join(' ')} \$tmp/${sample}/

    # 2) Noms de base des références déjà présentes
    gtf_basename=\$(basename ${gtf})
    fasta_basename=\$(basename ${fasta})

    # 3) Dossier de sortie local
    outdir=\$(pwd)/outdir
    mkdir -p \$outdir

    # 4) Appel du script Python
    python ${projectDir}/scripts/low_coverage/coverage2html_V2.py \\
          \$gtf_basename \\
          \$fasta_basename \\
          \$tmp \\
          \$outdir

    # 5) Fichier final renommé
    mv \$outdir/${sample}.html ${sample}_coverage.html
    """
}

