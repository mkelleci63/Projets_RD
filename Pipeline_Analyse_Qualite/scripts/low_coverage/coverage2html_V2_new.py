#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import json
from collections import defaultdict
import csv 
from pathlib import Path

###############################################################################
# Dictionnaire RefSeq -> chromosome (GRCh37/hg19)
###############################################################################
refseq_to_chrom = {
    "NC_000001.10": "chr1",
    "NC_000002.11": "chr2",
    "NC_000003.11": "chr3",
    "NC_000004.11": "chr4",
    "NC_000005.9":  "chr5",
    "NC_000006.11": "chr6",
    "NC_000007.13": "chr7",
    "NC_000008.10": "chr8",
    "NC_000009.11": "chr9",
    "NC_000010.10": "chr10",
    "NC_000011.9":  "chr11",
    "NC_000012.11": "chr12",
    "NC_000013.10": "chr13",
    "NC_000014.8":  "chr14",
    "NC_000015.9":  "chr15",
    "NC_000016.9":  "chr16",
    "NC_000017.10": "chr17",
    "NC_000018.9":  "chr18",
    "NC_000019.9":  "chr19",
    "NC_000020.10": "chr20",
    "NC_000021.8":  "chr21",
    "NC_000022.10": "chr22",
    "NC_000023.10": "chrX",
    "NC_000024.9":  "chrY"
}

###############################################################################
# 1) Parsing du GTF et FASTA
###############################################################################
def parse_gtf_pour_plusieurs_genes(gtf_file, genes_of_interest):
    """
    Lit le fichier GTF une seule fois pour tous les gènes de genes_of_interest.
    Retourne un dict { gene_id: { transcript_id: données_du_transcrit } }.
    """
    from collections import defaultdict
    import re

    # Préparer le conteneur par gène
    all_transcripts = {
        gene_id: defaultdict(lambda: {
            'strand': None,
            'chrom': None,
            'original_chrom': None,
            'exons': [],
            'CDS': [],
            'transcript_start': None,
            'transcript_end': None
        })
        for gene_id in genes_of_interest
    }

    # Lecture ligne à ligne du GTF
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_str, end_str, score, strand, frame, attrs = parts
            start = int(start_str)
            end   = int(end_str)

            # Parse des attributs (clé "valeur")
            attr_dict = {
                key: value
                for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attrs)
            }
            gene_id = attr_dict.get('gene') or attr_dict.get('gene_id')
            if gene_id not in all_transcripts:
                continue

            transcript_id = attr_dict.get('transcript_id')
            if not transcript_id or transcript_id.startswith(('XR','NR')):
                continue

            tr = all_transcripts[gene_id][transcript_id]

            # Enregistrer le chromosome original, puis le remapper si besoin
            tr['original_chrom'] = chrom
            tr['chrom'] = refseq_to_chrom.get(chrom, chrom)

            # Enregistrer le brin une fois
            if tr['strand'] is None:
                tr['strand'] = strand

            # Enregistrer selon le type de feature
            feat = feature.lower()
            if feat == 'transcript':
                tr['transcript_start'] = start
                tr['transcript_end']   = end
            elif feat == 'exon':
                exon_number = attr_dict.get('exon_number', '')
                tr['exons'].append({
                    'start': start,
                    'end': end,
                    'exon_number': exon_number
                })
            elif feat == 'cds':
                exon_number = attr_dict.get('exon_number', '')
                tr['CDS'].append({
                    'start': start,
                    'end': end,
                    'exon_number': exon_number
                })

    # Post‑traitement : complétion des bornes transcript et découpe UTR/CDS
    def merge_intervals(intervals):
        if not intervals:
            return []
        intervals = sorted(intervals, key=lambda x: x[0])
        merged = [intervals[0]]
        for curr in intervals[1:]:
            last = merged[-1]
            if curr[0] <= last[1] + 1:
                merged[-1] = (last[0], max(last[1], curr[1]))
            else:
                merged.append(curr)
        return merged

    for gene_id, transcripts in all_transcripts.items():
        for tid, tr in transcripts.items():
            # transcript_start/end manquants
            if tr['transcript_start'] is None or tr['transcript_end'] is None:
                if tr['exons']:
                    starts = [e['start'] for e in tr['exons']]
                    ends   = [e['end']   for e in tr['exons']]
                    tr['transcript_start'] = min(starts)
                    tr['transcript_end']   = max(ends)
                else:
                    tr['transcript_start'] = 0
                    tr['transcript_end']   = 0

            # regrouper CDS par exon_number
            cds_by_exon = {}
            for cds in tr['CDS']:
                key = cds['exon_number']
                cds_by_exon.setdefault(key, []).append((cds['start'], cds['end']))

            # pour chaque exon, construire la liste regions UTR/CDS
            for exon in tr['exons']:
                exon_num   = exon['exon_number']
                exon_start = exon['start']
                exon_end   = exon['end']
                merged_cds = merge_intervals(cds_by_exon.get(exon_num, []))
                regions = []
                cursor = exon_start
                for cs, ce in merged_cds:
                    if cursor < cs:
                        regions.append({'type':'UTR','start':cursor,'end':cs-1})
                    regions.append({'type':'CDS','start':cs,'end':ce})
                    cursor = ce + 1
                if cursor <= exon_end:
                    regions.append({'type':'UTR','start':cursor,'end':exon_end})
                exon['regions'] = regions

    return all_transcripts


def select_canonical_transcript(transcripts):
    best_tid, best_len = None, -1
    for tid, data in transcripts.items():
        cds_len = sum((cd["end"] - cd["start"] + 1) for cd in data["CDS"])
        data["cds_length"] = cds_len
        if cds_len > best_len or (cds_len == best_len and (best_tid is None or tid < best_tid)):
            best_tid, best_len = tid, cds_len
    return best_tid

def build_reference_exons(canonical_exons, strand):
    sorted_exons = sorted(canonical_exons, key=lambda x: x["start"])
    if strand == "-":
        sorted_exons.reverse()
    for exon in sorted_exons:
        exon["genomicWidth"] = exon["end"] - exon["start"] + 1
    return sorted_exons

def read_fasta_region(fasta_file, chromosome, start, end, strand):
    if start < 1:
        start = 1
    if end < start:
        return ""
    found_chrom = False
    seq_chunks = []
    with open(fasta_file, "r") as f:
        current_chrom = None
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                current_chrom = header.split()[0]
                if found_chrom and current_chrom != chromosome:
                    break
                found_chrom = (current_chrom == chromosome)
                seq_chunks = []
            else:
                if found_chrom:
                    seq_chunks.append(line.strip())
    if not found_chrom:
        return ""
    full_chrom_seq = "".join(seq_chunks)
    slice_part = full_chrom_seq[start-1:end] if end <= len(full_chrom_seq) else full_chrom_seq[start-1:]
    return slice_part

###############################################################################
# 1-b) Lecture des fichiers de coverage et extraction des intervalles low coverage
###############################################################################
def read_coverage_intervals(file_path):
    positions = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = re.split(r'\s+', line)
            if len(cols) < 3:
                continue
            try:
                pos = int(cols[1])
                depth = float(cols[2])
            except ValueError:
                continue
            if depth < 50:
                positions.append(pos)
    intervals = []
    if positions:
        positions.sort()
        start = positions[0]
        end = positions[0]
        for pos in positions[1:]:
            if pos == end + 1:
                end = pos
            else:
                intervals.append((start, end))
                start = pos
                end = pos
        intervals.append((start, end))
    return intervals

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        if current[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged

def read_variants_file(variants_file):
    """
    Lit le fichier de variants, n’extrait que la section [Small Variants],
    ignore la ligne d’en-tête, et renvoie une liste de dict { 'chrom': str, 'pos': int }.
    """
    variants = []
    in_section = False

    with open(variants_file, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # on saute les lignes vides
            if not row:
                continue

            # détection du début/fin de section
            if row[0].startswith("["):
                if row[0].strip() == "[Small Variants]":
                    in_section = True
                elif in_section:
                    # on quitte dès qu’on rencontre une autre rubrique
                    break
                continue

            if not in_section:
                continue

            # on ignore la ligne d’en-tête (la première ligne après [Small Variants])
            if row[0] == "Gene":
                continue

            # vérification qu’il y a au moins 3 colonnes
            if len(row) < 3:
                continue

            chrom = row[1]
            try:
                pos = int(row[2])
            except ValueError:
                # si la position n’est pas un entier, on passe
                continue

            variants.append({"chrom": chrom, "pos": pos})

    return variants

###############################################################################
# 2) Génération du HTML statique (vues globale, zoomée, détaillée)
###############################################################################
def generate_html(genes_data, output_file):
    json_data = json.dumps({"genes": genes_data}, indent=2)
    coverageColor = "#FFCC99"
    html = f"""<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8">
  <title>Visualisation des gènes</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
    body {{
      margin: 0;
      padding: 0;
      font-family: Arial, sans-serif;
    }}
    .fixed-header {{
      position: fixed;
      top: 0;
      left: 0;
      right: 0;
      height: 60px;
      background: #fff;
      border-bottom: 1px solid #ccc;
      display: flex;
      align-items: center;
      justify-content: center;
      font-size: 16px;
      font-weight: bold;
      z-index: 10;
    }}
    #geneSelector {{
      font-size: 16px;
      padding: 4px 8px;
    }}
    .header-right {{
      position: absolute;
      right: 20px;
      top: 15px;
    }}
    .header-right button {{
      font-size: 14px;
      padding: 4px 8px;
    }}
    .drawing-container, .sequence-container {{
      margin-top: 60px;
      height: calc(100vh - 60px);
      overflow: auto;
      background: #fafafa;
    }}
    .sequence-container {{
      display: none;
    }}
    .x.axis path,
    .x.axis line {{
      fill: none;
      stroke: #000;
      shape-rendering: crispEdges;
    }}
    .x.axis text {{
      fill: black;
      font-size: 10px;
    }}
    .clickable {{
      cursor: pointer;
    }}
    .arrow-lowcov {{
      fill: red;
      cursor: pointer;
    }}
    .close-btn {{
      float: right;
      cursor: pointer;
      font-weight: bold;
      margin-left: 10px;
    }}
    .tooltip {{
      position: absolute;
      padding: 8px;
      background: rgba(0,0,0,0.8);
      color: white;
      border-radius: 4px;
      font-size: 12px;
      pointer-events: none;
      opacity: 0;
      z-index: 100;
      max-width: 300px;
      line-height: 1.4em;
    }}
  </style>
</head>
<body>

<div class="fixed-header">
  <select id="geneSelector"></select>
  <div class="header-right" id="backButtonContainer" style="display: none;">
    <button id="backButton">Retour</button>
  </div>
</div>

<div class="drawing-container" id="chartContainer">
  <svg id="chart"></svg>
</div>

<div class="sequence-container" id="sequenceContainer">
  <svg id="sequenceSvg"></svg>
</div>

<script>
  var tooltip = d3.select("body").append("div")
      .attr("class", "tooltip");

  document.addEventListener("click", function(e) {{
    if (
      !e.target.classList.contains("base-rect") &&
      !e.target.classList.contains("tooltip")
    ) {{
      tooltip.style("opacity", 0);
    }}
  }});

  var marginLeft = 300, marginTop = 80, marginRight = 50, rowHeight = 60;
  var exonHeight = 30, intronWidth = 2;
  var coverageColor = "{coverageColor}";
  var currentView = "global";
  var lastZoomState = {{}};

  var xScale;
  var width = window.innerWidth;

  var svg = d3.select("#chart");
  svg.attr("width", width);
  var transcriptsGroup = svg.append("g").attr("id", "transcriptsGroup");

  var allGenesData = {json_data}.genes;
  var geneIds = Object.keys(allGenesData);
  var currentGeneId = geneIds[0];
  var currentGeneData = allGenesData[currentGeneId];
  var currentTranscript = currentGeneData.canonical;

  var geneSelector = d3.select("#geneSelector");
  geneIds.forEach(function(geneId) {{
    var gene = allGenesData[geneId];
    var optionText = "Gene: " + geneId +
                     " | Chr: " + gene.canonical.chrom +
                     " | Brin: " + gene.canonical.strand;
    geneSelector.append("option")
      .attr("value", geneId)
      .text(optionText);
  }});
  geneSelector.on("change", function() {{
    currentGeneId = this.value;
    currentGeneData = allGenesData[currentGeneId];
    currentTranscript = currentGeneData.canonical;
    currentView = "global";
    d3.select("#chartContainer").style("display", "block");
    d3.select("#sequenceContainer").style("display", "none");
    d3.select("#backButtonContainer").style("display", "none");
    renderGlobal();
  }});

  d3.select("#backButton").on("click", function() {{
    if (currentView === "detail") {{
       if (lastZoomState.type === "exon") {{
         renderZoomedExonView(lastZoomState);
       }} else if (lastZoomState.type === "feature") {{
         renderZoomedFeatureView(lastZoomState);
       }}
       currentView = "zoomed";
    }} else if (currentView === "zoomed") {{
      currentView = "global";
      d3.select("#chartContainer").style("display", "block");
      d3.select("#sequenceContainer").style("display", "none");
      d3.select("#backButtonContainer").style("display", "none");
      renderGlobal();
    }}
  }});

  function drawAxis(svgContainer, scale, yPos) {{
      /* axe (inchangé) */
      svgContainer.selectAll(".x.axis").remove();
      var axisG = svgContainer.append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(0," + yPos + ")")
          .call(
              d3.axisTop(scale)
                 .ticks(10)
                 .tickFormat(d3.format("d"))
          );
      axisG.selectAll("path").style("stroke", "black").style("fill", "none");
      axisG.selectAll("line").style("stroke", "black");

      /* légende “Position génomique” : alignée à droite, 10 px en dessous de l’axe */
      svgContainer.selectAll(".pos-label").remove();
      svgContainer.append("text")
          .attr("class", "pos-label")
          .attr("x", width - marginRight)   /* bord droit de l’échelle */
          .attr("y", yPos + 20)             /* 10 px en dessous de la ligne d’axe */
          .attr("text-anchor", "end")
          .attr("font-size", "10px")
          .attr("fill", "black")
          .text("Position génomique");
  }}
  
  function computeGlobalDomain() {{
    var allTranscripts = [currentGeneData.canonical].concat(currentGeneData.transcripts);
    var start = d3.min(allTranscripts, function(t) {{ return t.transcript_start; }});
    var end   = d3.max(allTranscripts, function(t) {{ return t.transcript_end; }});
    var span = end - start;
    var extraMargin = span * 0.1;
    return [start - extraMargin, end + extraMargin];
  }}

  function renderGlobal() {{
    currentView = "global";
    svg.selectAll("*").remove();
    xScale = d3.scaleLinear()
               .domain(computeGlobalDomain())
               .range([marginLeft, width - marginRight]);
    drawAxis(svg, xScale, 40);
    transcriptsGroup = svg.append("g").attr("id", "transcriptsGroup");
    var rowIndex = 0;
    transcriptsGroup.append("text")
       .attr("x", 10)
       .attr("y", marginTop + rowIndex * rowHeight + rowHeight/2 + 5)
       .attr("font-size", 16)
       .attr("font-weight", "bold")
       .text("Transcrit canonique");
    rowIndex++;
    renderTranscript(currentGeneData.canonical, rowIndex, true);
    rowIndex++;
    transcriptsGroup.append("text")
       .attr("x", 10)
       .attr("y", marginTop + rowIndex * rowHeight + rowHeight/2 + 5)
       .attr("font-size", 16)
       .attr("font-weight", "bold")
       .text("Transcrits alternatifs (" + currentGeneData.transcripts.length + ")");
    rowIndex++;
    currentGeneData.transcripts.forEach(function(tr) {{
      renderTranscript(tr, rowIndex, true);
      rowIndex++;
    }});
    svg.attr("height", marginTop + rowIndex * rowHeight + 100);
  }}

  function renderTranscript(transcript, rowIndex, showName) {{
    var y = marginTop + rowIndex * rowHeight;
    var currentDomain = xScale.domain();

    if (showName) {{
      var labelGroup = transcriptsGroup.append("g")
          .attr("transform", "translate(10," + (y + rowHeight/2 + 5) + ")");
      labelGroup.append("text")
         .attr("x", 0)
         .attr("y", 0)
         .attr("font-size", 16)
         .attr("font-weight", "normal")
         .text(transcript.transcript_id);
    }}

    var exons = transcript.exons.slice().sort(function(a, b) {{ return a.start - b.start; }});

    exons.forEach(function(ex, i) {{
      var clampStart = Math.max(ex.start, currentDomain[0]);
      var clampEnd = Math.min(ex.end, currentDomain[1]);
      if (clampEnd < clampStart) return;
      var px1 = xScale(clampStart);
      var px2 = xScale(clampEnd + 1);
      transcriptsGroup.append("rect")
         .attr("x", px1)
         .attr("y", y + (rowHeight - exonHeight)/2)
         .attr("width", px2 - px1)
         .attr("height", exonHeight)
         .attr("fill", "steelblue");

      if (i < exons.length - 1) {{
        var nextEx = exons[i+1];
        var intronStart = ex.end + 1;
        var intronEnd = nextEx.start - 1;
        var clampIStart = Math.max(intronStart, currentDomain[0]);
        var clampIEnd = Math.min(intronEnd, currentDomain[1]);
        if (clampIEnd > clampIStart) {{
          transcriptsGroup.append("line")
            .attr("x1", xScale(clampIStart))
            .attr("x2", xScale(clampIEnd))
            .attr("y1", y + rowHeight/2)
            .attr("y2", y + rowHeight/2)
            .attr("stroke", "black")
            .attr("stroke-width", intronWidth);
        }}
      }}

      var coverageIntervals = currentGeneData.coverage ? currentGeneData.coverage.intervals : null;
      var exonHasLowCov = false;
      if (coverageIntervals) {{
          coverageIntervals.forEach(function(interval) {{
              if (interval.end >= ex.start && interval.start <= ex.end) {{
                  exonHasLowCov = true;
              }}
          }});
      }}
      if (exonHasLowCov) {{
          var centralPos = (ex.start + ex.end) / 2;
          var arrowX = xScale(centralPos) - 5;
          var arrowY = y + (rowHeight - exonHeight) / 2 - 10;
          transcriptsGroup.append("path")
            .attr("d", "M0,0 L10,0 L5,10 Z")
            .attr("transform", "translate(" + arrowX + "," + arrowY + ")")
            .attr("class", "arrow-lowcov")
            .on("click", function(event) {{
              zoomToExon(ex, transcript);
              event.stopPropagation();
            }});
      }}
    }});
  }}

  function findContiguousBlock(clickedPos, dataArray) {{
    var block = {{start: clickedPos, end: clickedPos}};
    for (var i = 0; i < dataArray.length; i++) {{
      if (dataArray[i].pos === clickedPos) {{
        var j = i;
        while (j > 0 && dataArray[j].pos - dataArray[j-1].pos === 1) {{
          block.start = dataArray[j-1].pos;
          j--;
        }}
        j = i;
        while (j < dataArray.length - 1 && dataArray[j+1].pos - dataArray[j].pos === 1) {{
          block.end = dataArray[j+1].pos;
          j++;
        }}
        break;
      }}
    }}
    return block;
  }}

  function zoomToLowCoverageDetail(rangeStart, rangeEnd, transcript) {{
  currentView = "detail";
  var newStart = rangeStart - 20;
  var newEnd   = rangeEnd + 20;
  if (newStart < transcript.transcript_start) newStart = transcript.transcript_start;
  if (newEnd   > transcript.transcript_end)   newEnd   = transcript.transcript_end;

  var seqSVG = d3.select("#sequenceSvg");
  seqSVG.html("");
  seqSVG.attr("width", width).attr("height", window.innerHeight);

  var xScaleDetail = d3.scaleLinear()
       .domain([newStart, newEnd])
       .range([100, width - marginRight]);

  var axisY = marginTop + 20;
  drawAxis(seqSVG, xScaleDetail, axisY);

  var labelText = "Détail de coverage low (<50) pour " + transcript.transcript_id +
                  " | Positions : " + newStart + " - " + newEnd;
  seqSVG.append("text")
        .attr("x", 10)
        .attr("y", marginTop - 20)
        .attr("font-size", 16)
        .text(labelText);

  d3.select("#backButtonContainer").style("display", "block");
  d3.select("#chartContainer").style("display", "none");
  d3.select("#sequenceContainer").style("display", "block");

  var histMarginTop = marginTop + 60;
  var histHeight    = 300;
  drawCoverageHistogramDetail(seqSVG, xScaleDetail, newStart, newEnd, histMarginTop, histHeight);

  var squaresY = marginTop + histMarginTop + histHeight + 40 - rowHeight/2;
  window.squaresY = squaresY;

  var nucleotideGroup = seqSVG.append("g").attr("id", "nucleotideSquares");
  for (var pos = newStart; pos <= newEnd; pos++) {{
    (function(pos) {{
      var index = pos - transcript.transcript_start;
      var nuc   = transcript.full_sequence.charAt(index) || "N";
      var x     = xScaleDetail(pos);
      var nextX = xScaleDetail(pos + 1);
      var squareWidth = nextX - x;

      var regionType = 'intron';
      for (var i = 0; i < transcript.exons.length; i++) {{
        var exon = transcript.exons[i];
        if (pos >= exon.start && pos <= exon.end) {{
          exon.regions.forEach(function(region) {{
            if (pos >= region.start && pos <= region.end) {{
              regionType = region.type;
            }}
          }});
          break;
        }}
      }}

      var fillColor = regionType === 'CDS' ? '#ADD8E6'
                      : (regionType === 'UTR' ? '#FF0000' : 'gray');

      var handler = function(event) {{
        event.stopPropagation();
        var baseDepth = 0;
        if (currentGeneData.coverage_profile) {{
          for (var k = 0; k < currentGeneData.coverage_profile.length; k++) {{
            if (currentGeneData.coverage_profile[k].pos === pos) {{
              baseDepth = currentGeneData.coverage_profile[k].depth;
              break;
            }}
          }}
        }}
        var chrom = transcript.chrom;
        tooltip.html(
          "Chromosome: " + chrom + "<br>" +
          "Allèle: "       + nuc   + "<br>" +
          "Position: "    + pos   + "<br>" +
          "Profondeur: "  + baseDepth
        )
        .style("opacity", 1);
        var ttW = tooltip.node().getBoundingClientRect().width;
        var ttH = tooltip.node().getBoundingClientRect().height;
        var px  = event.pageX;
        var py  = event.pageY;
        if (px + ttW > window.innerWidth) {{
          px -= ttW + 20;
        }}
        tooltip.style("left", px + "px")
               .style("top", (py - ttH - 10) + "px");
      }};

      nucleotideGroup.append("rect")
        .attr("x", x)
        .attr("y", squaresY)
        .attr("width", squareWidth)
        .attr("height", exonHeight)
        .attr("fill", fillColor)
        .attr("class", "base-rect")
        .on("click", handler);

      nucleotideGroup.append("text")
        .attr("x", x + squareWidth/2)
        .attr("y", squaresY + exonHeight/2 + 4)
        .attr("text-anchor", "middle")
        .attr("font-size", "10px")
        .attr("fill", regionType === 'CDS' ? 'black' : 'white')
        .text(nuc)
        .on("click", handler);
    }})(pos);
  }}

  // ─── Affichage des variants (flèches violettes) ───
  var arrowYVariants = squaresY - 10;
  currentGeneData.variants
    .filter(v => v.pos >= newStart && v.pos <= newEnd)
    .forEach(function(v) {{
      var xBase   = xScaleDetail(v.pos);
      var xNext   = xScaleDetail(v.pos + 1);
      var centerX = xBase + (xNext - xBase) / 2;
      var vx      = centerX - 5;
      seqSVG.append("path")
        .attr("d", "M0,0 L10,0 L5,10 Z")
        .attr("transform", "translate(" + vx + "," + arrowYVariants + ")")
        .attr("fill", "purple");
    }});
  }}



  function drawCoverageHistogramDetail(svgContainer,
                                       localScale,
                                       seqStart,
                                       seqEnd,
                                       histMarginTop,
                                       histHeight) {{
      var profile = currentGeneData.coverage_profile;
      if (!profile) return;

      var filtered = profile.filter(function(d) {{ return d.pos >= seqStart && d.pos <= seqEnd; }});
      if (filtered.length === 0) return;

      var globalMax = d3.max(currentGeneData.coverage_profile, function(d) {{ return d.depth; }});
      // Étendre l’axe Y de 15px vers le bas
      var yScale = d3.scaleLinear()
          .domain([0, globalMax])
          .range([histMarginTop + histHeight + 35,   // +15px ici (était +20)
                  histMarginTop - 20]);

      var lowCoverageData = filtered.filter(function(d) {{ return d.depth < 50; }})
                                    .sort(function(a, b) {{ return a.pos - b.pos; }});

      /* barres */
      filtered.forEach(function(d) {{
          var x = localScale(d.pos);
          var barEnd = (d.pos + 1 > seqEnd) ? seqEnd : d.pos + 1;
          var barW = localScale(barEnd) - localScale(d.pos);
          var y = yScale(d.depth);
          // Étendre la hauteur des barres de 15px vers le bas
          var barH = (histMarginTop + histHeight + 35) - y;
          var color = (d.depth >= 50) ? "blue" : "#FFA500";

          var rect = svgContainer.append("rect")
              .attr("x", x)
              .attr("y", y)
              .attr("width", barW)
              .attr("height", barH)
              .attr("fill", color);

          if (d.depth < 50) {{
              rect.style("cursor", "pointer")
                  .on("click", function(event) {{
                      var block = findContiguousBlock(d.pos, lowCoverageData);
                      zoomToLowCoverageDetail(block.start, block.end, currentTranscript);
                      event.stopPropagation();
                  }});
          }}
      }});

      /* axe Y */
      svgContainer.append("g")
          .attr("class", "y axis")
          .attr("transform", "translate(50,0)")
          .call(d3.axisLeft(yScale).ticks(5));

      /* seuil 50 (inchangé) */  
      svgContainer.append("text")
          .attr("x", 40)
          .attr("y", yScale(50) + 3)
          .attr("text-anchor", "end")
          .attr("font-size", "10px")
          .attr("fill", "red")
          .attr("font-weight", "bold")
          .text("50");
      svgContainer.append("line")
          .attr("x1", 45).attr("x2", 50)
          .attr("y1", yScale(50)).attr("y2", yScale(50))
          .attr("stroke", "red").attr("stroke-width", 2);
      var scaleRange = xScale.range();
      svgContainer.append("line")
          .attr("x1", scaleRange[0]).attr("x2", scaleRange[1])
          .attr("y1", yScale(50)).attr("y2", yScale(50))
          .attr("stroke", "red").attr("stroke-dasharray", "4,4");

      /* légende axe Y : recentrée sur la nouvelle hauteur */
      svgContainer.append("text")
          .attr("class", "cov-label")
          .attr("transform", "rotate(-90)")
          .attr("x", -(histMarginTop + histHeight / 2 + 7.5))  // +7.5px pour recentrage
          .attr("y", 15)
          .attr("text-anchor", "middle")
          .attr("font-size", "10px")
          .attr("fill", "black")
          .text("Profondeur de couverture (X)");
  }}
    
  function zoomToExon(exon, transcript) {{
    currentView = "zoomed";
    d3.select("#backButtonContainer").style("display", "block");

    const PADDING = 100;
    let seqStart = exon.start - PADDING;
    let seqEnd   = exon.end   + PADDING;
    if (seqStart < transcript.transcript_start) seqStart = transcript.transcript_start;
    if (seqEnd   > transcript.transcript_end)   seqEnd   = transcript.transcript_end;

    xScale.domain([seqStart, seqEnd]).range([100, width - marginRight]);

    // bornes CDS du transcrit
    const cdsStart = d3.min(transcript.CDS, d => d.start);
    const cdsEnd   = d3.max(transcript.CDS, d => d.end);

    var seqSVG = d3.select("#sequenceSvg");
    seqSVG.html("");
    seqSVG.attr("width", width).attr("height", window.innerHeight);

    var axisY = marginTop + 20;
    drawAxis(seqSVG, xScale, axisY);

    var labelText = (transcript.transcript_id === currentGeneData.canonical.transcript_id ?
                     "Transcrit canonique : " : "Transcrit alternatif : ") +
                    transcript.transcript_id +
                    " | Positions : " + seqStart + " - " + seqEnd;
    seqSVG.append("text")
          .attr("x", 10)
          .attr("y", marginTop - 20)
          .attr("font-size", "16px")
          .text(labelText);

    var histMarginTop = marginTop + 60;
    var histHeight    = 300;
    drawCoverageHistogram(seqStart, seqEnd, histMarginTop, histHeight);

    var seqGroup = seqSVG.append("g").attr("id", "seqTranscriptGroup");

    // ─── Exons/UTR ───
    exon.regions.forEach(function(region) {{
        var rStart = Math.max(region.start, seqStart);
        var rEnd   = Math.min(region.end,   seqEnd);
        if (rStart > rEnd) return;
        var px1 = xScale(rStart);
        var px2 = xScale(rEnd + 1);
        var w   = px2 - px1;
        var color = region.type === "CDS" ? "#ADD8E6" : "#FF0000";
        seqGroup.append("rect")
                .attr("x", px1)
                .attr("y", marginTop + histMarginTop + histHeight + 40 - rowHeight/2)
                .attr("width",  w)
                .attr("height", exonHeight)
                .attr("fill",   color);

        if (w > 12) {{
            let label = region.type;
            if (region.type === "UTR") {{
                const fivePrime = (transcript.strand === "+" ?
                                   (region.end < cdsStart) :
                                   (region.start > cdsEnd));
                label = fivePrime ? "5'UTR" : "3'UTR";
            }}
            var fSize = Math.max(8, Math.min(14, w / label.length));
            seqGroup.append("text")
                    .attr("x", px1 + w/2)
                    .attr("y", marginTop + histMarginTop + histHeight + 40 - rowHeight/2 + exonHeight/2 + 4)
                    .attr("text-anchor", "middle")
                    .attr("font-size", fSize + "px")
                    .attr("fill", region.type === "CDS" ? "black" : "white")
                    .text(label);
        }}
    }});

    // ─── Introns + marqueurs 10/20/50 ───
    var baseY = marginTop + histMarginTop + histHeight + 40 - exonHeight/2;
    if (exon.start > transcript.transcript_start) {{
        seqGroup.append("line")
                .attr("x1", xScale(seqStart))
                .attr("x2", xScale(exon.start))
                .attr("y1", baseY)
                .attr("y2", baseY)
                .attr("stroke", "black")
                .attr("stroke-width", 2);
        drawMarkerAtBoundary("left", exon.start, transcript, baseY, seqGroup, seqSVG);
    }}
    if (exon.end < transcript.transcript_end) {{
        seqGroup.append("line")
                .attr("x1", xScale(exon.end + 1))
                .attr("x2", xScale(seqEnd))
                .attr("y1", baseY)
                .attr("y2", baseY)
                .attr("stroke", "black")
                .attr("stroke-width", 2);
        drawMarkerAtBoundary("right", exon.end, transcript, baseY, seqGroup, seqSVG);
    }}

    // ─── Affichage des variants (flèches violettes) ───
    var arrowY = marginTop + histMarginTop + histHeight + 40 - rowHeight/2 - 10;
    currentGeneData.variants
        .filter(v => v.pos >= seqStart && v.pos <= seqEnd)
        .forEach(function(v) {{
            var vx = xScale(v.pos) - 5;
            seqGroup.append("path")
                    .attr("d", "M0,0 L10,0 L5,10 Z")
                    .attr("transform", "translate(" + vx + "," + arrowY + ")")
                    .attr("fill", "purple");
        }});

    lastZoomState = {{
        type:       "exon",
        transcript: transcript,
        seqStart:   seqStart,
        seqEnd:     seqEnd,
        exon:       exon
    }};

    d3.select("#chartContainer").style("display", "none");
    d3.select("#sequenceContainer").style("display", "block");
  }}


  function zoomToFeature(feature, transcript) {{
    currentView = "zoomed";
    d3.select("#backButtonContainer").style("display", "block");

    var seqStart = feature.start - 100;
    var seqEnd   = feature.end   + 100;
    if (seqStart < transcript.transcript_start) seqStart = transcript.transcript_start;
    if (seqEnd   > transcript.transcript_end)   seqEnd   = transcript.transcript_end;

    xScale.domain([seqStart, seqEnd]).range([100, width - marginRight]);

    const cdsStart = d3.min(transcript.CDS, d => d.start);
    const cdsEnd   = d3.max(transcript.CDS, d => d.end);

    var seqSVG = d3.select("#sequenceSvg");
    seqSVG.html("");
    seqSVG.attr("width", width).attr("height", window.innerHeight);

    var axisY = marginTop + 20;
    drawAxis(seqSVG, xScale, axisY);

    var labelText = (transcript.transcript_id === currentGeneData.canonical.transcript_id ?
                     "Transcrit canonique : " : "Transcrit alternatif : ") +
                    transcript.transcript_id +
                    " | Positions : " + seqStart + " - " + seqEnd;

    seqSVG.append("text")
          .attr("x", 10)
          .attr("y", marginTop - 20)
          .attr("font-size", 16)
          .text(labelText);

    var rowY = marginTop + 60 + 300 + 40;
    var seqGroup = seqSVG.append("g").attr("id", "seqTranscriptGroup");

    var exons = transcript.exons.slice().sort(function(a, b) {{ return a.start - b.start; }});
    exons.forEach(function(ex, i) {{
      if (ex.end < seqStart || ex.start > seqEnd) return;

      ex.regions.forEach(function(region) {{
        var rStart = Math.max(region.start, seqStart);
        var rEnd   = Math.min(region.end,   seqEnd);
        if (rStart > rEnd) return;

        var px1 = xScale(rStart);
        var px2 = xScale(rEnd + 1);
        var w   = px2 - px1;

        var color = region.type === "CDS" ? "#ADD8E6" : "#FF0000";
        seqGroup.append("rect")
                .attr("x", px1)
                .attr("y", rowY + (rowHeight - exonHeight)/2)
                .attr("width",  w)
                .attr("height", exonHeight)
                .attr("fill",   color);

        if (w > 12) {{
          let label = region.type;
          if (region.type === "UTR") {{
            const fivePrime = (transcript.strand === "+" ?
                               (region.end < cdsStart) :
                               (region.start > cdsEnd));
            label = fivePrime ? "5'UTR" : "3'UTR";
          }}
          var fSize = Math.max(8, Math.min(14, w / label.length));
          seqGroup.append("text")
                  .attr("x", px1 + w / 2)
                  .attr("y", rowY + (rowHeight - exonHeight)/2 + exonHeight/2 + 4)
                  .attr("text-anchor", "middle")
                  .attr("font-size", fSize + "px")
                  .attr("fill", "white")
                  .text(label);
        }}
      }});

      /* intron + marqueurs inchangés */
      if (i < exons.length - 1) {{
        var nextEx = exons[i+1];
        var intronStart = ex.end + 1;
        var intronEnd   = nextEx.start - 1;
        var cStart = Math.max(intronStart, seqStart);
        var cEnd   = Math.min(intronEnd,   seqEnd);
        if (cEnd > cStart) {{
          seqGroup.append("line")
                  .attr("x1", xScale(cStart))
                  .attr("x2", xScale(cEnd))
                  .attr("y1", rowY + rowHeight/2)
                  .attr("y2", rowY + rowHeight/2)
                  .attr("stroke", "black").attr("stroke-width", 2);
        }}
        drawMarkerAtBoundary("right", ex.end, transcript, rowY, seqGroup, seqSVG);
        drawMarkerAtBoundary("left",  nextEx.start, transcript, rowY, seqGroup, seqSVG);
      }}
    }});

    lastZoomState = {{
      type: "feature",
      transcript: transcript,
      seqStart: seqStart,
      seqEnd: seqEnd
    }};

    d3.select("#chartContainer").style("display", "none");
    d3.select("#sequenceContainer").style("display", "block");
  }}

  function renderZoomedExonView(state) {{
    currentView = "zoomed";
    d3.select("#backButtonContainer").style("display", "block");

    xScale.domain([state.seqStart, state.seqEnd]).range([100, width - marginRight]);

    const cdsStart = d3.min(state.transcript.CDS, d => d.start);
    const cdsEnd   = d3.max(state.transcript.CDS, d => d.end);

    var seqSVG = d3.select("#sequenceSvg");
    seqSVG.html("");
    seqSVG.attr("width", width).attr("height", window.innerHeight);

    var axisY = marginTop + 20;
    drawAxis(seqSVG, xScale, axisY);

    var transcript = state.transcript;
    var labelText = (transcript.transcript_id === currentGeneData.canonical.transcript_id ?
                     "Transcrit canonique : " : "Transcrit alternatif : ") +
                    transcript.transcript_id +
                    " | Positions : " + state.seqStart + " - " + state.seqEnd;
    seqSVG.append("text")
          .attr("x", 10)
          .attr("y", marginTop - 20)
          .attr("font-size", "16px")
          .text(labelText);

    var histMarginTop = marginTop + 60;
    var histHeight    = 300;
    drawCoverageHistogram(state.seqStart, state.seqEnd, histMarginTop, histHeight);

    var seqGroup = seqSVG.append("g").attr("id", "seqTranscriptGroup");

    // ─── Exons/UTR ───
    state.exon.regions.forEach(function(region) {{
        var rStart = Math.max(region.start, state.seqStart);
        var rEnd   = Math.min(region.end,   state.seqEnd);
        if (rStart > rEnd) return;
        var px1 = xScale(rStart);
        var px2 = xScale(rEnd + 1);
        var w   = px2 - px1;
        var color = region.type === "CDS" ? "#ADD8E6" : "#FF0000";
        seqGroup.append("rect")
                .attr("x", px1)
                .attr("y", marginTop + histMarginTop + histHeight + 40 - rowHeight/2)
                .attr("width",  w)
                .attr("height", exonHeight)
                .attr("fill",   color);

        if (w > 12) {{
            let label = region.type;
            if (region.type === "UTR") {{
                const fivePrime = (transcript.strand === "+" ?
                                   (region.end < cdsStart) :
                                   (region.start > cdsEnd));
                label = fivePrime ? "5'UTR" : "3'UTR";
            }}
            var fSize = Math.max(8, Math.min(14, w / label.length));
            seqGroup.append("text")
                    .attr("x", px1 + w/2)
                    .attr("y", marginTop + histMarginTop + histHeight + 40 - rowHeight/2 + exonHeight/2 + 4)
                    .attr("text-anchor", "middle")
                    .attr("font-size", fSize + "px")
                    .attr("fill", region.type === "CDS" ? "black" : "white")
                    .text(label);
        }}
    }});

    // ─── Introns + marqueurs 10/20/50 ───
    var baseY = marginTop + histMarginTop + histHeight + 40 - exonHeight/2;
    if (state.exon.start > transcript.transcript_start) {{
        seqGroup.append("line")
                .attr("x1", xScale(state.seqStart))
                .attr("x2", xScale(state.exon.start))
                .attr("y1", baseY)
                .attr("y2", baseY)
                .attr("stroke", "black")
                .attr("stroke-width", 2);
        drawMarkerAtBoundary("left", state.exon.start, transcript, baseY, seqGroup, seqSVG);
    }}
    if (state.exon.end < transcript.transcript_end) {{
        seqGroup.append("line")
                .attr("x1", xScale(state.exon.end + 1))
                .attr("x2", xScale(state.seqEnd))
                .attr("y1", baseY)
                .attr("y2", baseY)
                .attr("stroke", "black")
                .attr("stroke-width", 2);
        drawMarkerAtBoundary("right", state.exon.end, transcript, baseY, seqGroup, seqSVG);
    }}

    // ─── Affichage des variants (flèches violettes) ───
    var arrowY = marginTop + histMarginTop + histHeight + 40 - rowHeight/2 - 10;
    currentGeneData.variants
        .filter(v => v.pos >= state.seqStart && v.pos <= state.seqEnd)
        .forEach(function(v) {{
            var vx = xScale(v.pos) - 5;
            seqGroup.append("path")
                    .attr("d", "M0,0 L10,0 L5,10 Z")
                    .attr("transform", "translate(" + vx + "," + arrowY + ")")
                    .attr("fill", "purple");
        }});
  }}

  function renderZoomedFeatureView(state) {{
    currentView = "zoomed";
    d3.select("#backButtonContainer").style("display", "block");

    xScale.domain([state.seqStart, state.seqEnd]).range([100, width - marginRight]);

    const cdsStart = d3.min(state.transcript.CDS, d => d.start);
    const cdsEnd   = d3.max(state.transcript.CDS, d => d.end);

    var seqSVG = d3.select("#sequenceSvg");
    seqSVG.html("");
    seqSVG.attr("width", width).attr("height", window.innerHeight);

    var axisY = marginTop + 20;
    drawAxis(seqSVG, xScale, axisY);

    var transcript = state.transcript;
    var labelText = (transcript.transcript_id === currentGeneData.canonical.transcript_id ?
                     "Transcrit canonique : " : "Transcrit alternatif : ") +
                    transcript.transcript_id +
                    " | Positions : " + state.seqStart + " - " + state.seqEnd;

    seqSVG.append("text")
          .attr("x", 10)
          .attr("y", marginTop - 20)
          .attr("font-size", 16)
          .text(labelText);

    var rowY = marginTop + 60 + 300 + 40;
    var seqGroup = seqSVG.append("g").attr("id", "seqTranscriptGroup");

    var exons = transcript.exons.slice().sort(function(a, b) {{ return a.start - b.start; }});
    exons.forEach(function(ex, i) {{
      if (ex.end < state.seqStart || ex.start > state.seqEnd) return;

      ex.regions.forEach(function(region) {{
        var rStart = Math.max(region.start, state.seqStart);
        var rEnd   = Math.min(region.end,   state.seqEnd);
        if (rStart > rEnd) return;

        var px1 = xScale(rStart);
        var px2 = xScale(rEnd + 1);
        var w   = px2 - px1;

        var color = region.type === "CDS" ? "#ADD8E6" : "#FF0000";
        seqGroup.append("rect")
                .attr("x", px1)
                .attr("y", rowY + (rowHeight - exonHeight)/2)
                .attr("width",  w)
                .attr("height", exonHeight)
                .attr("fill",   color);

        if (w > 12) {{
          let label = region.type;
          if (region.type === "UTR") {{
            const fivePrime = (transcript.strand === "+" ?
                               (region.end < cdsStart) :
                               (region.start > cdsEnd));
            label = fivePrime ? "5'UTR" : "3'UTR";
          }}
          var fSize = Math.max(8, Math.min(14, w / label.length));
          seqGroup.append("text")
                  .attr("x", px1 + w / 2)
                  .attr("y", rowY + (rowHeight - exonHeight)/2 + exonHeight/2 + 4)
                  .attr("text-anchor", "middle")
                  .attr("font-size", fSize + "px")
                  .attr("fill", "white")
                  .text(label);
        }}
      }});

      if (i < exons.length - 1) {{
        var nextEx = exons[i+1];
        var intronStart = ex.end + 1;
        var intronEnd   = nextEx.start - 1;
        var cStart = Math.max(intronStart, state.seqStart);
        var cEnd   = Math.min(intronEnd,   state.seqEnd);
        if (cEnd > cStart) {{
          seqGroup.append("line")
                  .attr("x1", xScale(cStart))
                  .attr("x2", xScale(cEnd))
                  .attr("y1", rowY + rowHeight/2)
                  .attr("y2", rowY + rowHeight/2)
                  .attr("stroke", "black").attr("stroke-width", 2);
        }}
        drawMarkerAtBoundary("right", ex.end,      transcript, rowY, seqGroup, seqSVG);
        drawMarkerAtBoundary("left",  nextEx.start, transcript, rowY, seqGroup, seqSVG);
      }}
    }});
  }}

  function drawMarkerAtBoundary(side, boundary, transcript, rowY, seqGroup, seqSVG) {{
    var distances = [10, 20, 50];
    distances.forEach(function(d) {{
      var markerX;
      if (side === "left") {{
        if (boundary - d < transcript.transcript_start) return;
        markerX = xScale(boundary - d);
      }} else {{
        if (boundary + d > transcript.transcript_end) return;
        markerX = xScale(boundary + d);
      }}
      var markerY1, markerY2, labelY;
      if (currentView === "detail" && window.squaresY !== undefined) {{
        // Vue détail
        markerY1 = window.squaresY - 50 - (exonHeight / 2);
        markerY2 = window.squaresY + 15 - (exonHeight / 2);
      }} else {{
        // Vue zoomée normale
        markerY1 = rowY + rowHeight/2 - 50 - (exonHeight / 2);
        markerY2 = rowY + rowHeight/2 + 25 - (exonHeight / 2);
      }}
      // On place le label 5px en dessous de la fin des pointillés
      labelY = markerY2 + 5;

      seqGroup.append("line")
              .attr("x1", markerX)
              .attr("x2", markerX)
              .attr("y1", markerY1)
              .attr("y2", markerY2)
              .attr("stroke", "red")
              .attr("stroke-dasharray", "4,4");
      var label = "";
      if (side === "left") {{
        label = (transcript.strand === "+" ? "-" : "+") + d;
      }} else {{
        label = (transcript.strand === "+" ? "+" : "-") + d;
      }}
      seqSVG.append("text")
            .attr("x", markerX)
            .attr("y", labelY)
            .attr("text-anchor", "middle")
            .attr("font-size", "10px")
            .attr("fill", "black")
            .text(label);
    }});
  }}

  function drawCoverageHistogram(seqStart, seqEnd, histMarginTop, histHeight) {{
      var profile = currentGeneData.coverage_profile;
      if (!profile) return;

      var filtered = profile.filter(function(d) {{ return d.pos >= seqStart && d.pos <= seqEnd; }});
      if (filtered.length === 0) return;

      var maxDepth = d3.max(filtered, function(d) {{ return d.depth; }});
      // Étendre l’axe Y de 15px vers le bas
      var yScale = d3.scaleLinear()
          .domain([0, maxDepth])
          .range([histMarginTop + histHeight + 35,   // +15px ici (était +20)
                  histMarginTop - 20]);

      var seqSVG = d3.select("#sequenceSvg");

      /* barres */
      var lowCoverageData = filtered.filter(function(d) {{ return d.depth < 50; }})
                                    .sort(function(a, b) {{ return a.pos - b.pos; }});
      filtered.forEach(function(d) {{
          var barX = xScale(d.pos);
          var barW = xScale(d.pos + 1) - xScale(d.pos);
          var barY = yScale(d.depth);
          // Étendre la hauteur des barres de 15px vers le bas
          var barH = (histMarginTop + histHeight + 35) - barY;
          var color = (d.depth >= 50) ? "blue" : "#FFA500";

          var rect = seqSVG.append("rect")
              .attr("x", barX)
              .attr("y", barY)
              .attr("width", barW)
              .attr("height", barH)
              .attr("fill", color);

          if (d.depth < 50) {{
              rect.style("cursor", "pointer")
                  .on("click", function(event) {{
                      var block = findContiguousBlock(d.pos, lowCoverageData);
                      zoomToLowCoverageDetail(block.start, block.end, currentTranscript);
                      event.stopPropagation();
                  }});
          }}
      }});

      /* axe Y */
      seqSVG.append("g")
          .attr("class", "y axis")
          .attr("transform", "translate(50,0)")
          .call(d3.axisLeft(yScale).ticks(5));

      /* seuil 50 (inchangé) */  
      seqSVG.append("text")
          .attr("x", 40)
          .attr("y", yScale(50) + 3)
          .attr("text-anchor", "end")
          .attr("font-size", "10px")
          .attr("fill", "red")
          .attr("font-weight", "bold")
          .text("50");
      seqSVG.append("line")
          .attr("x1", 45).attr("x2", 50)
          .attr("y1", yScale(50)).attr("y2", yScale(50))
          .attr("stroke", "red").attr("stroke-width", 2);
      var scaleRange = xScale.range();
      seqSVG.append("line")
          .attr("x1", scaleRange[0]).attr("x2", scaleRange[1])
          .attr("y1", yScale(50)).attr("y2", yScale(50))
          .attr("stroke", "red").attr("stroke-dasharray", "4,4");

      /* légende axe Y : recentrée sur la nouvelle hauteur */
      seqSVG.append("text")
          .attr("class", "cov-label")
          .attr("transform", "rotate(-90)")
          .attr("x", -(histMarginTop + histHeight / 2 + 7.5))  // +7.5px pour recentrage
          .attr("y", 15)
          .attr("text-anchor", "middle")
          .attr("font-size", "10px")
          .attr("fill", "black")
          .text("Profondeur de couverture (X)");
  }}
  
  renderGlobal();
</script>

</body>
</html>
"""
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html)

###############################################################################
# 3) Main : Traitement des échantillons, regroupement par gène, génération du HTML
###############################################################################
def main():
    import sys, os, re
    from collections import defaultdict

    if len(sys.argv) != 6:
        sys.stderr.write(
            "Usage: python coverage2html_v3.py "
            "<GTF_FILE> <FASTA_FILE> <VARIANTS_FILE> "
            "<SAMPLES_DIR> <OUTPUT_DIR>\n"
        )
        sys.exit(1)

    gtf_file      = sys.argv[1]
    fasta_file    = sys.argv[2]
    variants_file = sys.argv[3]
    samples_dir   = sys.argv[4]
    output_dir    = sys.argv[5]

    # Lecture unique des variants
    all_variants = read_variants_file(variants_file)

    # Liste des samples
    sample_names = [
        d for d in os.listdir(samples_dir)
        if os.path.isdir(os.path.join(samples_dir, d))
    ]

    for sample in sample_names:
        sample_path = os.path.join(samples_dir, sample)

        # 1) Collecte des fichiers *_depth.txt et des gènes associés
        gene_files = defaultdict(list)
        for fname in os.listdir(sample_path):
            if not fname.lower().endswith(".txt"):
                continue
            m = re.match(r'^([^_]+)_([^_]+)', fname)
            if not m:
                continue
            gene = m.group(2)
            gene_files[gene].append(os.path.join(sample_path, fname))

        if not gene_files:
            sys.stderr.write(f"Aucun fichier de couverture pour {sample}\n")
            continue

        genes = sorted(gene_files.keys())

        # 2) PARSING DU GTF EN UNE FOIS pour tous les gènes de cet échantillon
        all_transcripts = parse_gtf_pour_plusieurs_genes(gtf_file, genes)

        genes_data = {}
        for gene in genes:
            transcripts = all_transcripts.get(gene, {})
            if not transcripts:
                sys.stderr.write(f"Aucun transcrit pour {gene} (sample {sample})\n")
                continue

            canonical_tid = select_canonical_transcript(transcripts)
            if not canonical_tid:
                sys.stderr.write(f"Aucun canonique pour {gene} (sample {sample})\n")
                continue
            canonical_data = transcripts[canonical_tid]

            # TRI des exons/CDS et lecture FASTA (identique)
            for tid, tr in transcripts.items():
                tr['exons'].sort(key=lambda x: x['start'])
                tr['CDS'].sort(key=lambda x: x['start'])
                seq = read_fasta_region(
                    fasta_file,
                    tr['original_chrom'],
                    tr['transcript_start'],
                    tr['transcript_end'],
                    tr['strand']
                )
                tr['full_sequence'] = seq

            # Construire alt_transcripts
            alt_transcripts = []
            for tid, tr in transcripts.items():
                if tid == canonical_tid:
                    continue
                alt_transcripts.append({
                    'transcript_id':   tid,
                    'chrom':           tr['chrom'],
                    'strand':          tr['strand'],
                    'exons':           tr['exons'],
                    'CDS':             tr['CDS'],
                    'transcript_start':tr['transcript_start'],
                    'transcript_end':  tr['transcript_end'],
                    'cds_length':      tr.get('cds_length', 0),
                    'full_sequence':   tr.get('full_sequence', '')
                })

            # Calcul des intervalles low‑coverage et du profil
            all_intervals = []
            coverage_dict = {}
            for cov_file in gene_files[gene]:
                intervals = read_coverage_intervals(cov_file)
                all_intervals.extend(intervals)
                with open(cov_file, 'r') as cf:
                    for line in cf:
                        cols = re.split(r'\s+', line.strip())
                        if len(cols) < 3:
                            continue
                        try:
                            pos   = int(cols[1])
                            depth = float(cols[2])
                        except ValueError:
                            continue
                        coverage_dict[pos] = max(depth, coverage_dict.get(pos, 0.0))

            merged = merge_intervals(all_intervals)
            coverage_intervals = [{'start':s,'end':e} for s,e in merged]
            t0 = canonical_data['transcript_start']
            t1 = canonical_data['transcript_end']
            coverage_profile = [
                {'pos':p,'depth':coverage_dict[p]}
                for p in sorted(coverage_dict)
                if t0 <= p <= t1
            ]

            # Filtrer les variants pour ce gène
            gene_variants = [
                v for v in all_variants
                if v['chrom'] == canonical_data['chrom']
                   and t0 <= v['pos'] <= t1
            ]

            genes_data[gene] = {
                'canonical':       { **canonical_data },
                'transcripts':     alt_transcripts,
                'coverage':        {'intervals': coverage_intervals},
                'coverage_profile':coverage_profile,
                'variants':        gene_variants
            }

        if not genes_data:
            sys.stderr.write(f"Aucun gène valide pour {sample}\n")
            continue

        out_html = os.path.join(output_dir, f"{sample}.html")
        generate_html(genes_data, out_html)
        print(f"HTML généré pour {sample} -> {out_html}")

