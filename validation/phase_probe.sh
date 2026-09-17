#!/bin/bash
# phase_probe.sh - Determine haplotype phase (cis/trans) of the two pathogenic
# BUB1B alleles (chr15:40209701 T>G p.Leu737Ter and chr15:40220612 T>G p.Asn1002Lys)
# in the MVA challenge WGS.  Run on dichtator inside the `phase` conda env.
#
# Result: MEASURED-UNRESOLVABLE (physical limit) - see PHASE_PROBE.md
#
# Strategy: (1) triage reads mapping to chr15, keep BUB1B-window candidates;
# (2) re-align candidates to the WHOLE GRCh38 genome for honest MAPQ;
# (3) whatshap phase the window; (4) analyse the phase-block structure.
set -euo pipefail
export PATH="$HOME/.conda/envs/phase/bin:$PATH"
cd "$HOME/hf-rare/phase_probe"

DATA="$HOME/hf-rare/data"
REF="$HOME/hf-rare/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
WIN_START=40150000
WIN_END=40270000
LANES="1 2 3 4"
log(){ echo "[phase_probe] $*"; }

log "STEP 0: reference window + probes (chr15:40,150,000-40,270,000)"
samtools faidx "$REF" chr15 > chr15.fa
samtools faidx "$REF" chr15:40150000-40270000 > region_chr15_40150000_40270000.fa
python3 build_probes.py region_chr15_40150000_40270000.fa probes_31mer.txt

log "STEP 1: triage - map every read to chr15, keep BUB1B-window candidates"
: > all_names.txt
for L in $LANES; do for R in 1 2; do
  FQ="$DATA/WGS_EX2312012_HGWCNDSX7_S16_L00${L}_R${R}_001.fastq.gz"
  if [ ! -s "names_L00${L}_R${R}.txt" ]; then
    pigz -p 4 -dc "$FQ" | minimap2 -x sr -a -t 12 chr15.fa - 2>/dev/null \
      | python3 filter_sam_window.py > "names_L00${L}_R${R}.txt"
  fi
done; done
cat names_L00*_R*.txt | sort -u > all_names.txt
log "candidate reads: $(wc -l < all_names.txt)"

log "STEP 2: extract both mates (per-file, then concatenate; avoid interleaved append)"
for L in $LANES; do
  for R in 1 2; do
    pigz -p 4 -dc "$DATA/WGS_EX2312012_HGWCNDSX7_S16_L00${L}_R${R}_001.fastq.gz" \
      | seqtk subseq - all_names.txt > "ex_R${R}_L00${L}.fq"
  done
done
cat ex_R1_L00{1,2,3,4}.fq > cand_R1.fq
cat ex_R2_L00{1,2,3,4}.fq > cand_R2.fq
log "R1=$(grep -c ^@ cand_R1.fq) R2=$(grep -c ^@ cand_R2.fq)"

log "STEP 3: realign candidates to WHOLE GRCh38 (honest MAPQ), build BAM"
cat cand_R1.fq cand_R2.fq > cand_all.fq
minimap2 -x sr -a -t 24 "$REF" cand_all.fq 2>/dev/null \
  | samtools view -b -F 0x904 - 2>/dev/null > cand.unsorted.bam
samtools sort -@ 16 cand.unsorted.bam -o cand.bam
samtools index cand.bam
log "BAM: $(samtools view -c cand.bam) primary reads; window=$(samtools view -c cand.bam chr15:$WIN_START-$WIN_END)"

log "STEP 4: prepare VCF (chr-prefixed contig names to match GRCh38 ref)"
bcftools view "$DATA/WGS_EX2312012_HGWCNDSX7.vcf.gz" 15:40100000-40230000 > bub1b_region.vcf
bcftools annotate --rename-chrs chr_map.txt bub1b_region.vcf > bub1b_region_chr.vcf
bgzip -f bub1b_region_chr.vcf && tabix -p vcf bub1b_region_chr.vcf.gz

log "STEP 5: whatshap phase (BUB1B window)"
whatshap phase --reference "$REF" --ignore-read-groups --sample WGS_EX2312012 \
  -o phased_bub1b.vcf.gz bub1b_region_chr.vcf.gz cand.bam --chromosome chr15
bcftools view -H phased_bub1b.vcf.gz > phased_bub1b.rows.txt

log "STEP 6: analysis (blocks, target genotypes, AD sanity, insert size)"
python3 block_table.py phased_bub1b.rows.txt
python3 analyse_phase.py phased_bub1b.rows.txt cand.bam
python3 interleave.py cand_R1.fq cand_R2.fq cand_interleaved.fq
minimap2 -x sr -a -t 16 chr15.fa cand_interleaved.fq 2>/dev/null | python3 measure_insert.py

log "DONE -> PHASE_PROBE.md"
