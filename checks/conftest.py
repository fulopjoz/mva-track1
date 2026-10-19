"""Shared synthetic-reference fixtures for Track-1 annotator/pipeline tests.

All coordinates, genes and variants are invented; no case or clinical data.
The mini GTF builds three transcripts:
  GENEP (chr1, + strand, multi-exon), exons (100-199),(300-399), cds (110-199),(300-389)
  GENEM (chr2, - strand, multi-exon), same geometry
  GENES (chr3, + strand, single exon), exon (100-299), cds (150-250)
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent.parent / 'pipeline' / 'scripts'
sys.path.insert(0, str(SCRIPT_DIR))

GTF_LINES = []
for tid, gene, chrom, strand in (
        ('TXP.1', 'GENEP', '1', '+'),
        ('TXM.1', 'GENEM', '2', '-')):
    for feat, spans in (('exon', ((100, 199), (300, 399))),
                        ('CDS', ((110, 199), (300, 389)))):
        for s, e in spans:
            GTF_LINES.append(
                '\t'.join([chrom, 'syn', feat, str(s), str(e), '.', strand, '.',
                           f'gene_name "{gene}"; transcript_id "{tid}"']))
# single-exon transcript
for feat, s, e in (('exon', 100, 299), ('CDS', 150, 250)):
    GTF_LINES.append('\t'.join(['3', 'syn', feat, str(s), str(e), '.', '+', '.',
                                'gene_name "GENES"; transcript_id "TXS.1"']))
GTF_TEXT = '\n'.join(GTF_LINES) + '\n'


def write_gtf(tmp_path: Path) -> Path:
    p = tmp_path / 'mini.gtf.gz'
    with gzip.open(p, 'wt') as fh:
        fh.write(GTF_TEXT)
    return p


VCF_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')


def write_vcf(tmp_path: Path, rows: list[str]) -> Path:
    p = tmp_path / 'fixture.vcf'
    p.write_text(VCF_HEADER + '\n'.join(rows) + '\n')
    return p


def write_conf(tmp_path: Path, gtf: Path, vcf: Path, out_dir: Path | None = None) -> Path:
    conf = {
        'data': {'vcf': str(vcf), 'vcf_chr_prefix': False, 'sample': 'S1', 'build': 'GRCh38'},
        'refs': {
            'cadd_dir': [str(tmp_path / 'no_cadd')],
            'mane_gtf': str(gtf),
            'genome_fasta': str(tmp_path / 'absent.fna'),
            'hp_obo': str(tmp_path / 'absent.obo'),
            'hpo_map': str(tmp_path / 'absent.tsv'),
            'gnomad_v41_exomes_dir': str(tmp_path),
            'clinvar_vcf': str(tmp_path / 'absent_clinvar.vcf.gz'),
        },
        'variants': {
            'qc': {'min_gq': 20, 'min_dp': 10, 'min_alt_frac': 0.25},
            'af': {'max_rare': 0.01, 'max_ultrarare': 0.0001},
            'effect': {
                'lof_terms': ['stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor'],
                'missense_terms': ['missense_variant'],
                'inframe_terms': ['inframe_deletion', 'inframe_insertion'],
                'exploit_cadd_if_available': True,
            },
            'pairing': {'compound_het_rule': 'same_gene', 'same_gene_only': True,
                        'allow_unknown_phase': True},
        },
        'phenotype': {'hpo_terms': []},
        'scoring': {'hpo_weight': 2.0, 'founder_weight': 0.45,
                    'pathogenicity_weight': 0.60, 'rare_weight': 0.30},
        'output': {'write_dir': str(out_dir or tmp_path),
                   'track1_csv_columns': ['proband_id']},
    }
    p = tmp_path / 'conf.yml'
    yaml.safe_dump(conf, open(p, 'w'))
    return p


@pytest.fixture()
def mini_refs(tmp_path):
    gtf = write_gtf(tmp_path)
    return {'gtf': gtf, 'tmp': tmp_path}
