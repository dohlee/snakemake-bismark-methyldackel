import pandas as pd
from pathlib import Path

configfile: 'config.yaml'

include: 'rules/trim-galore.smk'
include: 'rules/bismark.smk'
include: 'rules/fastqc.smk'
include: 'rules/methyldackel.smk'
include: 'rules/sambamba.smk'

ruleorder: trim_galore_pe > trim_galore_se

manifest = pd.read_csv(config['manifest'])
DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir'])

SAMPLES = manifest.name.values
SAMPLE2LIB = {r.name:r.library_layout for r in manifest.to_records()}
SE_SAMPLES = manifest[manifest.library_layout == 'single'].name.values
PE_SAMPLES = manifest[manifest.library_layout == 'paired'].name.values

RAW_QC_SE = expand(str(DATA_DIR / '{sample}_fastqc.zip'), sample=SE_SAMPLES)
RAW_QC_PE = expand(str(DATA_DIR / '{sample}.read1_fastqc.zip'), sample=PE_SAMPLES)
TRIMMED_QC_SE = expand(str(RESULT_DIR / '01_trim-galore' / '{sample}.trimmed_fastqc.zip'), sample=SE_SAMPLES)
TRIMMED_QC_PE = expand(str(RESULT_DIR / '01_trim-galore' / '{sample}.read1.trimmed_fastqc.zip'), sample=PE_SAMPLES)
BEDGRAPHS = expand(str(RESULT_DIR / '03_methyldackel' / '{sample}_CpG.bedGraph'), sample=SAMPLES)

ALL = []
ALL.append(BEDGRAPHS)
ALL.append(RAW_QC_SE)
ALL.append(RAW_QC_PE)
ALL.append(TRIMMED_QC_SE)
ALL.append(TRIMMED_QC_PE)

rule all:
    input: ALL

rule clean:
    shell:
        "if [ -d {RESULT_DIR} ]; then rm -r {RESULT_DIR}; fi; "
        "if [ -d {DATA_DIR} ]; then rm -r {DATA_DIR}; fi; "
        "if [ -d logs ]; then rm -r logs; fi; "
        "if [ -d benchmarks ]; then rm -r benchmarks; fi; "
