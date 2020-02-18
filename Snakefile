import pandas as pd
from pathlib import Path

configfile: 'config.yaml'

include: 'rules/trim-galore.smk'
include: 'rules/bismark.smk'
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

BEDGRAPHS = expand(str(RESULT_DIR / '03_methyldackel' / '{sample}_CpG.bedGraph'), sample=SAMPLES)

ALL = []
ALL.append(BEDGRAPHS)

rule all:
    input: ALL

rule clean:
    shell:
        "if [ -d {RESULT_DIR} ]; then rm -r {RESULT_DIR}; fi; "
        "if [ -d {DATA_DIR} ]; then rm -r {DATA_DIR}; fi; "
        "if [ -d logs ]; then rm -r logs; fi; "
        "if [ -d benchmarks ]; then rm -r benchmarks; fi; "
