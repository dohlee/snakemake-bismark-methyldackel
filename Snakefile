import pandas as pd
from pathlib import Path

configfile: 'config.yaml'

include: 'rules/trim-galore.smk'
include: 'rules/bismark.smk'
include: 'rules/methyldackel.smk'
include: 'rules/sambamba.smk'

manifest = pd.read_csv(config['manifest'])
DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir'])

SE_SAMPLES = manifest[manifest.library_layout == 'single'].name.values
PE_SAMPLES = manifest[manifest.library_layout == 'paired'].name.values

SE_BEDGRAPHS = expand(str(RESULT_DIR / '03_methyldackel' / 'se' / '{sample}_CpG.bedGraph'), sample=SE_SAMPLES)
PE_BEDGRAPHS = expand(str(RESULT_DIR / '03_methyldackel' / 'pe' / '{sample}_CpG.bedGraph'), sample=PE_SAMPLES)

ALL = []
ALL.append(SE_BEDGRAPHS)
ALL.append(PE_BEDGRAPHS)

rule all:
    input: ALL

rule clean:
    shell:
        "if [ -d {RESULT_DIR} ]; then rm -r {RESULT_DIR}; fi; "
        "if [ -d {DATA_DIR} ]; then rm -r {DATA_DIR}; fi; "
        "if [ -d logs ]; then rm -r logs; fi; "
        "if [ -d benchmarks ]; then rm -r benchmarks; fi; "
