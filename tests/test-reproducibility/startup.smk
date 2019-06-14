REFERENCE = 'reference/hg38_chr20.fasta'
RRBS_SE = 'data/se.fastq.gz'
RRBS_PE = ['data/pe.read1.fastq.gz', 'data/pe.read2.fastq.gz']

ALL = []
ALL.append(REFERENCE)
ALL.append(RRBS_SE)
ALL.append(RRBS_PE)

rule all:
    input: ALL

rule clean:
    shell:
        "if [ -d data ]; then rm -r data; fi; "
        "if [ -d reference ]; then rm -r reference; fi; "
        "if [ -d logs ]; then rm -r logs; fi; "
        "if [ -d benchmarks ]; then rm -r benchmarks; fi; "
        "if [ -d result0 ]; then rm -r result0; fi; "
        "if [ -d result1 ]; then rm -r result1; fi; "

rule reference:
    output: REFERENCE
    wrapper: 'http://dohlee-bio.info:9193/test/reference'

rule rrbs_se:
    output: RRBS_SE
    wrapper: 'http://dohlee-bio.info:9193/test/rrbs/se'

rule rrbs_pe:
    output: RRBS_PE
    wrapper: 'http://dohlee-bio.info:9193/test/rrbs/pe'
