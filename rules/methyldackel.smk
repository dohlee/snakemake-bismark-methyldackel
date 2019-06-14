from pathlib import Path
RESULT_DIR = Path(config['result_dir'])

MES = config['methyldackel_extract_se']
rule methyldackel_extract_se:
    input:
        bam = RESULT_DIR / '02_bismark' / 'se' / '{sample}' / '{sample}.bismark.sorted.bam',
        bam_index = RESULT_DIR / '02_bismark' / 'se' / '{sample}' / '{sample}.bismark.sorted.bam.bai',
        reference = config['reference']['fasta'],
    output:
        # NOTE: {sample}_CpG.meth.bedGraph will automatically use --fraction option,
        # {sample}_CpG.counts.bedGraph will automatically use --counts option, and
        # {sample}_CpG.logit.bedGraph will automatically use --logit option.
        RESULT_DIR / '03_methyldackel' / 'se' / '{sample}_CpG.bedGraph',
    params:
        # Optional parameters. Omit if unneeded.
        extra = MES['extra'],
        # Minimum MAPQ threshold to include an alignment.
        # Default: 10
        q = MES['q'],
        # Minimum Phred threshold to include a base. This must be >0.
        # Default: 5
        p = MES['p'],
        # Maximum per-base depth.
        # Default: 2000
        D = MES['D'],
        # Minimum per-base depth for reporting output. If you use --mergeContext, this then applies to the
        # merged CpG/CHG.
        # Default: 1
        d = MES['d'],
        # Region string in which to extract methylation.
        # Default: False
        r = MES['r'],
        # A BED file listing regions for inclusion.
        # Default: False
        l = MES['l'],
        # If a BED file is specified, then this option will cause the strand column (column 6) to be utilized,
        # if present. Thus, if a region has a '+' in this column, then only metrics from the top strand will 
        # be output. Note that the -r option can be used to limit the regions of -l.
        keepStrand = MES['keepStrand'],
        # The size of the genome processed by a single thread at a time.
        # The default is 1000000 bases. This value MUST be at least 1.
        # Default: 1000000
        chunkSize = MES['chunkSize'],
        # Merge per-Cytosine metrics from CpG and CHG contexts into per-CPG or per-CHG metrics.
        # Default: False,
        mergeContext = MES['mergeContext'],
        # By default, any alignment marked as a duplicate is ignored. This option causes them to be incorporated.
        # This will unset bit 0x400 in --ignoreFlags.
        # Default: False
        keepDupes = MES['keepDupes'],
        # By default, if only one read in a pair aligns (a singleton) then it's ignored.
        # Default: False
        keepSingleton = MES['keepSingleton'],
        # By default, paired-end alignments with the properly-paired bit unset in the FLAG field are ignored.
        # Note that the definition of concordant and discordant is based on your aligner settings.
        # Default: False
        keepDiscordant = MES['keepDiscordant'],
        # By default, any alignment marked as secondary (bit 0x100), failing QC (bit 0x200), a PCR/optical
        # duplicate (0x400) or supplemental (0x800) is ignored. This equates to a value of 0xF00 or 3840 in
        # decimal. If you would like to change that, you can specify a new value here.
        # Default: False
        ignoreFlags = MES['ignoreFlags'],
        # Require each alignment to have all bits in this value present, or else the alignment is ignored.
        # This is equivalent to the -f option in samtools. The default is 0, which includes all alignments.
        # Default: False
        requireFlags = MES['requireFlags'],
        # Do not output CpG context methylation metrics.
        # Default: False
        noCpG = MES['noCpG'],
        # Output CHG context methylation metrics.
        # Default: False
        CHG = MES['CHG'],
        # Output CHH context methylation metrics.
        # Default: False
        CHH = MES['CHH'],
        # If you would like to exclude sites that likely contain SNPs, then specifying a value greater than 0
        # here will indicate the minimum depth required on the strand opposite of a C to look for A/T/C bases. 
        # The fraction of these necessary to exclude a position from methylation extraction is specified by 
        # the --maxVariantFrac parameter. The default is 0, which means that no positions will be excluded.
        # Note that the -p and -q apply here as well. Note further that if you use --mergeContext that a
        # merged site will be excluded if either of its individual Cs would be excluded.
        # Default: False
        minOppositeDepth = MES['minOppositeDepth'],
        # The maximum fraction of A/T/C base calls on the strand opposite of a C to allow before a position
        # is declared a variant (thereby being excluded from output). The default is 0.0.
        # See also --minOppositeDepth.
        # Default: False
        maxVariantFrac = MES['maxVariantFrac'],
        # Output in the format required by methylKit. Note that this is incomaptible with --mergeContext,
        # --fraction and --counts.
        # Default: False
        methylKit = MES['methylKit'],
        # A per-base exhaustive report comparable to that produced with the same option in Bismark's methylation
        # extractor. The output is a tab-separated file with the following columns:
        # chromosome, position (1-based), strand, number of alignments supporting methylation, number of
        # alignments supporting unmethylation, CG/CHG/CHH, trinucleotide context. This is not compatible with
        # --fraction, --counts, --methylKit, or --mergeContext. This produces a single file with a
        # .cytosine_report.txt extension. Note that even bases with no coverage will be included in the output.
        # Default: False
        cytosine_report = MES['cytosine_report'],
        # Inclusion bounds for methylation calls from reads/pairs originating from the original top strand.
        # Suggested values can be obtained from the MBias program. Each integer represents a 1-based position
        # on a read. For example --OT A,B,C,D translates to, "Include calls at positions from A through B on 
        # read #1 and C through D on read #2". If a 0 is used at any position then that is translated to 
        # mean star/end of the alignment, as appropriate. For example, --OT 5,0,0,0 would include all but
        # the first 4 bases on read #1. Users are strongly advised to consult a methylation bias plot,
        # for example by using the MBias program.
        # Default: False
        # Example: OT = '5,5,5,5'
        OT = MES['OT'],
        # As with --OT, but for the original bottom, complementary to the original top, and complementary to
        # the original bottom strands, respectively.
        # Default: False for all
        OB = MES['OB'],
        CTOT = MES['CTOT'],
        CTOB = MES['CTOB'],
        # Like --OT, but always exclude INT bases from a given end from inclusion, regardless of the length of
        # an alignment. This is useful in cases where reads may have already been trimmed to different
        # lengths, but still non-the-less contain a certain length bias at one or more ends.
        # Default: False
        nOT = MES['nOT'],
        # As with --nOT, but for the original bottom, complementary to the original top, and complementary to
        # the original bottom strands, respectively.
        nOB = MES['nOB'],
        nCTOT = MES['nCTOT'],
        nCTOB = MES['nCTOB'],
    threads: config['threads']['methyldackel_extract_se']
    log: 'logs/methyldackel_extract/{sample}.log'
    benchmark: 'benchmarks/methyldackel_extract/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/methyldackel/extract'

MEP = config['methyldackel_extract_pe']
rule methyldackel_extract_pe:
    input:
        bam = RESULT_DIR / '02_bismark' / 'pe' / '{sample}' / '{sample}.bismark.sorted.bam',
        bam_index = RESULT_DIR / '02_bismark' / 'pe' / '{sample}' / '{sample}.bismark.sorted.bam.bai',
        reference = config['reference']['fasta'],
    output:
        # NOTE: {sample}_CpG.meth.bedGraph will automatically use --fraction option,
        # {sample}_CpG.counts.bedGraph will automatically use --counts option, and
        # {sample}_CpG.logit.bedGraph will automatically use --logit option.
        RESULT_DIR / '03_methyldackel' / 'pe' / '{sample}_CpG.bedGraph',
    params:
        # Optional parameters. Omit if unneeded.
        extra = MEP['extra'],
        # Minimum MAPQ threshold to include an alignment.
        # Default: 10
        q = MEP['q'],
        # Minimum Phred threshold to include a base. This must be >0.
        # Default: 5
        p = MEP['p'],
        # Maximum per-base depth.
        # Default: 2000
        D = MEP['D'],
        # Minimum per-base depth for reporting output. If you use --mergeContext, this then applies to the
        # merged CpG/CHG.
        # Default: 1
        d = MEP['d'],
        # Region string in which to extract methylation.
        # Default: False
        r = MEP['r'],
        # A BED file listing regions for inclusion.
        # Default: False
        l = MEP['l'],
        # If a BED file is specified, then this option will cause the strand column (column 6) to be utilized,
        # if present. Thus, if a region has a '+' in this column, then only metrics from the top strand will 
        # be output. Note that the -r option can be used to limit the regions of -l.
        keepStrand = MEP['keepStrand'],
        # The size of the genome processed by a single thread at a time.
        # The default is 1000000 bases. This value MUST be at least 1.
        # Default: 1000000
        chunkSize = MEP['chunkSize'],
        # Merge per-Cytosine metrics from CpG and CHG contexts into per-CPG or per-CHG metrics.
        # Default: False,
        mergeContext = MEP['mergeContext'],
        # Output filename prefix. CpG/CHG/CHH context metrics will be output to STR_CpG.bedGraph and so on.
        # Default: False
        opref = MEP['opref'],
        # By default, any alignment marked as a duplicate is ignored. This option causes them to be incorporated.
        # This will unset bit 0x400 in --ignoreFlags.
        # Default: False
        keepDupes = MEP['keepDupes'],
        # By default, if only one read in a pair aligns (a singleton) then it's ignored.
        # Default: False
        keepSingleton = MEP['keepSingleton'],
        # By default, paired-end alignments with the properly-paired bit unset in the FLAG field are ignored.
        # Note that the definition of concordant and discordant is based on your aligner settings.
        # Default: False
        keepDiscordant = MEP['keepDiscordant'],
        # By default, any alignment marked as secondary (bit 0x100), failing QC (bit 0x200), a PCR/optical
        # duplicate (0x400) or supplemental (0x800) is ignored. This equates to a value of 0xF00 or 3840 in
        # decimal. If you would like to change that, you can specify a new value here.
        # Default: False
        ignoreFlags = MEP['ignoreFlags'],
        # Require each alignment to have all bits in this value present, or else the alignment is ignored.
        # This is equivalent to the -f option in samtools. The default is 0, which includes all alignments.
        # Default: False
        requireFlags = MEP['requireFlags'],
        # Do not output CpG context methylation metrics.
        # Default: False
        noCpG = MEP['noCpG'],
        # Output CHG context methylation metrics.
        # Default: False
        CHG = MEP['CHG'],
        # Output CHH context methylation metrics.
        # Default: False
        CHH = MEP['CHH'],
        # If you would like to exclude sites that likely contain SNPs, then specifying a value greater than 0
        # here will indicate the minimum depth required on the strand opposite of a C to look for A/T/C bases. 
        # The fraction of these necessary to exclude a position from methylation extraction is specified by 
        # the --maxVariantFrac parameter. The default is 0, which means that no positions will be excluded.
        # Note that the -p and -q apply here as well. Note further that if you use --mergeContext that a
        # merged site will be excluded if either of its individual Cs would be excluded.
        # Default: False
        minOppositeDepth = MEP['minOppositeDepth'],
        # The maximum fraction of A/T/C base calls on the strand opposite of a C to allow before a position
        # is declared a variant (thereby being excluded from output). The default is 0.0.
        # See also --minOppositeDepth.
        # Default: False
        maxVariantFrac = MEP['maxVariantFrac'],
        # Output in the format required by methylKit. Note that this is incomaptible with --mergeContext,
        # --fraction and --counts.
        # Default: False
        methylKit = MEP['methylKit'],
        # A per-base exhaustive report comparable to that produced with the same option in Bismark's methylation
        # extractor. The output is a tab-separated file with the following columns:
        # chromosome, position (1-based), strand, number of alignments supporting methylation, number of
        # alignments supporting unmethylation, CG/CHG/CHH, trinucleotide context. This is not compatible with
        # --fraction, --counts, --methylKit, or --mergeContext. This produces a single file with a
        # .cytosine_report.txt extension. Note that even bases with no coverage will be included in the output.
        # Default: False
        cytosine_report = MEP['cytosine_report'],
        # Inclusion bounds for methylation calls from reads/pairs originating from the original top strand.
        # Suggested values can be obtained from the MBias program. Each integer represents a 1-based position
        # on a read. For example --OT A,B,C,D translates to, "Include calls at positions from A through B on 
        # read #1 and C through D on read #2". If a 0 is used at any position then that is translated to 
        # mean star/end of the alignment, as appropriate. For example, --OT 5,0,0,0 would include all but
        # the first 4 bases on read #1. Users are strongly advised to consult a methylation bias plot,
        # for example by using the MBias program.
        # Default: False
        # Example: OT = '5,5,5,5'
        OT = MEP['OT'],
        # As with --OT, but for the original bottom, complementary to the original top, and complementary to
        # the original bottom strands, respectively.
        # Default: False for all
        OB = MEP['OB'],
        CTOT = MEP['CTOT'],
        CTOB = MEP['CTOB'],
        # Like --OT, but always exclude INT bases from a given end from inclusion, regardless of the length of
        # an alignment. This is useful in cases where reads may have already been trimmed to different
        # lengths, but still non-the-less contain a certain length bias at one or more ends.
        # Default: False
        nOT = MEP['nOT'],
        # As with --nOT, but for the original bottom, complementary to the original top, and complementary to
        # the original bottom strands, respectively.
        nOB = MEP['nOB'],
        nCTOT = MEP['nCTOT'],
        nCTOB = MEP['nCTOB'],
    threads: config['threads']['methyldackel_extract_pe']
    log: 'logs/methyldackel_extract/{sample}.log'
    benchmark: 'benchmarks/methyldackel_extract/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/methyldackel/extract'
