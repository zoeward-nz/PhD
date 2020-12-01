import glob, os, os.path, datetime, sys, csv
from os.path import join

## Workflow starts with trimmed data! Trimming is done separately!!! 
## ============================================================================= 
## LOAD VARIABLES FROM CONFIGFILE
## config-file needs to be submitted on command-line via --configfile

## General
LOGDIR = os.path.abspath(config["logs"])
BENCHMARKDIR = os.path.abspath(config["benchmarks"])
ENVDIR = os.path.abspath(config["envdir"])
SCRIPTS = os.path.abspath(config["scriptdir"])
DATA = os.path.abspath(config["datadir"])
RES = os.path.abspath(config["results"])
THREADS = config["threads"]


## sample inputs
SAMPLEDIR = os.path.abspath(config["samples"]["dir"])
SAMPLESHEET = os.path.abspath(config["samples"]["samplesheet"])
R1 = config["samples"]["r1"]
R2 = config["samples"]["r2"]
FASTQ = os.path.abspath(config["fastq"])

## genome
GENOME_DIR = os.path.abspath(config["genome"]["dir"])
GENOME = join(GENOME_DIR, config["genome"]["genome"])
ANNOTATION = join(GENOME_DIR, config["genome"]["gtf"])
ANNOTATIONBED = join(GENOME_DIR, config["genome"]["bed"])

## TOOLS
## STAR

#STAR_THREADSIDX = config["star"]["threads_idx"]
STAR_SJDBOVERHANG = config["star"]["sjdbOverhang"]
STAR_OUTFILTERTYPE = config["star"]["outFilterType"]
STAR_OUTSAMSTRANDFIELD = config["star"]["outSAMstrandField"]
STAR_OUTFILTERMULTIMAPNMAX = config["star"]["outFilterMultimapNmax"]
STAR_ALIGNSJOVERHANGMIN = config["star"]["alignSJoverhangMin"]
STAR_ALIGNSJDBOVERHANGMIN = config["star"]["alignSJDBoverhangMin"]
STAR_OUTFILTERMISMATCHNMAX = config["star"]["outFilterMismatchNmax"]
STAR_ALIGNINTRONMIN = config["star"]["alignIntronMin"]
STAR_ALIGNINTRONMAX = config["star"]["alignIntronMax"]
STAR_ALIGNMATESGAPMAX = config["star"]["alignMatesGapMax"]
STAR_CHIMSEGMENTMIN = config["star"]["chimSegmentMin"]
STAR_LIMITBAMSORTRAM = config["star"]["limitBAMsortRAM"]
STAR_SCOREGENOMICLENGTHLOG2SCALE = config["star"]["scoreGenomicLengthLog2scale"]

# junction selection
SJ_NUM = config["junctions"]["num_samples_with_sj"]
SJ_NUM_UM = config["junctions"]["num_uniqmappers"]

# stringtie
STRINGTIE_J = config["stringtie"]["j"]
STRINGTIE_F = config["stringtie"]["F"]
STRINGTIE_THREADSMERGE = config["stringtie"]["threads_merge"]
STRINGTIE_STRAND= config["stringtie"]["strand"]

# CPAT
CPAT_RDATA = config["cpat"]["rdata"]
CPAT_TSV = config["cpat"]["tsv"]
CPAT_CP = config["cpat"]["cp"]

# SALMON
SALMON_EXTRA = config["salmon"]["extra"]

# CIRCS
DOCIRCS = config["circs"]

# multiqc
MULTIQC_FILE = os.path.abspath(config["multiqc"]["configfile"])

## outputs
MAPPING = join(RES, config["outputs"]["mapping"])
MAPPED1 = join(MAPPING, "mapped_1stpass")
MAPPED2 = join(MAPPING, "mapped_2ndpass")
STRINGTIE = join(RES, config["outputs"]["stringtie"])
GFF = join(RES, config["outputs"]["gffcompare"])
NOVEL = join(RES, config["outputs"]["novel"])
COUNTS = join(RES, config["outputs"]["counts"])
COUNTS2 = join(RES, config["outputs"]["counts2"])
SALMON_TRANS_INDEX = join(RES, config["outputs"]["salmon_trans_index"])
CIRCEXP = join(RES, config["outputs"]["circs"])
QC = join(RES, config["outputs"]["qc"])
MULTIQC     = join(RES, config["outputs"]["multiqc"])
###-------------------------------------------------------------------------------

#SAMPLES, = glob_wildcards(join(FASTQ,'{sample}_fwd_paired.fastq'))
## reading samplename from samplesheet
print("Reading samples from samplesheet...")
SAMPLES = []
reader = csv.reader(open(SAMPLESHEET), delimiter="\t")
for a in reader:
	SAMPLES.append(a[0])
# test if sample in dir
for fname in expand(join(SAMPLEDIR, R1), sample=SAMPLES):
	if not os.path.isfile(fname):
		print("File '{}' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
		sys.exit()

NUM_SAMPLES = len(SAMPLES)
print('{} samples to process'.format(NUM_SAMPLES))



#####################################################################################
## SET UP TARGETS
#####################################################################################
TARGETS = [expand(join(COUNTS2, "{sample}"), sample=SAMPLES),
		join(NOVEL, "final.novel.cpat.transcript.list.bed.gz"),  # TX-ID workflow done
		join(STRINGTIE, 'tx_gene_map.txt'),
		expand(join(QC, "{sample}.junction.annot.done"),sample=SAMPLES),
                expand(join(QC, "{sample}.junction.sat.done"),sample=SAMPLES),
                expand(join(QC, "{sample}.infer.exp.txt"), sample=SAMPLES),
                expand(join(QC, "{sample}.read.dist.txt"), sample=SAMPLES),
		join(CIRCEXP, 'all_circular_sorted_bed.gz'),
		MULTIQC]
		#join(COUNTS2, "results-tx.txt.gz"),
		#join(COUNTS2, "results-gene.txt.gz"),
		#join(COUNTS2, "results-tx.stats.txt.gz"),
		#join(COUNTS2, "results-tx.novels.txt.gz"),
		#join(COUNTS2, "results-tx.novels.stats.txt.gz")]

# Are we doing circ identification?
#if DOCIRCS:
#	TARGETS += [join(CIRCEXP, 'all_circular_sorted_bed.gz')  # CIRC workflow done



#####################################################################################
## RULES
#####################################################################################


rule all:
	input:
		TARGETS

#####################################################################################
## get data

## first I need to create the bed file for rseqc to work
#rule gtf2bed12:
#	input:
#		ANNOTATION
#	output:
#		ANNOTATIONBED
#	params:
#		script = join(SCRIPTS, "gtf2bed.pl")
#	shell:
#		"""
#		perl {params.script} {input} > {output}
#		"""


### if the library is unstranded then we must put in --outSAMstrandField intronMotif
### if library is stranded change config file --outSAMstrandField 
rule star_mapping1:
	input:
		r1=join(FASTQ, R1),
		r2=join(FASTQ, R2)
	output:
		touch(join(MAPPED1, '{sample}/mapping.done'))
	log:
		log1=join(LOGDIR, "star_map/{sample}.stdout"),
		log2=join(LOGDIR, "star_map/{sample}.stderr")
	benchmark:
		join(BENCHMARKDIR,"{sample}_star_map.txt")
	params:
		ref=GENOME_DIR,
		gtf=ANNOTATION,
		outFileNamePrefix=join(MAPPED1,'{sample}')
	#threads: THREADS
	conda:
		join(ENVDIR, "star2.yaml")
	shell:
		"""
		nice STAR --runThreadN 8 --genomeDir {GENOME_DIR} \
		--readFilesIn {input.r1} {input.r2} \
		--sjdbGTFfile {ANNOTATION} \
		--outFileNamePrefix {params.outFileNamePrefix} \
		--outFilterType {STAR_OUTFILTERTYPE} \
		--outSAMstrandField {STAR_OUTSAMSTRANDFIELD} \
		--outFilterMultimapNmax {STAR_OUTFILTERMULTIMAPNMAX} \
		--alignSJoverhangMin {STAR_ALIGNSJOVERHANGMIN}\
		--alignSJDBoverhangMin {STAR_ALIGNSJDBOVERHANGMIN} \
		--outFilterMismatchNmax {STAR_OUTFILTERMISMATCHNMAX} \
		--alignIntronMin {STAR_ALIGNINTRONMIN} \
		--alignIntronMax {STAR_ALIGNINTRONMAX} \
		--alignMatesGapMax {STAR_ALIGNMATESGAPMAX} \
		--limitBAMsortRAM {STAR_LIMITBAMSORTRAM} \
		--scoreGenomicLengthLog2scale {STAR_SCOREGENOMICLENGTHLOG2SCALE} \
		> {log.log1} 2> {log.log2}
		"""

rule splice_junctions:
	input:
		# Here I make a target input that combines all mapping1 done files,
		expand(join(MAPPED1, '{sample}/mapping.done'), sample=SAMPLES)
	output:
		SJ=join(MAPPED1, 'SJ.out')
	log:
		log = join(LOGDIR, 'splice_junctions.log')
	benchmark:
		join(BENCHMARKDIR,"splice_junctions.txt")
	params:
		script=join(SCRIPTS, 'sjCollapseSamples.awk')
	shell:
		# Get junction stats using Alex's script
		# Remove ChrM
		# Sub-select SJs:
		# total number of uniq-mappers = SJ_NUM_UM
		# number of samples supporting the splice junction = JS_NUM
		# Only consider novel SJs, i.e. unannotated col6 ==0
		"""
		awk -f {params.script} {MAPPED1}/*SJ.out.tab | \
		sort -k1,1V -k2,2n -k3,3n | awk '{{if($7>={SJ_NUM_UM} && $10 >={SJ_NUM})print$0}}' | grep -v '^chrM' > {output} 2> {log.log}
		"""

rule star_mapping2:
	input:
		r1=join(FASTQ, R1),
		r2=join(FASTQ, R2),
		sjdb=join(MAPPED1, 'SJ.out') # this should make sure that rule "splice_junctions" has been run before
	output:
		touch(join(MAPPED2, '{sample}/mapping2.done'))
	log:
		log1=join(LOGDIR, "star_map2/{sample}.stdout"),
		log2=join(LOGDIR, "star_map2/{sample}.stderr")
	benchmark:
		join(BENCHMARKDIR,"{sample}_star_map2.txt")
	params:
		ref=GENOME_DIR,
		gtf=ANNOTATION,
		outFileNamePrefix=join(MAPPED2,'{sample}')
	#threads: THREADS
	conda:
		join(ENVDIR, "star2.yaml")
	shell:
		"""
		nice STAR --runThreadN 10 --genomeDir {GENOME_DIR} \
		--readFilesIn {input.r1} {input.r2} \
		--sjdbGTFfile {ANNOTATION} \
		--sjdbFileChrStartEnd {input.sjdb} \
		--outFileNamePrefix {params.outFileNamePrefix} \
		--outFilterType {STAR_OUTFILTERTYPE} \
		--outSAMstrandField {STAR_OUTSAMSTRANDFIELD} \
		--outFilterMultimapNmax {STAR_OUTFILTERMULTIMAPNMAX} \
		--alignSJoverhangMin {STAR_ALIGNSJOVERHANGMIN}\
		--alignSJDBoverhangMin {STAR_ALIGNSJDBOVERHANGMIN} \
		--outFilterMismatchNmax {STAR_OUTFILTERMISMATCHNMAX} \
		--alignIntronMin {STAR_ALIGNINTRONMIN} \
		--alignIntronMax {STAR_ALIGNINTRONMAX} \
		--alignMatesGapMax {STAR_ALIGNMATESGAPMAX} \
		--chimSegmentMin {STAR_CHIMSEGMENTMIN} \
		--limitBAMsortRAM {STAR_LIMITBAMSORTRAM} \
		--scoreGenomicLengthLog2scale {STAR_SCOREGENOMICLENGTHLOG2SCALE} \
		--outReadsUnmapped Fastx \
		> {log.log1} 2> {log.log2}
		"""

rule unique_bam:
	input:
		map=join(MAPPED2, '{sample}/mapping2.done')
	output:
		join(MAPPED2, "{sample}.bam")
	log:
		join(LOGDIR, "unique_bam/{sample}.log")
	benchmark:
		join(BENCHMARKDIR,"{sample}_unique_bam.txt")
	params:
		bam=join(MAPPED2, '{sample}Aligned.out.sam')
	conda:
		join(ENVDIR, "samtools.yaml")
	shell:
		"""
		nice samtools view -@ 5 -h -q 255 -b {params.bam} > {output} 2> {log}
		"""

rule bam_idx:
	input:
		join(MAPPED2, "{sample}.bam")
	output:
		join(MAPPED2, "{sample}.bai")
	log:
		join(LOGDIR, "bam_idx/{sample}.stderr")
	benchmark:
		join(BENCHMARKDIR, "bam_idx/{sample}.txt")
	threads:
		THREADS
	conda:
		join(ENVDIR, "samtools.yaml")
	params:
		params = r''	
	shell:
		"""
		nice samtools index -@ {threads} {input} {output} 2> {log}
		"""


### when the library is stranded (for our library) we put --rf into the shell command
rule stringtie:
	input:
		string_in=join(MAPPED2, '{sample}.bam')
	output:
		join(STRINGTIE, "{sample}.gtf")
	log:
		join(LOGDIR, "stringtie/{sample}.log")
	benchmark:
		join(BENCHMARKDIR,"{sample}_stringtie.txt")
	params:
		#strand=STRINGTIE_STRAND,
		#script = join(SCRIPTS, "fix_gtf.py")  # for fixing stringtie gene_id
	#threads: THREADS
	conda:
		join(ENVDIR, "stringtie.yaml")
	shell:
		"""
		nice stringtie {STRINGTIE_STRAND} -p 10 -G {ANNOTATION} {input.string_in} 2> {log}  > {output}
		"""

rule stringtie_merge:
	input:
		expand(join(STRINGTIE, "{sample}.gtf"), sample=SAMPLES)
	output:
		join(STRINGTIE, "stringtie_merged_temp.gtf")
	log:
		join(LOGDIR, "stringtie_merge.log")
	benchmark:
		join(BENCHMARKDIR,"stringtie_merge.txt")
	#threads: THREADS
	params:
		script = join(SCRIPTS, "fix_gtf.py")  # for fixing stringtie gene_id
	conda:
		join(ENVDIR, "stringtie.yaml")
	shell:
		"""
		nice stringtie --merge -p 10 -F {STRINGTIE_F} -G {ANNOTATION} {input} 2> {log}  > {output}
		"""

rule gffcompare:
	## use gff compare to get the class codes see https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
	input:
		join(STRINGTIE, "stringtie_merged_temp.gtf")
	output:
		touch(join(GFF, 'gffcompare.done'))
	log:
		join(LOGDIR, "gffcompare.log")
	benchmark:
		join(BENCHMARKDIR,"gffcompare.txt")
	params:
		#extras = r'-D',
		#script = join(SCRIPTS, "gffcompare"),  # latest version not on bioconda use executable 10.6
		gtf=ANNOTATION,
		outprefix=join(GFF, "gff")
	shell:
		"""
		nice gffcompare -D -r {params.gtf} -o {params.outprefix} {input} 2> {log}
		"""

rule get_tx2gene_map:
	# extract gtf entries for particular classes. 
	# Clean the gtf for gene classes, that is replace the gene_if with gene_name field
	input:
		join(GFF, 'gffcompare.done')
	output:
		join(STRINGTIE, 'tx_gene_map.txt'),
		join(STRINGTIE, "stringtie_merged.gtf")
	log:
		join(LOGDIR, "get_tx2gene_map.stderr")
	benchmark:
		join(BENCHMARKDIR, "get_tx2gene_map.txt")
	params:
		extra  = '-n "i,x,u,y" -g "=,j,k,c,o"',
		script = join(SCRIPTS, "extract_txmap_from_gtf.py"),
		gffanno = join(GFF, "gff.annotated.gtf")
	shell:
		"""
		python {params.script} {params.extra} --gtf {output[1]} {params.gffanno} {ANNOTATION} > {output[0]} 2> {log}
		"""


rule make_salmon_index:
	input:
		join(STRINGTIE, 'stringtie_merged.gtf')
	output:
		join(STRINGTIE, 'transcripts.fa')
	log:
		join(LOGDIR, 'make_salmon_index.txt')
	benchmark:
		join(BENCHMARKDIR,"make_salmon_index.txt")
	conda:
		join(ENVDIR, "gffread.yaml")
	shell:
		"""
		nice gffread -w {output} -g {GENOME} {input} 2> {log} 
		"""

rule salmon_index:
	input:
		join(STRINGTIE, 'transcripts.fa')
	output:
		directory(SALMON_TRANS_INDEX)
	log:
		log1=join(LOGDIR, "salmon_index.err.log"),
		log2=join(LOGDIR, "salmon_index.out.log")
	benchmark:
		join(BENCHMARKDIR,"salmon_idx.txt")
	#threads: THREADS
	conda:
		join(ENVDIR, "salmon.yaml")
	## keep duplicates https://github.com/COMBINE-lab/salmon/issues/214
	shell:
		"""
		nice salmon index -p 10 --keepDuplicates -t {input} -i {output} --type quasi -k 31 > {log.log2} 2> {log.log1}
		"""

rule salmon_quant:
	input:
		r1=join(FASTQ, R1),
		r2=join(FASTQ, R2),
		idx = directory(SALMON_TRANS_INDEX)
	output:
		directory(join(COUNTS2, "{sample}")),
	log:
		log1=join(LOGDIR, "{sample}_salmon_quant.err.log"),
		log2=join(LOGDIR, "{sample}_salmon_quant.out.log")
	benchmark:
		join(BENCHMARKDIR,"{sample}_salmon_pe.txt")
	params:
		extra=SALMON_EXTRA,
	#threads: THREADS
	conda:
		join(ENVDIR, "salmon.yaml")
	shell:
		"""
		nice salmon quant -i {input.idx} -p 10 {params.extra} -1 {input.r1} -2 {input.r2} -o {output[0]} > {log.log2} 2>{log.log1} 
		"""

rule novel1:
	input:
		join(STRINGTIE, 'tx_gene_map.txt'),
		join(STRINGTIE, "stringtie_merged.gtf")
	output:
		ixsuy = temp(join(NOVEL, "ixsuy.novel.gtf")),
		classes = join(NOVEL, 'class_code_list.txt')
	log:
		join(LOGDIR, 'novel1.log')
	benchmark:
		join(BENCHMARKDIR,"novel1.txt")
	params:
		#gffcmp = join(STRINGTIE, "gff.stringtie_merged_temp.gtf.tmap"),
		awk = join(SCRIPTS, "awk-novel1.awk")
	shell:
		## we only want to take class codes 'i', 'u', 'x', 's', 'y' for potential novels
		"""
		tail -n+2 {STRINGTIE}/gff.stringtie_merged_temp.gtf.tmap | awk 'BEGIN{{OFS="\\t";}} {{print$3,$1,$5,$6}}' > {output.classes} 2> {log};
		grep -F -f <(awk -f {params.awk} {input[0]}) {input[1]} > {output.ixsuy} 2>> {log}
		"""

rule novel_gffread:
	# use gffread to turn into fasta for input to CPAT
	input:
		gtf=join(NOVEL, "ixsuy.novel.gtf")
	output:
		temp(join(NOVEL, 'novel.fa'))
	log:
		join(LOGDIR, 'novel_get_fasta.log')
	benchmark:
		join(BENCHMARKDIR, "novel_get_fasta.txt")
	conda:
		join(ENVDIR, "gffread.yaml")
	shell:
		"""
		nice gffread -w {output} -g {GENOME} {input.gtf} 2> {log}
		"""

rule novel_subselect_len:
	input:
		join(NOVEL, 'novel.fa')
	output:
		join(NOVEL, 'novel.gt200.fa')
	log:
		join(LOGDIR, 'subselect_len.log')
	benchmark:
		join(BENCHMARKDIR,"subselect_len.txt")
	conda:
		join(ENVDIR, "bioawk.yaml")
	shell:
		"""
		nice bioawk -c fastx '{{if(length($seq)>200) {{print">"$name;print$seq}}}}' {input} > {output} 2> {log}
		"""

rule novel_cpat:
	input:
		join(NOVEL, 'novel.gt200.fa')
	output:
		join(NOVEL, 'cpat.out')
	log:
		join(LOGDIR, 'cpat.log')
	benchmark:
		join(BENCHMARKDIR,"novel_cpat.txt")
	params:
		cpat1=CPAT_RDATA,
		cpat2=CPAT_TSV
	conda:
		join(ENVDIR, "cpat.yaml")
	shell:
		"""
		nice cpat.py  -g {input} -d {params.cpat1} -x {params.cpat2} -o {output} 2> {log}
		"""

rule novel2:
	input:
		in1=join(NOVEL, 'cpat.out'),
		in2=join(NOVEL,"ixsuy.novel.gtf")
	output:
		out1=join(NOVEL, 'final.novel.cpat.gtf.gz'),
		out2=join(NOVEL, 'novel.cpat.transcript.list')
	log:
		log1=join(LOGDIR, 'novel2.log')
	benchmark:
		join(BENCHMARKDIR,"novel2.txt")
	shell:
		# 1 only take the entries with a coding prediction score of < CPAT_CP
		# 2 now take the ids from the novel.longerThanThreshold.cpat.CP.out file and extract these from
		#   the ixsuy.novel.gtf (the file with only i,x,s,u,y class codes) 
		# 3 extract relevant information to merge the table later in R
		"""
		awk -F "\\t" '{{if ($6<{CPAT_CP}) print $0}}' {input.in1} > {NOVEL}/cpat.CP.out 2> {log};
		awk 'FNR==NR {{f1[$1];next}} $2 in f1' {NOVEL}/cpat.CP.out FS='"' {input.in2} 2>> {log} | gzip > {output.out1};
		zcat {output.out1} | awk 'BEGIN {{OFS="\\t";}} {{if ($3=="transcript" && $13 =="xloc") {{ print $1,$4,$5,$7,$10,$16,"-" }} else if ($3=="transcript" && $13 =="gene_name"){{ print $1,$4,$5,$7,$10,$20,$14 }}}}' 2>> {log} | sed -e s'/"//g' -e 's/;//g' > {output.out2}
		"""

rule novel_gtf_2_bed:
	# convert to bed12-format
	input:
		join(NOVEL, "final.novel.cpat.gtf.gz")
	output:
		join(NOVEL, "final.novel.cpat.transcript.list.bed.gz")
	log:
		join(LOGDIR, 'novel2bed.stderr')
	benchmark:
		join(BENCHMARKDIR,"nonovel2bed.txt")
	params:
		script = join(SCRIPTS, "gtf2bed.pl")
	shell:
		"""
		perl {params.script} {input} | gzip > {output}
		"""

### still leave this in with the -B Ballgown table maker so we can extract intron read coverage info
rule stringtie_counts:
	input:
		string_in=join(MAPPED2, '{sample}.bam'),
		strg=join(STRINGTIE, 'stringtie_merged.gtf')
	output:
		join(COUNTS, '{sample}/{sample}.gtf')
	log:
		join(LOGDIR, '{sample}_stringtie_counts.log')
	benchmark:
		join(BENCHMARKDIR,"{sample}_stringtie_counts.txt")
	#threads: THREADS
	#params: STRINGTIE_STRAND
	conda:
		join(ENVDIR, "stringtie.yaml")
	shell:
		"""
		nice stringtie -B -p 10 -G {input.strg} -o {output} {input.string_in} > {log}
		"""

## ## ############################################################## 
## ## EXPRESSION

## ## a rule to create samplesheet with quant.sf path for tximport
	# we attach the path to the quant.sf file to the sample
rule get_sample_quant_path:
	input:
		sheet=SAMPLESHEET,
		samples=expand(join(COUNTS2, "{sample}"), sample=SAMPLES) # to have it run after salmon has run
	output:
		join(COUNTS2, "sample_list.txt")
	log:
		join(LOGDIR, "get_sample_quant_path.stderr")
	params:
		countdir = COUNTS2,
		script = join(SCRIPTS, "sample_list.py")
	shell:
		"python {params.script} {input.sheet} {params.countdir} > {output} 2> {log}"



## ## a rule to create a matrix of all samples and counts

rule tximport_get_exp:
	input:
		samplelist = join(COUNTS2, "sample_list.txt"),
		t2g = join(STRINGTIE, 'tx_gene_map.txt'),
		design = join(DATA, 'design_matrix_outliers_excluded.csv')
	output:
		join(COUNTS2, "results-tx.txt.gz"),
		join(COUNTS2, "results-gene.txt.gz")
	log:
		join(LOGDIR, "tximport_get_tpm.stderr")
	benchmark:
		join(BENCHMARKDIR,"tximport_get_tpm.txt")
	conda:
		join(ENVDIR, "dexpr.yaml")
	params:
		script = join(SCRIPTS, "build_tx_matrix.R")
	shell:
		"Rscript --vanilla --slave {params.script} {input.samplelist} {input.t2g} {input.design} {output[0]} {output[1]} {THREADS} 2> {log}"

rule get tx_stats:
	input:
		join(COUNTS2, "results-tx.txt.gz"),
		join(GFF, 'gffcompare.done')
	output:
		join(COUNTS2, "results-tx.stats.txt.gz")
	log:
		join(LOGDIR, "get_tx_stats.stderr")
	benchmark:
		join(BENCHMARKDIR,"get_tx_stats.txt")
	conda:
		join(ENVDIR, "dexpr.yaml")
	params:
		script = join(SCRIPTS, "tx_stats.py"),
		gffcmp = join(STRINGTIE, "gff.stringtie_merged_temp.gtf.tmap")
	shell:
		"python {params.script} {input[0]} {params.gffcmp} 2> {log} | gzip > {output}"

rule get_subset_exp_novels:
	input:
		join(NOVEL, "final.novel.cpat.transcript.list.bed.gz"),
		join(COUNTS2, "results-tx.txt.gz")
	output:
		join(COUNTS2, "results-tx.novels.txt.gz")
	log:
		join(LOGDIR, "get_subset_novels.stderr")
	params:
		script = join(SCRIPTS, "subselect.py")
	shell:
		"nice python {params.script} --header -f1 1 -f2 4 {input[1]} {input[0]} 2> {log} | gzip > {output}"

# rule to make some stats for the novel tx
rule get_tx_stats_novels:
	input:
		join(COUNTS2, "results-tx.novels.txt.gz"),
		join(GFF, 'gffcompare.done')
	output:
		join(COUNTS2, "results-tx.novels.stats.txt.gz")
	log:
		join(LOGDIR, "get_tx_stats_novels.stderr")
	benchmark:
		join(BENCHMARKDIR,"get_tx_stats_novels.txt")
	conda:
		join(ENVDIR, "dexpr.yaml")
	params:
		script = join(SCRIPTS, "tx_stats.py"),
		gffcmp = join(STRINGTIE, "gff.stringtie_merged_temp.gtf.tmap")
	shell:
		"python {params.script} {input[0]} {params.gffcmp} 2> {log} | gzip > {output}"


###################
###### circs ######
###################


rule gtf2refFlat:
	## we have to make a geneprediction format annotation file for CIRCexplorer2 from stringtie_merged.gtf
	## this is to catch circs from any novel transcripts picked up in my stringtie_merged.gtf
	input:
		strg=join(STRINGTIE, 'stringtie_merged.gtf')
	output:
		temp(join(CIRCEXP, "tmp.refFlat.txt")),
	log:
		join(LOGDIR, 'gtf2refFlat.log')
	benchmark:
		join(BENCHMARKDIR, "gtf2reFlat.txt")
	conda:
		join(ENVDIR, "gtftogenepred.yaml")
	shell:
		"""
		gtfToGenePred -genePredExt -geneNameAsName2 {input} {output} 2> {log};
		"""

rule combine_anno:
	input:
		join(CIRCEXP, "tmp.refFlat.txt")
	output:
		join(CIRCEXP, 'all.anno.txt')
	log:
		join(LOGDIR, 'combine.anno.stderr')
	params:
		hg38ref=join(DATA, "hg38_ref.txt"),  ##TODO: somewhere there needs to be an explanation where these files come from
		hg38kg=join(DATA, "hg38_kg.txt"),
		awk = join(SCRIPTS, "awk-combine_anno.awk")
	shell:
		"""
		cat {params.hg38ref} {params.hg38kg} <(awk -f {params.awk} {input} | cut -f 1-11) 2> {log} > {output};
		"""


rule circexplorer_parse:
	input:
		join(MAPPED2, '{sample}/mapping2.done')
	output:
		join(CIRCEXP, '{sample}/back_spliced_junction.bed')
	log:
		log1=join(CIRCEXP, '{sample}/CIRCexplorer2_parse.log')
	benchmark:
		join(BENCHMARKDIR, "circexplorer_parse/{sample}.txt")
	conda:
		join(ENVDIR, "circexplorer2.yaml")
	params:
		chim=join(MAPPED2, '{sample}Chimeric.out.junction'),
	shell:
		"""
		CIRCexplorer2 parse -t STAR {params.chim} -b {output} > {log.log1}
		"""

rule circexplorer_annotate:
	# CIRCexplorer2 annotate creates temp files
	# thus only run this rule once at a time as the temp files
	# otherwise get overwritten by concurrent runs
	input:
		bed=join(CIRCEXP, '{sample}/back_spliced_junction.bed'),
		all_anno=join(CIRCEXP, 'all.anno.txt')
	output:
		join(CIRCEXP, '{sample}/circularRNA_known.bed')
	log:
		join(CIRCEXP, 'circexplorer_annotate/{sample}.stdout'),
		join(CIRCEXP, 'circexplorer_annotate/{sample}.stderr')
	benchmark:
		join(BENCHMARKDIR,"circexplorer_annotate/{sample}.txt")
	conda:
		join(ENVDIR, "circexplorer2.yaml")
	# threads: 1000000 really just makes sure that only one is run at a time as I only supply max 10 cores on snakemake.
	# it should be changed to something more general or get the shadow stuff working. Terrible hack....
	threads: 1000000
	shell:
		"""
		CIRCexplorer2 annotate -r {input.all_anno} -g {GENOME} -b {input.bed} -o {output} > {log[0]} 2> {log[1]}
		"""

rule collect_circs:
	input:
		# here I make sure all annotate runs are done before this rule is run...
		expand(join(CIRCEXP, '{sample}/circularRNA_known.bed'), sample=SAMPLES)
	output:
		join(CIRCEXP, 'all_circular_sorted_bed.gz')
	log:
		join(LOGDIR, "collect_circs.log")
	benchmark:
		join(BENCHMARKDIR, "collect_circs.txt")
	shell:
		"cat {input} | sort -uk1,1 -k2,2n -k3,3n 2> {log} | gzip > {output}" 
############
#### QC ####
############

rule rseqc:
	input:
		join(MAPPED2, "{sample}.bam"),
		ANNOTATIONBED
	output:
		touch(join(QC, '{sample}.junction.annot.done')),
		touch(join(QC, '{sample}.junction.sat.done')),
		out2=join(QC, "{sample}.infer.exp.txt"),
		out3=join(QC, "{sample}.read.dist.txt")
	priority: 1
	log:
		log1 = join(LOGDIR, "rseqc_junction_annotation/{sample}.stdout"),
		log2 = join(LOGDIR, "rseqc_junction_annotation/{sample}.stderr"),
		std_out=join(LOGDIR, "{sample}.rseqc.out.stdout"),
		std_err=join(LOGDIR, "{sample}.rseqc.out.stderr")
	benchmark:
		join(BENCHMARKDIR, "rseqc_junction_annotation/{sample}.txt")
	params:
		extra = r'-q 255', # STAR uses 255 as a scrore for uniq mappers
		out = join(QC, '{sample}')
	conda:
		join(ENVDIR, "rseqc.yaml")
	shell:
		"""
		junction_annotation.py {params.extra}  -i {input[0]} -o {params.out} -r {input[1]} > {log.log1} 2> {log.log2};
		junction_saturation.py {params.extra} -i {input[0]} -r {input[1]} -o {params.out} > {log.std_out} 2> {log.std_err};
		infer_experiment.py -r {input[1]} -i {input[0]} > {output.out2} 2> {log.std_err};
		read_distribution.py -i {input[0]} -r {input[1]} > {output.out3} 2> {log.std_err};
		"""

rule mutiqc_samplenames:
	input:
		SAMPLESHEET
	output:
		temp(join(DATA, "multiqc-samplenames.txt"))
	shell:
		"""
		echo -e "samplename" > {output};
		cat {input} | cut -f 1 >> {output}
		"""

rule multiqc:
	input:
		join(DATA, "multiqc-samplenames.txt"),
		expand(join(MAPPED2, '{sample}/mapping2.done'), sample=SAMPLES),
		expand(join(QC, '{sample}.junction.annot.done'), sample=SAMPLES),
		expand(join(QC, '{sample}.junction.sat.done'), sample=SAMPLES),
		expand(join(COUNTS2, "{sample}"), sample=SAMPLES)
	output:
		directory(MULTIQC)
	log:
		log1 = join(LOGDIR, "multiqc.stdout"),
		log2 = join(LOGDIR, "multiqc.stderr")
	conda:
		join(ENVDIR, "multiqc.yaml")
	params:
		configfile = MULTIQC_FILE,
		RSEQC_LOGS = join(LOGDIR, "rseqc_junction_annotation")
	shell:
		"multiqc -f -c {params.configfile} --sample-names {input[0]} -o {output} {MAPPING} {params.RSEQC_LOGS} {QC} {COUNTS2} > {log.log1} 2> {log.log2}"


