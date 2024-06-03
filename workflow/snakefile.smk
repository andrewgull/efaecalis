# snakefile for E.faecalis assembly and annotation
from snakemake.io import touch, directory, expand, temp

SRC=config["SRC"]

# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/final/{strain}_all.done", strain=config['strains'])

rule merge_nanopore:
    input: SRC + "{strain}/Nanopore"
    output: temp("resources/{strain}/Nanopore/merged_reads/{strain}_all.fastq.gz")
    threads: 10
    log: "results/logs/{strain}_zcat.log"
    shell: "zcat {input}/*.gz | pigz -c -p {threads} 1> {output} 2> {log}"
# copy BGI using script copy_bgi.sh
# because you need to run QC before the pipeline
# and some params may require adjustment

rule length_filter_nanopore:
    input: "resources/{strain}/Nanopore/merged_reads/{strain}_all.fastq.gz"
    output: temp("results/reads_filtered/{strain}/Nanopore/{strain}_all.fastq.gz")
    message: "executing filtlong on {wildcards.strain} long reads"
    log: "results/logs/{strain}_filtlong.log"
    conda: "filtlong-env"
    threads: 10
    params: min_len=config["min_nanopore_length"]
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule chop_nanopore:
    input: "results/reads_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output: temp("results/reads_filtered/{strain}/Nanopore/{strain}_all_chopped.fastq.gz")
    threads: 10
    log: "results/logs/{strain}_porechop.log"
    conda: "porechop-env"
    shell: "porechop -i {input} -o {output} --threads {threads} &> {log}"

rule trim_reads:
    # overall quality is good
    # per base sequence content in the first 10-15 bases is bad
    input:
        r1 = "resources/{strain}/Illumina/{strain}_1.fq.gz",
        r2 = "resources/{strain}/Illumina/{strain}_2.fq.gz"
    output:
        r1 = temp("results/reads_filtered/{strain}/Illumina/{strain}_1.fq.gz"),
        r2 = temp("results/reads_filtered/{strain}/Illumina/{strain}_2.fq.gz")
    threads: 10
    message: "trimming front ends of the {wildcards.strain} reads"
    log: "results/logs/{strain}_illum_trimming.log"
    conda: "fastp-env"
    params: f = config["trim_front"], adapter1 = config["adapter1"], adapter2 = config["adapter2"]
    shell: "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --thread {threads} "
           "--trim_front1 {params.f} --trim_front2 {params.f} --adapter_sequence {params.adapter1} "
           "--adapter_sequence_r2 {params.adapter2} &> {log}"

rule adaptive_hybrid_assembly:
    input:
        short_reads_1 = "results/reads_filtered/{strain}/Illumina/{strain}_1.fq.gz",
        short_reads_2 = "results/reads_filtered/{strain}/Illumina/{strain}_2.fq.gz",
        long_reads = "results/reads_filtered/{strain}/Nanopore/{strain}_all_chopped.fastq.gz"
    output:
        assembly_dir = directory("results/assemblies/{strain}"),
        draft_dir = directory("results/drafts/{strain}"),
        polish_dir = directory("results/polished/{strain}")
    threads: 18
    message: "executing assembly script with {threads} threads on {wildcards.strain} reads"
    log: "results/logs/{strain}_assembly.log"
    conda: "uni+fmp-env"
    params: basecaller=config["basecaller"], genome_size=config["genome_size"], coverage=config["coverage"], genome_length=config["genome_length"], cov_threshold=config["cov_threshold"]
    script: "scripts/adaptive_hybrid_assembly.py"
    
rule assembly_summary:
    input: "results/assemblies/{strain}"
    output: "results/summaries/{strain}_assembly_summary.tsv"
    threads: 1
    message: "summarizing unicycler and SPAdes assemblies of strain {wildcards.strain}"
    log: "results/logs/{strain}_assembly_summary.log"
    script: "scripts/assembly_summary.py"

rule qc_assembly:
    input: "results/assemblies/{strain}", "resources/busco_downloads"
    output: directory("results/qualcheck_assembly/{strain}")
    threads: 18
    message: "executing BUSCO and QUAST with {threads} threads on {wildcards.strain} assembly"
    conda: "busco-quast-env"
    params: tax_dataset=config["tax_dataset"]
    script: "scripts/QC_assembly.py"

rule annotation:
    # proteins (prodigal) and rRNA (barrnap)
    input: "results/assemblies/{strain}"
    output: directory("results/annotations/{strain}/prokka")
    threads: 10
    message: "executing PROKKA with {threads} threads on full assembly of {wildcards.strain}"
    log: "results/logs/{strain}_prokka.log"
    conda: "prokka-env"
    params: centre=config["centre"], minlen=config["minlen"], genus=config["genus"], species=config["species"]
    shell:
        # skip tRNAs search?
        "prokka --addgenes --addmrna --compliant --notrna --outdir {output} --prefix {wildcards.strain}_genomic --centre {params.centre} --genus {params.genus} "
        "--species {params.species} --strain {wildcards.strain} --kingdom Bacteria --cpus {threads} "
        "--mincontiglen {params.minlen} {input}/assembly.fasta &> {log}"

rule final:
    input: busco="results/qualcheck_assembly/{strain}",
           prk="results/annotations/{strain}/prokka",
           sum="results/summaries/{strain}_assembly_summary.tsv"

    output: touch("results/final/{strain}_all.done")
    shell: "echo 'DONE'"