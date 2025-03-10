# command to run the program: snakemake -c8

# Genome source directory 
source_directory = "/staging/leuven/stg_00079/teaching/1000genomes/"
# Reference genome database
genome_db = "/git@github.com:Kodai0531/Snakemake-project.gitustre1/project/stg_00079/teaching/hg38_21/chr21.fa"
# Path to the snpEff database folder.
snpeff_jar = ("/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar")
# Extracts sample names from fastq files.
snpeff_db_folder = "/staging/leuven/stg_00079/teaching/snpeff_db/"

# Data file name pattern: HG0???8.*.fq.gz 
# Capture the target file names ends with '8' in the directory.
patterns, = glob_wildcards(source_directory + "{pattern}8.GRCh38DH.exome.chr21.fq.gz")
sample_names = [pattern + '8.GRCh38DH.exome.chr21' for pattern in patterns]

print("-"*50)
print("Samples:\n", sample_names)
print("The number of sample files:", len(sample_names))
print("-"*50)

rule all:
    input:
        expand("010.fastqc/{sample}_fastqc.zip", sample=sample_names), # rule fastqc
        "genes.vcf" # rule extract_genes

# 00.Copy & Unzip sample files
rule copy:
    input:
        genome = expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz",sample=sample_names)
    output:
        data = expand("000.fastq/{sample}.fq",sample=sample_names)
    params:
        source = source_directory
    log:
        err="copy.err"
    shell:
        """
        # Ensure the directory exists
        mkdir -p 000.fastq        
        # Copy genome data
        cp -v {input.genome} 000.fastq/

        gunzip -v 000.fastq/*.fq.gz

        # List the first 10 files in the directory
        ls -lh 000.fastq/ | head -n 10

        # Checks if file successfully copied 
		if [ ! -s {output.data} ]; then
			echo "Error in rule copy: genome copy {output.data} don't exist" >> {log.err}
		    exit 1
		fi
        """

# 01.Quality Check
rule fastqc: 
    input: 
        fq = "000.fastq/{sample}.fq"
    output: 
        zip = "010.fastqc/{sample}_fastqc.zip",
        html = "010.fastqc/{sample}_fastqc.html",
        summarytxt = "010.fastqc/{sample}_fastqc/summary.txt",
        basequal = "010.fastqc/{sample}_fastqc/Images/per_base_quality.png"
    shell:
        """
        echo "Running FastQC..."
        fastqc -o 010.fastqc {input.fq} --extract 
        
        echo "## Quality Control for fastqc ##"  
        # Check if the line count is divisible by 4; proper fastq file formatting
        if [ $(( $(wc -l < {input.fq}) % 4)) -eq 0 ]; then
            echo "The number of reads in {input.fq} file is divisible by four." 
        else
            echo "The number of reads in {input.fq} file is not divisible by four."
            exit 1
        fi
        """
 
# 02.BWA
# Align the reads to the reference genome(genome_db) using BWA.
rule bwa_map:
    input: 
        fq = "000.fastq/{sample}.fq"
    output: 
        bam = "020.bwa/{sample}.bam",
        bai = "020.bwa/{sample}.bam.bai" 
    params:
        db = genome_db
    shell:
        """
        bwa mem {params.db} {input.fq} | samtools sort -o {output.bam}
        samtools index {output.bam}

        # Check number of 'primary mapped' reads
        echo "## Quality control for BWA,samtools ##"
        echo "How many 'primary mapped' reads?:"
        samtools flagstat {output.bam} | grep 'primary mapped'
        """

# 03.SNP Calling: bcftools
rule snp_calling:
    input: 
       bam = expand("020.bwa/{sample}.bam",sample=sample_names),
       bai = expand("020.bwa/{sample}.bam.bai", sample=sample_names)
    output:
        vcf = "030.vcf/snp_raw.vcf",
        stats = "030.vcf/snp_raw.stats" # stats report
    params:
        db = genome_db    
    log:
        err="030.vcf/error_log.err"
    shell:
        """
        # Generates genotype likelihoods from BAM & Call variants (determine the most likely genotypes and variant sites.)
        bcftools mpileup -f {params.db} {input.bam} \
        | bcftools call -mv -Ov - > {output.vcf}

        echo "## Quality control for snp_calling##"
        # Produce VCF stats
        bcftools stats {output.vcf} > {output.stats}

        # Produce visual plots
        plot-vcfstats -P -p 030.vcf/Visual_plots/ {output.stats}
       
		# Checks if file exist 
		if [ ! -s {output.vcf} ]; then
			echo "Error in snp_calling: VCF file {output.vcf} don't exist" >> {log.err}
		    exit 1
		fi
        """
 
# 04.Normalizing & Filtering
rule cleanup:
    input:
        vcf="030.vcf/snp_raw.vcf"
    output:
        vcf_cleaned="040.cleaned/snp_cleaned.vcf"
    params:
        db=genome_db
    log:
        cleanup_log="040.cleaned/cleanup.log" # standard log(.log)
    shell:
        """
        set -euo pipefail
        echo "processing cleanup..." >> {log.cleanup_log}

        # Decompose & Normalize & Filter 
        cat {input.vcf} \
        | vt decompose - \
        | vt normalize -n -r {params.db} - \
        | vt uniq - \
        | vt view -f "QUAL>20" -h - \
        > {output.vcf_cleaned}
        
        # Check if output file was created
        echo "## Quality control for cleanup ##"  >> {log.cleanup_log}
        if [ ! -s {output.vcf_cleaned} ]; then
            echo "Error: The output VCF file is empty or does not exist." >> {log.cleanup_log}
            exit 1
        fi
        
        # Verify file structure: Ensure file has content and proper header
        if ! grep -q '^#' {output.vcf_cleaned}; then
            echo "Error: VCF header missing, file structure is incorrect." >> {log.cleanup_log}
            exit 1
        fi
        """

# 05.Annotation
rule snpeff: 
    input: 
        vcf = "040.cleaned/snp_cleaned.vcf"
    output: 
        vcf_annotated = "050.snpeff/snp_annotated.vcf",
        html = "050.snpeff/snpEff_summary.html",
        genetxt = "050.snpeff/snpEff_genes.txt"
    params:
        snpeff_db_folder = snpeff_db_folder,
        snpeff_jar = snpeff_jar
    log: 
        err="050.snpeff/snp_annotated.err" #Error Output Log (.err)
    shell:
        """
        # Annoatate SNP using SnpEff and capture standard output to a file
        java -Xmx3400m -jar \
            {params.snpeff_jar} eff hg38 \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.vcf_annotated}
        
        # Move output files to the porper folder so not to cause a problem
        mv snpEff_genes.txt snpEff_summary.html 050.snpeff

        # Check if output file was created
        echo "## Quality control for annotation ##"  >> {log.err}
        if [ ! -s {output.vcf_annotated} ]; then
            echo "Error: The output annotated VCF file is empty or does not exist." >> {log.err}
            exit 1
        fi
        """

# Rule for extracting SNPs of genes of interest: "APP", "SOD1" and "DYRK1A"
rule extract_genes:
    input:
        vcf = "050.snpeff/snp_annotated.vcf"
    output:
        genes_vcf = "genes.vcf"
    shell:
        """
        # Extract lines with gene names "APP", "SOD1" and "DYRK1A" from annotated vcf file
        grep '#' {input} > {output.genes_vcf}
        grep -E '\\|APP\\||\\|SOD1\\||\\|DYRK1A\\|'\
        {input.vcf} > {output.genes_vcf}
        
        # Checks if file exist and has a size greater than 0
		if [ ! -s {output.genes_vcf} ]; then 
			echo 'Error in extract_genes rule: VCF file {output.genes_vcf} is not created or is empty' >> error_log.txt
			exit 1 
		fi
        """

# rule heatmap:
