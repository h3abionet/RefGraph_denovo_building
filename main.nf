#!/usr/bin/env nextflow

params.data_dir = false
params.out_dir = false
params.ref_seq = false
params.gen_ref = false
params.genomeSize = false
params.time = false
params.partition = false
params.sample_prefix = false
params.reads_corr = false
params.file_type = false
params.canu_options = false
params.genome_large = true /*Default*/
params.genome_type = true /*Default*/
params.quast_options = false
params.assembler = false
params.no_cpus = false

/*Prepare input*/

input_ch = Channel.fromPath("${params.data_dir}/*.bam")
ref_seq = Channel.fromPath(params.ref_seq).tolist()
gen_ref = Channel.fromPath(params.gen_ref).tolist()


//Validate inputs

if ( params.data_dir == false) {
    exit 1, "Must specify path to directory containing bam files"
}
if ( params.out_dir == false) {
    exit 1, "Must specify path to output directory"
}
if ( params.ref_seq == false) {
    exit 1, "Must specify path to the reference sequence"
}

if ( params.genomeSize == false) {
    exit 1, "Please specify the size of your genome as <integer>[g,m,k] e.g. 2.8g"
}
if ( params.sample_prefix == false) {
    exit 1, "Must specify sample prefix"
}
if ( params.file_type == false) {
    exit 1, "Must specify file/data type e.g. -pacbio-hifi"
}

	
		
/* Extract fastq files from bam files*/


process ExtractFastq {
    publishDir "${params.out_dir}/fastq-raw"
	
    input:
    file input from input_ch

    output:
    file "${input.baseName}.fq" into fastq_ch
    file "$input" into orig_ch
    
    script:
    """
    bedtools bamtofastq -i $input -fq "${input.baseName}.fq"
    """
}

if (params.assembler == "canu") {

	/* Run the canu genome assembler. A sample prefix, output directory and reference sequence 
	   must be specified in the config file in order for canu to run. It will use default 
	   parameters beyond that but additional options can be set in the .config file*/
	 
	process runCanu {
		publishDir "${params.out_dir}/canu-out"
		
		input:
		file "${input.baseName}.fq" from fastq_ch
		
		output:
		file "*.contigs.fasta" into assembly_ch
	 
		script:	
		"""
		canu -p ${params.sample_prefix} \
		     -d "${params.out_dir}/canu-out" \
		     genomeSize=${params.genome_size} \
		     gridOptions="--time=${params.time} --partition=${params.partition} \
		     ${params.file_type} ${params.reads_corr} "${input.baseName}.fq"\
		     ${params.canu_options}	
		"""
	}
	/* Run quast on the resulting assembly. Additional parameters may be added using the 
	   config file. --large and --eukaryote is enabled by default. Please change this if 
	   you do not want this enabled.*/

	process quast {
	        publishDir "${params.out_dir}/quast-canu-out"

	        input:
	        file "*.contigs.fasta" from assembly_ch 
	        file $ref_seq
	        file $gen_ref

	        Output:
	        file "${params.out_dir}/quast-canu-out/*" into quastc_ch

	        script:
	        """""
	        quast.py "*.contigs.fasta" \
		       -r ${ref_seq} \
		       -g ${gen_ref} \
		       -o ${params.out_dir}/quast-canu-out/ \
		       ${params.genome_large} ${params.genome_type} ${params.quast_options}
	        """""
	}
}

else if (params.assembler == "hifiasm") {

	/* will concatenate and zip fq files into .fq.gz format required for hifiasm*/

	process comp_fq {
		publishDir "${params.out_dir}/fq_zip"
		
		input
		file "${input.baseName}.fq" from fastq_ch
		
		output
		file "*.fq.gz" into zip_ch
		
		script
		"""
		#!/usr/bin/env bash
		
		for { name in ./*.fq; do {
			if [[ "$name" = ([0-9][a-z]+)_.*(..)\.fq ]]; then
			outfile="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}.fq"
			cat "$name" >> "$outfile"
			fi
			}
			then {
			
			gzip "$outfile"
			}
		}
		"""
		
	}

	/* will run the hifiasm assembler*/
	/* run hifiasm for pacbio-hifi reads only!*/
	
	process runhifi {
		publishDir "${params.out_dir}/hifiasm"
		
		input:
		file "*.fq.gz" from zip_ch
		
		output:
		file *p_ctg.gfa into contig_ch
		
		script:
		"""
		hifiasm -o "${params.sample_prefix}_hifi.asm" -t ${params.no_cpus} "*.fq.gz" /*need to add t to 
to the list of parameters: specified number of cpus*/
		"""
	}
	
	process quast {
	        publishDir "${params.out_dir}/quast-hifi-out"

	        input:
	        file *p_ctg.gfa from contig_ch
	        file $ref_seq
	        file $gen_ref

	        Output:
	        file "${params.out_dir}/quast-hifi-out/*" into quasth_ch

	        script:
	        """""
	        quast.py "*p_ctg.gfa" \
		       -r ${ref_seq} \
		       -g ${gen_ref} \
		       -o "${params.out_dir}/quast-hifi-out/" \
		       ${params.genome_large} ${params.genome_type} ${params.quast_options}
	        """""
	}
}




/* To add:
  * Need to add a loop script to allow this to run across multiple folders/files. 
  * By default the script will run for any bam files contained within the input directory.
  * It then needs to feed into minigraph (requires at least 2 genomes + a reference to run).
  * We may also need to add an additional process to phase in haplotypes - currently
    the script builds a genome with many repeats resulting in an assembly much larger 
    than the reference*/
