#!/usr/bin/env nextflow

params.data_dir = false
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

/*Prepare input*/

input_ch = Channel.fromPath("${params.data_dir}/*.bam")
ref_seq = Channel.fromPath(params.ref_seq).tolist()
gen_ref = Channel.fromPath(params.gen_ref).tolist()


//Validate inputs

if ( params.data_dir == false) {
    exit 1, "Must specify path to directory containing bam files"
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

/* Run the canu genome assembler. A sample prefix, output directory and reference sequence 
   must be specified in the config file in order for canu to run. It will use default 
   parameters beyond that but additional options can be set in the .config file*/
 
process runCanu {
    input:
    file "${input.baseName}.fq" from fastq_ch
    file orig from orig_ch
    
	output:
    file "${orig/Assembly/*}" into assembly_ch
  
    script:	
    """
    canu -p ${params.sample_prefix} \
	 -d "$orig/Assembly/" \
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
    input:
    file "${orig/Assembly/}*.contigs.fasta" from assembly_ch 
    file $ref_seq
    file $gen_ref
	
    Output:
    file "${orig/quast/}" into quast_ch
    
    Script:
    """""
    quast.py "${orig/Assembly/}*.contigs.fasta" \
        -r ${ref_seq} \
        -g ${gen_ref} \
        -o $orig/quast/ \
	${params.genome_large} ${params.genome_type} ${params.quast_options}
    """""
}

/* To add:
  * Need to add a loop script to allow this to run across multiple folders/files. 
  * By default the script will run for any bam files contained within the input directory.
  * It then needs to feed into minigraph (requires at least 2 genomes + a reference to run).
  * We may also need to add an additional process to phase in haplotypes - currently
    the script builds a genome with many repeats resulting in an assembly much larger 
    than the reference*/
