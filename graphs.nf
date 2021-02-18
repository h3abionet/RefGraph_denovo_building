#!/usr/bin/env nextflow


params.data_dir = false
params.ref_seq = false
params.out_dir = false
params.assembler = false


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
if ( params.assembler == false) {
    exit 1, "Must specify which assembler to use"
}



/*Prepare input*/

assembly_ch = Channel.fromPath("${params.data_dir}/*.contigs.fasta")
ref_seq = Channel.fromPath(params.ref_seq).tolist()



if(params.assembler == "minigraph") {


/*Run minigraph*/	

process minigraph{
	        publishDir "${params.out_dir}/graphs-out"

	        input:
	        file "*.contigs.fasta" from assembly_ch
			file $ref_seq

	        Output:
	        file "${params.out_dir}/graphs-out/*" into graphs_ch

	        script:
	        """""
			
	        minigraph -xggs -t${task.cpus} $ref_seq "*.contigs.fasta" > out.gfa
			
			gfatools gfa2fa -s out.gfa > out.fa

	        """""
	}


} else if(params.assembler == "novograph"){


/*Run novograph*/

process novograph{	
            publishDir "${params.out_dir}/novographs-out"
			
			input:
	        file "*.contigs.fasta" from assembly_ch
			file $ref_seq

	        Output:
	        file "${params.out_dir}/novographs-out/*" into novographs_ch

	        script:
	        """""
			
            cat *.contigs.fasta > allContigs.fa
			
			perl suggestCommands.pl --inputContigs allContigs.fa --referenceGenome $ref_seq --outputDirectory novographs-out > commands_for_myGraph.txt

	        """""
	}

}

