


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
