#!/usr/bin/env nextflow

params.read1
params.read2
params.read_1_primers
params.read_2_primers
params.output
params.nproc = 2
params.pipeline_log = 'brc_assemble.log'

process read_assemble {
    output:
    file "full_seq_file.fasta" into full_fasta	

    """
    python3 brc_assemble.py --read1 $params.read1 --read2 $params.read2 --output full_seq_file.fasta
    """
}

process mask_forward_primer {
    input:
    file f_seq_file from full_fasta
    
    output:
    file "forward_cut_seq_primers-pass.fasta" into forward_fasta    

    """
    MaskPrimers.py align -s $f_seq_file -p $params.read_2_primers --mode cut --pf VPRIMER --outname forward_cut_seq --fasta --log MPV.log --nproc $params.nproc >> $params.pipeline_log
    """
}

process mask_reverse_primer {
   input:
   file r_seq_file from forward_fasta

   """
   MaskPrimers.py align -s $r_seq_file -p $params.read_1_primers \
	--mode cut  --revpr --pf CPRIMER \
	--outname $params.output --log MPC.log --nproc $params.nproc >> $params.pipeline_log
   """

}
