import pysam
import sys
import re
from Bio.Seq import Seq
import subprocess
import os
import glob

def MateDiscord(rec, mRNA_len):
    if rec.next_reference_name != rec.reference_name:
        return(True)
    elif rec.next_reference_name == rec.reference_name and (rec.next_reference_start-rec.reference_start >= mRNA_len or rec.reference_start-rec.next_reference_start > 0): # if on the same chromosome and distance > transcript length, count as discordant pair
        return(True)
    else:
        return(False)

def SaveClipped(in_bam, out_fasta, utr_combined, min_overlap, min_clip, qname_anno_table, gene_info):

    # including primary assembly chromosomes only
    chrs = set(["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"])

    input = pysam.AlignmentFile(in_bam, 'rb')
    output = open(out_fasta, mode='w')
    mRNA_len, side, chrom, start, end = utr_combined[1], utr_combined[2], utr_combined[3], utr_combined[4], utr_combined[5]
    d = int(min_overlap)
    n = int(min_clip)

    qname_meta = open(qname_anno_table, mode = "a")
    info =re.findall(r"^(.*)_(.*)_(.*)_.*fasta$" , os.path.basename(out_fasta)) 
    sample_name = info[0][0]
    gene_name = info[0][1]
    utr_name = info[0][2]

    utr_len = end - start + 1 # in case the length of the utr is too short, i.e. utr_len < 2*d, we do not require minimal covered bases
    if utr_len > 2 * d: 
        start = start + d
        end = end - d
    
    for read in input.fetch(chrom, start, end):
        #if not read.is_unmapped and not read.mate_is_unmapped:
        if not read.is_unmapped and read.reference_name in chrs: # mate not necessarily unmapped
            if len(read.cigarstring) == 0: continue

            utr_read_strand = "+"
            soft_l = re.findall(r"^(\d+)[S]", read.cigarstring)
            soft_r = re.findall(r"(\d+)[S]$", read.cigarstring)
            if side == "right":
                if len(soft_r) != 0 and int(soft_r[0]) >= n and MateDiscord(read, mRNA_len):
                    pos = len(read.query_sequence)-int(soft_r[0])
                    fastr_r = read.query_sequence[pos:] # output right side soft-clipped sequences

                    utr_pos = "{}:{}-{}".format(chrom, str(read.reference_start+1), str(read.reference_end))
                    utr_read_seq = read.query_sequence[:pos]

                    #if mapped to minus strand, reverse complement the sequence
                    if read.is_reverse: 
                       fastr_r = str(Seq(fastr_r).reverse_complement())
                       utr_read_strand = "-"

                    output.write(">%s\n%s\n" % (read.query_name, fastr_r))
                    qname_meta.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.query_name, sample_name, gene_name, gene_info, utr_name, side, utr_read_strand, utr_pos, utr_read_seq, fastr_r))

        ########################################################################## only consider left side for 3UTR
            elif side == "left":
                if len(soft_l) != 0 and int(soft_l[0]) >= n and MateDiscord(read, mRNA_len):
                    pos = int(soft_l[0])
                    fastr_l = read.query_sequence[:pos]

                    utr_pos = "{}:{}-{}".format(chrom, str(read.reference_start+1), str(read.reference_end))
                    utr_read_seq = read.query_sequence[pos:]

                    #if mapped to minus strand, reverse complement the sequence
                    if read.is_reverse: 
                       fastr_l = str(Seq(fastr_l).reverse_complement())
                       utr_read_strand = "-"
                    
                    output.write(">%s\n%s\n" % (read.query_name, fastr_l))
                    qname_meta.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.query_name, sample_name, gene_name, gene_info, utr_name, side, utr_read_strand, utr_pos, utr_read_seq, fastr_l))


    input.close()
    output.close()
    qname_meta.close()
    if os.stat(out_fasta).st_size == 0: os.remove(out_fasta) # remove fasta output if it is empty