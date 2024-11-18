#!/usr/bin/env python
# Main script to identify UTR-spanning insertion sites; it iterates through a list of target genes and checks one sample at a time; multiple samples will be combined and parallelized by GNU parallel

# Created by Yunjiao Zhu 3-7-2020

import gffutils
import statistics
import argparse
import sys
import subprocess as sp
#from tqdm import *
import os
import re
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
import functools
import gffutils.pybedtools_integration
import glob
import pysam
import multiprocessing
from joblib import Parallel, delayed
from Extract_clipped_bases_FUN import MateDiscord, SaveClipped
from Filter_mate_read_FUN import check_clipping_and_edit_disc, FilterMate

os.chdir(sys.path[0])

def reduce_concat(x, sep=""):
    return functools.reduce(lambda x, y: str(x) + sep + str(y), x)

def paste(*lists, sep=" ", collapse=None):
    result = map(lambda x: reduce_concat(x, sep=sep), zip(*lists))
    if collapse is not None:
        return reduce_concat(result, sep=collapse)
    return list(result)

def Extract_clipping(gene_name, query_chrom, query_start, query_end, utr_out, in_bam, path_sample_name, qname_anno_table, min_overlap, min_clip):
    # call Extract_clipped_bases.py
    count = 1
    for i in utr_out:
        out_fasta = paste([path_sample_name], [gene_name], [i[0]], [count], sep="_", collapse=None)[0] + ".fasta"
        gene_info = "{}:{}-{}".format(query_chrom, query_start, query_end)
        SaveClipped(in_bam, out_fasta, i, min_overlap, min_clip, qname_anno_table+".dup", gene_info)
        count += 1

def Get_discord_reads(in_bam, path_sample_name, outSam_all_discord ,out_folder, sample_name):
    fasta_wildcard = paste([path_sample_name], ["*"], sep="_", collapse=None)[0] + ".fasta"

    # store file status: if empty, store 0
    fastas = glob.glob(fasta_wildcard)
    test_fasta = []
    for file in fastas:
        if os.path.getsize(file) == 0:
            test_fasta.append(0)
        else:
            test_fasta.append(1)

    if not any(test_fasta):    
        print("{} completed: no UTR clipped reads found!\n".format(sample_name))
        with open(os.path.join(out_folder,"progress.track"),'a') as p:
            p.write("{} completed: no UTR clipped reads found!\n".format(sample_name))
        sys.exit()
        
    else:
        cmd_ids = paste(["cat"], [fasta_wildcard], ["| grep -E '^>.*$' | awk '{print $1}' | tr -d '>'"], sep=" ",collapse=None)[0] # skip everything after the first space
        process = sp.Popen(cmd_ids, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
        ids = process.communicate()[0]
        ids = set(ids.rstrip('\n').split('\n'))

        # filter bam file by qnames to get all discordant pairs for the clipped reads
        input = pysam.AlignmentFile(in_bam, 'rb')
        output = pysam.AlignmentFile(outSam_all_discord, 'w', template=input)
        for rec in input.fetch(until_eof=True):
            if rec.query_name in ids: output.write(rec)
        input.close()
        output.close()


def Filter_discord_mates(outSam_all_discord, path_sample_name, gene_name, gene_loc):
        parse = re.findall(r"^(.*):(\d+)-(\d+)", gene_loc)
        chrom = parse[0][0]
        start = int(parse[0][1])
        end = int(parse[0][2])

        outSam_mate_discord = paste([path_sample_name], [gene_name], ["filtered_mate_discord"], sep="_", collapse=None)[0] + ".sam"
        outSam_mate_discord_lowQUAL = paste([path_sample_name], [gene_name], ["lowQUAL_mate_discord"], sep="_", collapse=None)[0] + ".sam"
        FilterMate(outSam_all_discord, outSam_mate_discord, outSam_mate_discord_lowQUAL, chrom, start, end, cutoff=0.05)

def chunks(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def main():
    parser = argparse.ArgumentParser(description='Identify UTR-insertion site gencDNA reads')
    parser.add_argument('Input_Bam_file')
    parser.add_argument('-o', '--output_folder', required=True)
    parser.add_argument('-f', '--gene_list', required=True, help='A text file containing gene names; each gene occupies a line.')
    parser.add_argument('-g', '--utr_gap', type=int, default=15)
    parser.add_argument('-m', '--min_overlap', type=int, default=15, help='Num. of bases a read much overlap the UTR to be considered.')
    parser.add_argument('-c', '--min_clip', type=int, default=20)
    parser.add_argument('-t', '--num_threads', type=int, default=8)
    parser.add_argument('-s', '--chunk_size', type=int, default=16)
    parser.add_argument('--database', default="/tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/GRCh38_GCF_000001405.40-RS_2023_10_Prim_gene_mRNA_UTR_chrName_UCSC.db")

    args = vars(parser.parse_args())

    db = gffutils.FeatureDB(args['database'], keep_order=True)
    in_bam = os.path.abspath(args['Input_Bam_file'])
    out_folder = os.path.abspath(args['output_folder'])
    #if not os.path.exists(out_folder): os.makedirs(out_folder)

    sample_name = os.path.basename(in_bam)
    path_sample_name = os.path.join(out_folder,sample_name)

    qname_anno_table = path_sample_name+"_qname_anno.table"
    outSam_all_discord = paste([path_sample_name], ["UTR_unfiltered_discord"], sep="_", collapse=None)[0] + ".sam"

    blast_DB = "/tscc/nfs/home/yzhu/blast_ncbi_DBs_V5/human_genome/GCF_000001405.39_top_level"
    blast_header = 'read_id\tref_id\tquery_length\talignment_length\tprct_exact_match\tquery_start\tquery_end\tref_start\tref_end\tEvalue\tBitscore\tblast_strand'
    
    print("Begin processing "+sample_name+".")

    genelist = [line.rstrip('\n') for line in open(args['gene_list'])]
    num_genes = len(genelist)

    gene_chunks = list(chunks(genelist, args['chunk_size']))
    for chunk in gene_chunks:
        chunk_utr_list = []
        chunk_gene_list = []
        for gene_name in chunk:
            utr_out = []
            ## if the output already generated, skip the current iteration; use to resume unfinished job
            #if os.path.isfile(paste([path_sample_name], [gene_name], ["filtered_mate_discord"], sep="_", collapse=None)[0]): continue
            query = 'gene-'+gene_name

            # check if gene has the "mRNA" feature
            ##############################
            gene_features = set([i.featuretype for i in db.children(query)])
            if "mRNA" not in gene_features:
                #print(gene_name+" has no annotated mRNA. Skip!")
                continue
            else:
                if ("three_prime_UTR" not in gene_features) and ("five_prime_UTR" not in gene_features):
                    #print(gene_name+" has no annotated 3' or 5' UTRs. Skip!")
                    continue
            ##############################

            # get median transcript length
            mRNA_len = []
            for i in db.children(query, featuretype='mRNA', order_by='start'):
                mRNA_len.append(len(i))
            median_mRNA_len = int(statistics.median(mRNA_len))

            # determine which side of the UTR to look for insertions
            gene = db[query]
            query_chrom = gene.seqid
            query_start = gene.start
            query_end = gene.end

            if gene.strand == "+":
                take_clip_5utr = "left"
                take_clip_3utr = "right"
            elif gene.strand == "-":
                take_clip_5utr = "right"
                take_clip_3utr = "left"
                
            chunk_gene_list.append((gene_name, query_chrom, query_start, query_end))

            # collapse multiple UTRs
            five_utrs = db.children(query, featuretype='five_prime_UTR', order_by='start')
            bed_five_utrs = gffutils.pybedtools_integration.to_bedtool(five_utrs)
            five_utr_reduced = bed_five_utrs.merge(d=args['utr_gap'], s=True) # allowing 15bp gap; merge in a strad-specific way
            for i in five_utr_reduced:
                #utr_out.append( ("5utr", "%s %s %s:%s-%s" % (median_mRNA_len, take_clip_5utr, i[0], str(int(i[1])+1), i[2])) ) # convert 0-based coordinate to 1-based
                utr_out.append(("5utr", median_mRNA_len, take_clip_5utr, i[0], int(i[1])+1, int(i[2])))

            three_utrs = db.children(query, featuretype='three_prime_UTR', order_by='start')
            bed_three_utrs = gffutils.pybedtools_integration.to_bedtool(three_utrs)
            three_utr_reduced = bed_three_utrs.merge(d=args['utr_gap'], s=True) # allowing 15bp gap; merge in a strad-specific way
            for i in three_utr_reduced:
                #utr_out.append( ("3utr", "%s %s %s:%s-%s" % (median_mRNA_len, take_clip_3utr, i[0], str(int(i[1])+1), i[2])) )
                utr_out.append(("3utr", median_mRNA_len, take_clip_3utr, i[0], int(i[1])+1, int(i[2])))
            chunk_utr_list.append(utr_out)

        Parallel(n_jobs=args['num_threads'], verbose=2)(delayed(Extract_clipping)(chunk_gene_list[i][0], chunk_gene_list[i][1], chunk_gene_list[i][2], chunk_gene_list[i][3], chunk_utr_list[i], in_bam, path_sample_name, qname_anno_table, args['min_overlap'], args['min_clip']) for i in range(len(chunk_utr_list)))

    ### dedup qname_anno_table
    with open(qname_anno_table+".dup", "r") as input:
        lines = set(input.readlines())
    with open(qname_anno_table, "w") as output:
        for line in lines:
            output.write(line)
    os.remove(qname_anno_table+".dup")

    # filter bam file to get all raw discordant reads
    Get_discord_reads(in_bam, path_sample_name, outSam_all_discord ,out_folder, sample_name)

    # extract and quality filter mates
    gene_info = pd.read_csv(qname_anno_table, sep='\t', header=None)
    genes = gene_info.iloc[:,2].unique().tolist()
    locs = gene_info.iloc[:,3].unique().tolist()

    Parallel(n_jobs=args['num_threads'], verbose=5)(delayed(Filter_discord_mates)(outSam_all_discord, path_sample_name, genes[i], locs[i]) for i in range(len(genes)))

    # consolidate the BAM files generated for each gene into one BAM file/sample

    ### for mapped mate reads ###
    filtered_mate_sams = glob.glob(path_sample_name+"_*_filtered_mate_discord.sam")
    lowQUAL_mate_sams = glob.glob(path_sample_name+"_*_lowQUAL_mate_discord.sam")

    if len(filtered_mate_sams) or len(lowQUAL_mate_sams) != 0:
        qnames = []
        highQ_mate = 0
        lowQ_mate = 0

        if len(filtered_mate_sams) != 0:
            highQ_mate = 1
            print("Comparing clipped sequences with high-quality mate reads in "+sample_name+".")

            ##### export headless sam file (fields 1-10 + strand) and deduped qnames ############
            merge_para = ["-f", path_sample_name+"_merged_UTR_discord_mates.bam"] + filtered_mate_sams
            pysam.merge(*merge_para)
            input = pysam.AlignmentFile(path_sample_name+"_merged_UTR_discord_mates.bam", 'rb')
            out_headless_sam = open(path_sample_name+"_merged_UTR_discord_mates_noheader.temp.sam", mode='w')
            for rec in input.fetch(until_eof=True):
                qnames.append(rec.query_name)
                strand = "-" if rec.is_reverse else "+"
                out_headless_sam.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rec.query_name, rec.flag, rec.reference_name, rec.reference_start, rec.mapping_quality, rec.cigarstring, rec.next_reference_name, rec.next_reference_start, rec.template_length, rec.query_sequence, strand))
            
            input.close()
            out_headless_sam.close()

    ### for unmapped or poorly mapped mate reads ###

        if len(lowQUAL_mate_sams) != 0:
            lowQ_mate = 1
            print("Comparing clipped sequences with poorly mapped mate reads (for polyA search) in "+sample_name+".")

            ##### export headless sam file (fields 1-10 + strand) and deduped qnames ############
            merge_para_low = ["-f", path_sample_name+"_merged_UTR_lowQUAL_mates.bam"] + lowQUAL_mate_sams
            pysam.merge(*merge_para_low)
            input = pysam.AlignmentFile(path_sample_name+"_merged_UTR_lowQUAL_mates.bam", 'rb')
            out_headless_sam = open(path_sample_name+"_merged_UTR_lowQUAL_mates_noheader.temp.sam", mode='w')
            for rec in input.fetch(until_eof=True):
                qnames.append(rec.query_name)
                if rec.is_unmapped: strand = "NA"
                if not rec.is_unmapped: strand = "-" if rec.is_reverse else "+"
                out_headless_sam.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rec.query_name, rec.flag, rec.reference_name, rec.reference_start, rec.mapping_quality, rec.cigarstring, rec.next_reference_name, rec.next_reference_start, rec.template_length, rec.query_sequence, strand))
            input.close()
            out_headless_sam.close()

        qnames = set(qnames)
        with open(path_sample_name+"_mate.temp.qnames", mode='w') as out_qnames:
            for i in qnames:
                out_qnames.write("%s\n" % i)
        ##########################

        sed_para = '"'+'1i'+blast_header+'"'
        cmd_call_blast = paste(["cat"], [path_sample_name+"_*.fasta"], ["| seqtk subseq -"], [path_sample_name+"_mate.temp.qnames"], ["| blastn -db"], [blast_DB], ["-query - -word_size 10 -max_hsps 5 -num_threads 2 -outfmt '6 qseqid sseqid qlen length pident qstart qend sstart send evalue bitscore sstrand' -out"], [path_sample_name+"_clipped.blast"], sep=" ", collapse=None)[0]
        sp.call(cmd_call_blast, shell=True)

        if os.path.isfile(path_sample_name+"_clipped.blast") and os.stat(path_sample_name+"_clipped.blast").st_size == 0:
            os.remove(path_sample_name+"_clipped.blast") # remove blast output file if it is empty

        if os.path.isfile(path_sample_name+"_clipped.blast"):
            sp.call(paste(["sed -i"], [sed_para], [path_sample_name+"_clipped.blast"], sep=" ", collapse=None)[0], shell=True)


        else:
            print(sample_name+" completed: no blast hits found!\n")
            with open(os.path.join(out_folder,"progress.track"),'a') as p:
                p.write("{} completed: no blast hits found!\n".format(sample_name))
            sys.exit()

        if highQ_mate:
            cmd_call_R_processor = paste(["Rscript /tscc/nfs/home/yzhu/projects/gencDNA/scripts/Merge_blast_X_sam_UTR_discor_pairs.R -i"], [path_sample_name+"_clipped.blast"], ["-s"], [path_sample_name+"_merged_UTR_discord_mates_noheader.temp.sam"], ["-a"], [qname_anno_table], ["-o"], [path_sample_name+"_blast_vs_sam.highQual.raw"], sep=" ", collapse=None)[0]

            sp.call(cmd_call_R_processor, shell=True)

        if lowQ_mate:
            cmd_call_R_processor = paste(["Rscript /tscc/nfs/home/yzhu/projects/gencDNA/scripts/Merge_blast_X_sam_UTR_discor_pairs_lowQual_polyA.R -i"], [path_sample_name+"_clipped.blast"], ["-s"], [path_sample_name+"_merged_UTR_lowQUAL_mates_noheader.temp.sam"], ["-a"], [qname_anno_table], ["-o"], [path_sample_name+"_blast_vs_sam.lowQual.raw"], sep=" ", collapse=None)[0]
            sp.call(cmd_call_R_processor, shell=True)


        print(sample_name+" completed!\n")
        with open(os.path.join(out_folder,"progress.track"),'a') as p:
            p.write("{} completed!\n".format(sample_name))
    else:
        print(sample_name+" completed: no confident mate reads found!\n")

        with open(os.path.join(out_folder,"progress.track"),'a') as p:
            p.write("{} completed: no confident mate reads found!\n".format(sample_name))


if __name__ == "__main__":
    main()