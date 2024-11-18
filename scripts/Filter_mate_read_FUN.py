#!/usr/bin/env python

#updated 3/10/2020
import pysam
import re
import sys
import os

def check_clipping_and_edit_disc(rec, cutoff):
    #check if num(hard+soft clipping + edit_distance) < read_length * cutoff
    #return TRUE if yes otherwise False
    readl = rec.infer_read_length() #get read length from CIGAR
    nm = rec.get_tag("NM")
    match = re.findall(r'(\d+)[SH]', rec.cigarstring)
    if len(match) == 0:
        return(float(nm)/readl <= cutoff)
    else:
        num_clip = sum([int(i) for i in match])
        out = (float(num_clip) + float(nm))/readl <= cutoff
        return(out)

def FilterMate(insam, outsam, outsam_lowQUAL, chrom, start, end, cutoff):
    chrs = set(["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"])
    count = 0
    cout_lowq = 0
    input = pysam.AlignmentFile(insam, 'r')
    filtered = pysam.AlignmentFile(outsam, 'w', template=input)
    lowqual = pysam.AlignmentFile(outsam_lowQUAL, 'w', template=input)

    for rec in input.fetch(until_eof=True):
        if rec.is_unmapped:
            cout_lowq += 1
            lowqual.write(rec)

        if not rec.is_unmapped:
            if rec.reference_name in chrs:
                if rec.reference_name != chrom:
                    if check_clipping_and_edit_disc(rec, cutoff):
                        count += 1
                        filtered.write(rec)
                    else:
                        cout_lowq += 1
                        lowqual.write(rec)

                elif rec.reference_name == chrom and (rec.reference_start + rec.infer_read_length() < start  or rec.reference_start > end): # currently not considering reads that's inserted back to the same gene
                    if check_clipping_and_edit_disc(rec, cutoff):
                        count += 1
                        filtered.write(rec)
                    else:
                        cout_lowq += 1
                        lowqual.write(rec)
    input.close()
    filtered.close()
    lowqual.close()
    if count == 0: os.remove(outsam)
    if cout_lowq == 0: os.remove(outsam_lowQUAL)
















