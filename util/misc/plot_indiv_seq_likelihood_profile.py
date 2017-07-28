#!/usr/bin/env python

import os,sys
import re
import matplotlib.pyplot as plt
import argparse
import subprocess

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description="plot likelihood profile for sequence ")

parser.add_argument("--orf_id", type=str, required=True, help="orf accession")

parser.add_argument("--longest_orfs_cds", type=str, required=True, help="long orfs cds file")

parser.add_argument("--kmer_scores", type=str, required=True, help= "kmer likelihood score file")

args = parser.parse_args()



def main():
    
    seq = get_seq(args.orf_id, args.longest_orfs_cds)

    framed_kmers_to_likelihoods = parse_kmer_likelihoods(args.kmer_scores)
    

    score_vec = score_seq(seq, framed_kmers_to_likelihoods)

    plt.plot(range(1,len(score_vec)+1), score_vec, marker ='+')
    plt.show()


def score_seq(seq, framed_kmer_likelihoods):

    score_vec = []

    for i in range(1, len(seq)):
        frame = i % 3
        markov_use = min(i, 5)
        kmer = seq[i-markov_use:i+1]
        #print("i:{}, markov_use:{}, kmer:{}".format(i, markov_use, kmer))

        framed_kmer = "{}-{}".format(kmer, frame)

        loglikelihood = framed_kmer_likelihoods[framed_kmer]

        #print("i:{}, likelihood: {}".format(i, loglikelihood))
        score_vec.append(loglikelihood)

    return score_vec
    



def parse_kmer_likelihoods(kmer_scores_file):

    framed_kmers_to_likelihoods = {}

    with open(kmer_scores_file) as fh:
        for line in fh:
            if re.search("^#", line): continue
            line = line.rstrip()
            (framed_kmer, count, countkmerminus1, likelihood) = line.split("\t")
            framed_kmers_to_likelihoods[framed_kmer] = float(likelihood)


    return framed_kmers_to_likelihoods



def get_seq(orf_id, fasta_file):

    cmd = "samtools faidx {} \"{}\"".format(fasta_file, orf_id)

    fasta_entry = subprocess.check_output(cmd, shell=True)
    print(fasta_entry)

    lines = fasta_entry.split("\n")
    
    header = lines.pop(0)

    seq = "".join(lines)
    seq = seq.replace(" ", "")

    return seq






if __name__ == '__main__':
    main()

