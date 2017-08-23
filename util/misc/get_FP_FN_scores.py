#!/usr/bin/env python

import os, sys, re

usage = "usage: {} longest_orfs.cds.scores longest_orfs.cds cds.gff.scored\n\n".format(sys.argv[0])

if len(sys.argv) < 4:
    sys.stderr.write(usage)
    sys.exit(1)



def main():

    cds_scores_file = sys.argv[1]
    cds_fasta = sys.argv[2]
    scored_orfs_file = sys.argv[3]

    # link orf name to transcript and 3' coordinate
    trans_end_to_orf = {}
    with open(cds_fasta) as fh:
        for line in fh:
            if re.search(">", line):
                line = line[1:]
                m = re.search("^(\S+).*:(\d+)-(\d+)\([+-]\)", line)
                if m:
                    orf_name = m.group(1)
                    lend = int(m.group(2))
                    rend = int(m.group(3))
                    (lend, rend) = sorted((lend, rend))

                    pts = re.split("[\|\:]+", orf_name)
                    trans_id = pts[1]
                    trans_token = trans_id + "|" + str(rend)
                    trans_end_to_orf[trans_token] = orf_name
                    #print("{} -> {}".format(trans_token, orf_name))
                else:
                    raise RuntimeError("Error, no regex match for extracting model coords from : {}".format(line))


    # get scores info for orf
    orf_scores = {}
    with open(cds_scores_file) as fh:
        for line in fh:
            line = line.rstrip()
            x = line.split("\t")
            orf_name = x[0]
            len_n_scores = x[2:]
            orf_scores[orf_name] = len_n_scores


    # write FP FN score info
    with open(scored_orfs_file) as fh:
        for line in fh:
            line = line.rstrip()
            x = line.split("\t")
            c = x[0]
            pred_rend = x[7]
            if c not in ("FP", "FN"):
                continue
            
            trans_name = x[1]
            if c == "FN":
                pred_rend = x[3]

            token = trans_name + "|" + pred_rend
            if token in trans_end_to_orf:
                pred_orf_name = trans_end_to_orf[token]
                scores = orf_scores[pred_orf_name]
                x += scores
            else:
                x.append("__MISSING_INFO__")
            
            print("\t".join(x))



    sys.exit(0)




if __name__ == '__main__':
    main()
