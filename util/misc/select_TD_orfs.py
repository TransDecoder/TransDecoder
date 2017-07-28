#!/usr/bin/env python

import sys, os, re


usage = "\n\n\tusage: {} longest_orfs.cds longest_orfs.cds.scored\n\n"
if len(sys.argv) < 3:
    sys.stderr.write(usage)
    sys.exit(1)

long_orfs_cds_file = sys.argv[1]
long_orfs_scored_file = sys.argv[2]

def main():

    predicted_orf_coords = retrieve_orf_coords(long_orfs_cds_file)

    select(long_orfs_scored_file, long_orfs_scored_file + ".markov_val_any.gff", predicted_orf_coords)

    select(long_orfs_scored_file, long_orfs_scored_file + ".markov_2.gff", predicted_orf_coords, markov_val="2")

    select(long_orfs_scored_file, long_orfs_scored_file + ".markov_3.gff", predicted_orf_coords, markov_val="3")

    select(long_orfs_scored_file, long_orfs_scored_file + ".markov_4.gff", predicted_orf_coords, markov_val="4")

    select(long_orfs_scored_file, long_orfs_scored_file + ".markov_5.gff", predicted_orf_coords, markov_val="5")

    
    
    sys.exit(0)



def retrieve_orf_coords(long_orfs_cds_file):

    orf_acc_to_coord_info = {}

    with open(long_orfs_cds_file) as fh:
        for line in fh:
            if not re.search("^>", line):
                continue
            match = re.search(" (\S+):(\d+)-(\d+)\(([+-])\)$", line)
            if not match:
                raise RuntimeError("Error, cannot extract orf coord info from line: {}".format(line))
            (transcript_id, lend, rend, orient) = (match.group(1), match.group(2), match.group(3), match.group(4))
            (lend, rend) = sorted((lend, rend))

            match = re.search("^>(\S+)", line)
            orf_id = match.group(1)

            orf_acc_to_coord_info[orf_id] = {
                'transcript_id' : transcript_id,
                'lend' : lend,
                'rend' : rend,
                'orient' : orient }


    return orf_acc_to_coord_info



def select(input_file, output_file, predicted_orf_coords,
           markov_val=None,
           require_pos_F1 = True,
           require_F1_max = True):

    sys.stderr.write("-writing {}\n".format(output_file))

    fh = open(input_file)
    ofh = open(output_file, 'w')

    seen_orf_ids = set()

    for line in fh:
        if re.search("^#", line): continue
        line = line.rstrip()
        (orf_id, markov_order, seq_length, score_1, score_2, score_3, score_4, score_5, score_6)  = line.split("\t")

        if orf_id in seen_orf_ids:
            continue

        ## apply filtering criteria
        
        if markov_val and markov_val != markov_order:
            continue

        if require_pos_F1 and score_1 <= 0:
            continue

        if require_F1_max:
            nonF1_scores = (score_2, score_3, score_4, score_5, score_6)
            max_nonF1_score = max(nonF1_scores)
            if max_nonF1_score > score_1:
                continue

        ## passed filters, report it
        orf_struct = predicted_orf_coords[orf_id]
        
        ofh.write("\t".join([orf_struct['transcript_id'], "selected", "CDS",
                        orf_struct['lend'], orf_struct['rend'], '.', orf_struct['orient'], '.', '.']) + "\n")
        
        seen_orf_ids.add(orf_id)
        

    return




if __name__ == '__main__':
    main()
