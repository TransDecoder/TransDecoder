#!/usr/bin/env python

import os, re, sys

usage = "\n\n\tusage: {} longest_orfs.gff3 longest_orfs.cds.scores\n\n".format(sys.argv[0])

if len(sys.argv) < 3:
    sys.stderr.write(usage)
    sys.exit(1)


longest_orfs_gff3 = sys.argv[1]
longest_orfs_scores = sys.argv[2]


def main():

    ref_orfs_dict = define_ref_orfs(longest_orfs_gff3)

    
    with open(longest_orfs_scores) as fh:
        header = fh.next()
        header = header.rstrip()
        print(header + "\treforf")
        for line in fh:
            line = line.rstrip()
            x = line.split("\t")
            orf_name = x[0]
            class_val = "NO"
            if orf_name in ref_orfs_dict:
                class_val = "YES"
            print(line + "\t" + class_val)


    sys.exit(0)
    
                          
            
            


def define_ref_orfs(gff3_file):
    ref_orfs_dict = {}

    with open(gff3_file) as fh:
        for line in fh:
            x = line.split("\t")
            if len(x) < 8:
                continue
            trans_name = x[0]
            feat_type = x[2]
            lend = x[3]
            rend = x[4]
            info = x[8]
            if feat_type != 'CDS':
                continue
            (ref_lend, ref_rend) = extract_ref_coords(trans_name)
            

            if ref_rend == rend:
                orf_name = extract_orf_name(info)
                ref_orfs_dict[orf_name] = True

    return ref_orfs_dict


def extract_ref_coords(trans_name):

    m = re.search("\|CDS:(\d+)-(\d+)\|", trans_name)
    if not m:
        raise RuntimeError("Not able to extract cds coords from trans name: {}".format(trans_name))

    lend = m.group(1)
    rend = m.group(2)

    return(lend, rend)


def extract_orf_name(info):

    m = re.search("Parent=([^;\s]+)", info)
    if not m:
        raise RuntimeError("Error, cannot extract orf ID from info: {}".format(info))

    orf_name = m.group(1)
    return orf_name


if __name__ == '__main__':
    main()
