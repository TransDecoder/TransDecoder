#!/usr/bin/env python

import sys, os, re
import collections
    
ORIG_TRANSDECODER_FLAG = False

def main():


    import argparse

    parser = argparse.ArgumentParser(description="transdecoder orf selection algorithm", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--long_orfs_cds", dest="long_orfs_cds_filename", type=str, default="", required=True, help="longest_orfs.cds file")
    parser.add_argument("--long_orfs_scores", dest="long_orfs_scores_filename", type=str, default="", required=True, help="longest_orfs.cds.scores file")
    parser.add_argument("--single_best", dest="single_best_flag", action='store_true', default=False, help="select only single best orf")
    parser.add_argument("--all_good_orfs", dest="all_good_orfs", action='store_true', default=False, help="retain all good orfs allowing overlaps")
    parser.add_argument("--td_orig", dest="td_orig_flag", action='store_true', default=False, help='use the original TD algorithm')

    args = parser.parse_args()
    
    long_orfs_cds_file = args.long_orfs_cds_filename
    long_orfs_scored_file = args.long_orfs_scores_filename
        
    predicted_orf_coords = retrieve_orf_coords(long_orfs_cds_file)

    prediction_list = parse_predictions_and_scores(long_orfs_scored_file, predicted_orf_coords)


    def orig_transdecoder_algorithm(seq_length, frame_scores):
        return (fst_gt_zero(frame_scores) and fst_is_max_all(frame_scores))

    def custom_td_alg(orf_length, frame_scores, long_orf_size=1000, moderate_orf_size=500):

        pass_orf = True
        
        if orf_length >= long_orf_size:
            pass_orf = True # stays true
        elif orf_length >= moderate_orf_size:
            if not (fst_is_max3(frame_scores) and fst_gt_zero(frame_scores)):
                pass_orf = False
        else:
            if not (fst_is_max3(frame_scores) and fst_gt_zero(frame_scores)):
                pass_orf = False

        return pass_orf


    select_alg = custom_td_alg
    if args.td_orig_flag:
        select_alg = orig_transdecoder_algorithm
    

    # all good orfs
    trans_to_preds_list = select(prediction_list, predicted_orf_coords, select_alg)

    if not args.all_good_orfs:
        # by default, we prune overlaps
        trans_to_preds_list = select_best_non_overlapping_orfs(trans_to_preds_list)
        
    elif args.single_best_flag:
        trans_to_preds_list = select_single_orf_per_transcript(trans_to_preds_list)
    

    write_preds_to_file(trans_to_preds_list, sys.stdout)
    
    
    
    sys.exit(0)


def fst_is_max_all(frame_scores):
    return (frame_scores[0] > max(frame_scores[1:]))

def fst_gt_zero(frame_scores):
    return (frame_scores[0] > 0)

def fst_is_max3(frame_scores):
    return (frame_scores[0] > max(frame_scores[1:3]))



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
            (lend, rend) = sorted((int(lend), int(rend)))

            match = re.search("^>(\S+)", line)
            orf_id = match.group(1)

            orf_acc_to_coord_info[orf_id] = {
                'transcript_id' : transcript_id,
                'lend' : lend,
                'rend' : rend,
                'orient' : orient }


    return orf_acc_to_coord_info


def parse_predictions_and_scores(long_orfs_scored_file, predicted_orf_coords):
        
    prediction_list = []
    
    fh = open(long_orfs_scored_file)

    for line in fh:
        if re.search("^#", line): continue
        line = line.rstrip()
        (orf_id, markov_order, orf_length, score_1, score_2, score_3, score_4, score_5, score_6)  = line.split("\t")

        orf_length = int(orf_length)
        score_1 = float(score_1)
        score_2 = float(score_2)
        score_3 = float(score_3)
        score_4 = float(score_4)
        score_5 = float(score_5)
        score_6 = float(score_6)


        orf_struct = predicted_orf_coords[orf_id]

        prediction = { 'orf_id' : orf_id,
                       'markov_order' : markov_order,
                       'orf_length' : orf_length,
                       'frame_scores' : (score_1, score_2, score_3, score_4, score_5, score_6),
                       'orf_struct' : orf_struct
                       }

        prediction_list.append(prediction)



    prediction_list.sort(key=lambda x: x['orf_length'])

    prediction_list.reverse()

    return prediction_list



def select(prediction_list, predicted_orf_coords, selection_alg):
    
    
    transcript_to_selected_orfs = collections.defaultdict(list)

    for prediction in prediction_list:

        orf_id = prediction['orf_id']
        markov_order = prediction['markov_order']
        orf_length = prediction['orf_length']
        frame_scores = prediction['frame_scores']
        (score_1, score_2, score_3, score_4, score_5, score_6) = frame_scores

        orf_struct = prediction['orf_struct']

        transcript_id = orf_struct['transcript_id']
        lend = orf_struct['lend']
        rend = orf_struct['rend']
        orient = orf_struct['orient']


        ###########################
        ## apply filtering criteria
        pass_orf = selection_alg(orf_length, frame_scores)

        if pass_orf:
            transcript_to_selected_orfs[transcript_id].append(prediction)


    
    return transcript_to_selected_orfs



def select_single_orf_per_transcript(transcript_to_selected_orfs):


    non_overlapping_orfs = select_best_non_overlapping_orfs(transcript_to_selected_orfs)

    ret_transcript_to_selected_orfs = collections.defaultdict(list)
    

    for transcript_id in non_overlapping_orfs:

        prediction_list = non_overlapping_orfs[transcript_id]
                
        prediction_list.sort(key=lambda x: x['orf_length'], reverse=True)
        #prediction_list.sort(key=lambda x: x['frame_scores'][0], reverse=True)

        top_selected_pred = prediction_list[0]
        
        ret_transcript_to_selected_orfs[transcript_id].append(top_selected_pred)
        
        
    return ret_transcript_to_selected_orfs


def select_best_non_overlapping_orfs(transcript_to_selected_orfs):


    ret_transcript_to_non_overlapping_orfs = collections.defaultdict(list)
    

    for transcript_id in transcript_to_selected_orfs:

        prediction_list = transcript_to_selected_orfs[transcript_id]

        
        if ORIG_TRANSDECODER_FLAG:
            prediction_list.sort(key=lambda x: x['orf_length'], reverse=True)
        else:
            prediction_list.sort(key=lambda x: x['frame_scores'][0], reverse=True)
        
        selected_preds = []
        
        for prediction in prediction_list:
        
            ## passed filters, report it

                                
            # ensure it doesn't overlap a selected one.
            found_overlap = False
            for pred in selected_preds:
                if prediction['orf_struct']['lend'] < pred['orf_struct']['rend'] and \
                   prediction['orf_struct']['rend'] > pred['orf_struct']['lend']:
                    found_overlap = True
                    break
            
            if not found_overlap:
                selected_preds.append(prediction)
        
        selected_preds.sort(key=lambda x: x['orf_length'], reverse=True)

        ret_transcript_to_non_overlapping_orfs[transcript_id].extend(selected_preds)


    return ret_transcript_to_non_overlapping_orfs



def write_preds_to_file(transcript_to_selected_orfs, ofh):


    for transcript_id in transcript_to_selected_orfs:

        pred_list = transcript_to_selected_orfs[transcript_id]

        for pred in pred_list:
            orf_struct = pred['orf_struct']
                        
            ofh.write("\t".join([transcript_id, "selected", "CDS",
                                 str(orf_struct['lend']), str(orf_struct['rend']), '.',
                                 orf_struct['orient'], str(pred['frame_scores'][0]), pred['orf_id']]) + "\n")
            
           



    return




if __name__ == '__main__':

    main()
    
