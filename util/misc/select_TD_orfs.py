#!/usr/bin/env python

import sys, os, re


usage = "\n\n\tusage: {} longest_orfs.cds longest_orfs.cds.scores\n\n"
if len(sys.argv) < 3:
    sys.stderr.write(usage)
    sys.exit(1)

long_orfs_cds_file = sys.argv[1]
long_orfs_scored_file = sys.argv[2]

ORIG_TRANSDECODER_FLAG = False


def main():

    predicted_orf_coords = retrieve_orf_coords(long_orfs_cds_file)

    prediction_list = parse_predictions_and_scores(long_orfs_scored_file, predicted_orf_coords)


    def default_frame_analysis_func(frame_scores):
        frame_1_score = frame_scores[0]
        if frame_1_score > 0 and frame_1_score > max(frame_scores[1:]):
            return True
        else:
            return False

    #select(prediction_list, long_orfs_scored_file + ".def_all_best_orfs.gff", predicted_orf_coords, default_frame_analysis_func)

    #select(prediction_list, long_orfs_scored_file + ".def_single_best_orf.gff", predicted_orf_coords, default_frame_analysis_func, longest_single_orf=True)

    #select(prediction_list, long_orfs_scored_file + ".def_single_best_orf_c700.gff", predicted_orf_coords, default_frame_analysis_func, longest_single_orf=True, capture_long_orfs_size=700)
    #select(prediction_list, long_orfs_scored_file + ".def_single_best_orf_c1000.gff", predicted_orf_coords, default_frame_analysis_func, longest_single_orf=True, capture_long_orfs_size=1000)

    select_custom(prediction_list, long_orfs_scored_file + ".def_single_custom.gff", predicted_orf_coords)
    

    
    
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


def select(prediction_list, output_file, predicted_orf_coords,
           frame_analysis_func,
           markov_val=None,
           longest_single_orf = False,
           capture_long_orfs_size=-1):

    sys.stderr.write("-writing {}\n".format(output_file))

    
    ofh = open(output_file, 'w')

    seen_orf_ids = set() 
    seen_transcript_ids = set()

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

        
        if orf_id in seen_orf_ids:
            # in case already selected under a different Markov model
            continue

        if longest_single_orf and transcript_id in seen_transcript_ids:
            continue

        ###########################
        ## apply filtering criteria
        pass_orf = True
        
        if pass_orf and \
               markov_val is not None and \
               markov_val != markov_order:

            pass_orf = False


        if pass_orf and not frame_analysis_func(frame_scores):
            pass_orf = False
        

        if capture_long_orfs_size > 0 and orf_length >= capture_long_orfs_size:
            # free pass
            pass_orf = True
        
        if pass_orf:
            ## passed filters, report it
        
            ofh.write("\t".join([transcript_id, "selected", "CDS",
                                 str(lend), str(rend), '.',
                                 orient, '.', orf_id]) + "\n")
            
            seen_orf_ids.add(orf_id)
            seen_transcript_ids.add(transcript_id)

    return


def select_custom(prediction_list, output_file, predicted_orf_coords):

    long_orf_size = 1000
    moderate_orf_size = 500
    
    sys.stderr.write("-writing {}\n".format(output_file))
    
    ofh = open(output_file, 'w')

    seen_orf_ids = set() 
    seen_transcript_ids = set()

    longest_single_orf = True


    import collections
    
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

        
        if longest_single_orf and transcript_id in seen_transcript_ids:
            continue

        ###########################
        ## apply filtering criteria
        pass_orf = True

        fst_is_max = (frame_scores[0] > max(frame_scores[1:]))
        fst_gt_zero = (frame_scores[0] > 0)


        if ORIG_TRANSDECODER_FLAG:
            # standard alg
            
            if not (fst_is_max and fst_gt_zero):
                pass_orf = False
        else:
            
            
            fst_max3 = (frame_scores[0] > max(frame_scores[1:3]))

            if orf_length >= long_orf_size:
                pass_orf = True
            elif orf_length >= moderate_orf_size:
                if not (fst_max3 and fst_gt_zero):
                    pass_orf = False
            else:
                if not (fst_max3 and fst_gt_zero):
                    pass_orf = False

        if pass_orf:

            transcript_to_selected_orfs[transcript_id].append(prediction)






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

        top_selected_pred = selected_preds[0]
        orf_struct = top_selected_pred['orf_struct']
        
        ofh.write("\t".join([transcript_id, "selected", "CDS",
                             str(orf_struct['lend']), str(orf_struct['rend']), '.',
                             orf_struct['orient'], str(top_selected_pred['frame_scores'][0]), top_selected_pred['orf_id']]) + "\n")
            
           



    return




if __name__ == '__main__':
    main()
