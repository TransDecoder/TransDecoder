## building the PWM:

The start codons and flanking sequence are extracted from the non-redundant long-orf file, and these are used to construct the (+) PWM.

All start codons found downstream of the (+) selected start (and regardless of frame) are extracted to represent the null model, and the (-) PWM is constructed from them.

Outputs:   'atg.+.pwm', 'atg.-.pwm' correspond to the initial PWMs for the positive and null models, and there's a corresponding .features file that contains the sequence features extracted and used to construct the PWMs.

By default, 20 bases upstream and 10 bases downstream from the start codon are extracted, but this is configurable above.



##  Feature scoring, computing ROC and AUC

For scoring, the +.features and -.features are taken and go through 5 rounds of the following:

For each round, 75% of features are selected for building a PWM, and the remaining 25% of features are used for scoring the PWM.   This is done separately for the (+) and the (-) model.  The range of bases within the PWM used for scoring is explored across all combinations of extensions upstream and downstream from the start codon (ATG) within the length of the PWM.  Scores are written to an output file '.scores'.   (This is sort of like cross-validation, but not exactly....)

ROC curve:  True positives, false positives, true negatives, and false negatives are scored according to minimum score thresholds.  Ten uniformly distributed minimum score threshold intervals are chosen based on the min and max of scores.   A '.roc' file is generated including the minimum threshold used, the number of TP, FP, TN, FN, etc.   Note, although we do several rounds of selection (cross-validation style) above in generating 'truth' and 'false' sets that are scored as TP, FP, etc.,  all the resulting TPs, etc.,  are grouped together and treated as a single resulting data set when computing the true positive and false positive rates (for the ROC).  I don't do any variance calculation here....  just assume that the several rounds of selection & scoring will balance out any biases.

The area under the ROC curve (AUC) is computed, generating a '.auc' file.
