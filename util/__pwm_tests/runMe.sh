#!/bin/bash

set -ev

gunzip -c longest_orfs.cds.top_longest_5000.nr80.gz > test.train_orfs.fa

gunzip -c pasa_assemblies.fasta.gz > test.transcripts.fa

../PWM/build_atgPWM_+-.pl --transcripts test.transcripts.fa --selected_orfs test.train_orfs.fa --out_prefix test

../PWM/feature_scoring.+-.pl --features_plus test.+.features --features_minus test.-.features > test.features.scores

../PWM/feature_scores_to_ROC.pl test.features.scores > test.features.scores.ROC

../PWM/plot_ROC.Rscript test.features.scores.ROC

../PWM/compute_AUC.Rscript test.features.scores.ROC


