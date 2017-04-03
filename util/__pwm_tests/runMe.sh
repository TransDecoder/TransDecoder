#!/bin/bash

set -ev

../PWM/simulate_feature_seq_from_PWM.pl  --pwm atg.-.pwm --num_features 1000 > rand.features
../PWM/build_pwm.pl  --features rand.features --pwm_out rand.pwm
../PWM/optimize_atgPWM_+-.predict.pl --features rand.features --atg_pwm_minus rand.pwm 
xpdf rand.features.+-.pwm_range_scores.boxplot.pdf
