#!/bin/bash

set -ev

gunzip -c longest_orfs.cds.top_longest_5000.nr80.gz > test.train_orfs.fa

gunzip -c pasa_assemblies.fasta.gz > test.transcripts.fa

../train_start_PWM.pl --transcripts test.transcripts.fa --selected_orfs test.train_orfs.fa --out_prefix test

