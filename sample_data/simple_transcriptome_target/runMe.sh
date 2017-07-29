#!/bin/bash

if [ ! -e Trinity.fasta ]; then
    gunzip -c Trinity.fasta.gz > Trinity.fasta
fi

../../TransDecoder.LongOrfs -t Trinity.fasta 

../../TransDecoder.Predict -t Trinity.fasta  $*

exit 0
