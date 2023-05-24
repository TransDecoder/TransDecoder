#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build transdecoder.v${VERSION}.simg docker://trinityrnaseq/transdecoder:$VERSION

singularity exec -e transdecoder.v${VERSION}.simg TransDecoder.LongOrfs

ln -sf  transdecoder.v${VERSION}.simg  transdecoder.simg

