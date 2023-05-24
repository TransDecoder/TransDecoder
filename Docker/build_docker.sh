#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

rm -f ./*simg


docker build -t trinityrnaseq/transdecoder:$VERSION .
docker build -t trinityrnaseq/transdecoder:latest .

