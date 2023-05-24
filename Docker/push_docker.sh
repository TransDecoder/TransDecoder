#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker push trinityrnaseq/transdecoder:$VERSION 
docker push trinityrnaseq/transdecoder:latest

