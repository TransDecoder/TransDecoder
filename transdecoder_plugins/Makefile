SHELL := /bin/bash

CDHIT="cd-hit-v4.6.1-2012-08-27"

all: cdhit


cdhit:
	tar -xvf ${CDHIT}.tgz \
	&& cd ${CDHIT} \
    && $(MAKE) openmp=yes && $(MAKE) install PREFIX=../../util/bin
	mv ${CDHIT} cdhit


clean:
	rm -rf ./cdhit
	rm -f ../util/bin/*
