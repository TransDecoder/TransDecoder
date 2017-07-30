SHELL := /bin/bash

all:
	@echo Nothing to build.  Run \'make test\' to run example executions on each of the sample data sets provided.

clean:
	cd ./sample_data && make clean

test:
	cd ./sample_data/ && make test


