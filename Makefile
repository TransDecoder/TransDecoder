SHELL := /bin/bash

all:
	cd ./transdecoder_plugins/ && $(MAKE) all

clean:
	cd ./transdecoder_plugins/ && $(MAKE) clean
	cd ./sample_data && ./cleanme.pl

test:
	cd ./sample_data/ && ./runMe.sh


