
DIRS=cufflinks_example pasa_example simple_transcriptome_target stringtie_example supertranscripts_example

test:
	@for i in $(DIRS); do \
	echo "Running example in $$$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done


test_docker:
	cd simple_transcriptome_target && ./runMe.Docker.sh

clean:
	@for i in $(DIRS); do \
	echo "Running example in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done

