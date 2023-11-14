##
# TOV
#
# @file
# @version 0.1

.PHONY: init-files
init-files:
	mkdir out -p
	touch out/p0MR1.csv
	touch out/p0MR2.csv
	touch out/MR1.csv
	touch out/MR2.csv

.PHONY: enter-repl
enter-repl:
	julia --threads 12 --project=.

# end
