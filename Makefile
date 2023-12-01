##
# TOV
#
# @file
# @version 0.1

.PHONY: init-files
init-files:
	mkdir out -p
	touch out/seqtest_p0MR1.csv
	touch out/seqtest_p0MR2.csv
	touch out/seqtest_MR1.csv
	touch out/seqtest_MR2.csv
	touch out/deftest_MR.csv

.PHONY: enter-repl
enter-repl:
	julia --threads 12 --project=.

# end
