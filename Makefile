.PHONY: enter-repl
enter-repl: init-out
	julia --threads=auto --project=. -i

.PHONY: init-out
init-out:
	mkdir -p test/out