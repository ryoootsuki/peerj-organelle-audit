SHELL := /usr/bin/env bash
RUN := ./run

.PHONY: plan autotune test doctor local-check references pilot full silene chicken audit downstream report status readiness performance clean-tmp clean-intermediate

plan:
	$(RUN) plan both
autotune:
	$(RUN) autotune both
test:
	python tests/offline_smoke.py
	python tests/results_numeric_regression.py
	python tests/public_safety_audit.py
doctor:
	$(RUN) doctor both
local-check:
	$(RUN) local-check both
references:
	$(RUN) references both
pilot:
	$(RUN) pilot both
full:
	$(RUN) full both
silene:
	$(RUN) full silene
chicken:
	$(RUN) full chicken
audit:
	$(RUN) audit both
downstream:
	$(RUN) downstream both
report:
	$(RUN) report both
status:
	$(RUN) status both
readiness:
	$(RUN) readiness both
performance:
	$(RUN) performance both
clean-tmp:
	$(RUN) clean both --clean-level tmp
clean-intermediate:
	$(RUN) clean both --clean-level intermediate
