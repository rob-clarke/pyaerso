test: pyaerso.so
	python3 test.py

pyaerso.so: ../pyaerso/src/*
	cd .. && cargo build
	
.PHONY: test