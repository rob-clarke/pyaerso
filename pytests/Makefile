compare: out_rs.csv out_py.csv
	python ../../../plot_compare.py out_rs.csv out_py.csv Rust Python

plot_rs: out_rs.csv
	python ../../../plot_compare.py out_rs.csv out_rs.csv Rust Rust

plot_bix: out_bixler.csv
	python ../../../plot_compare.py out_bixler.csv out_bixler.csv Rust Rust

out_bixler.csv: ../target/release/libpyaerso.so bixler.py
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 bixler.py out_bixler.csv

# Generate output for comparison
out_rs.csv: ../target/release/libpyaerso.so mxs.py
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py out_rs.csv

# Generate Python output
out_py.csv: ../pyaerso_py/** mxs.py
	-rm pyaerso.so
	ln -fsT ../pyaerso_py pyaerso
	python3 mxs.py out_py.csv

# Rust-backed python version
test: ../target/release/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py

# Python version with re-implemented cross-product
test-python: ../pyaerso_py/**
	-rm pyaerso.so
	ln -fsT ../pyaerso_py pyaerso
	python3 mxs.py

# Python version with numpy's cross-product
test-python-np: ../pyaerso_py
	-rm pyaerso.so
	ln -fsT ../pyaerso_py pyaerso
	NUMPY_XPROD=true python3 mxs.py

# Rust-back version (debug build)
test-debug: ../target/debug/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py


SOURCES:=../pyaerso/src/** ../aerso/src/** ../pyaerso/Cargo.toml ../aerso/Cargo.toml
../target/debug/libpyaerso.so: $(SOURCES)
	cd .. && cargo build

../target/release/libpyaerso.so: $(SOURCES)
	cd .. && cargo build --release

.PHONY: test test-python test-mxs test-debug
