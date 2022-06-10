# Generate output for comparison
out.csv: ../target/release/libpyaerso.so mxs.py
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py out.csv

# Rust-backed python version
test: ../target/release/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py

# Python version with re-implemented cross-product
test-python: ../pyaerso_py
	-rm pyaerso.so
	ln -fsT $< pyaerso
	python3 mxs.py

# Python version with numpy's cross-product
test-python-np: ../pyaerso_py
	-rm pyaerso.so
	ln -fsT $< pyaerso
	NUMPY_XPROD=true python3 mxs.py

# Rust-back version (debug build)
test-debug: ../target/debug/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py


SOURCES:=../pyaerso/src/* ../aerso/src/* ../pyaerso/Cargo.toml ../aerso/Cargo.toml
../target/debug/libpyaerso.so: $(SOURCES)
	cd .. && cargo build

../target/release/libpyaerso.so: $(SOURCES)
	cd .. && cargo build --release

.PHONY: test test-python test-mxs test-debug
