out.csv: ../target/release/libpyaerso.so mxs.py
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py out.csv

test: ../target/release/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 test.py

test-python: ../pyaerso_py
	-rm pyaerso.so
	ln -fsT $< pyaerso
	python3 mxs.py

test-mxs: ../target/release/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 mxs.py

test-debug: ../target/debug/libpyaerso.so
	-rm pyaerso
	ln -fs $< pyaerso.so
	python3 test.py


SOURCES:=../pyaerso/src/* ../aerso/src/* ../pyaerso/Cargo.toml ../aerso/Cargo.toml
../target/debug/libpyaerso.so: $(SOURCES)
	cd .. && cargo build

../target/release/libpyaerso.so: $(SOURCES)
	cd .. && cargo build --release

.PHONY: test test-python test-mxs test-debug
