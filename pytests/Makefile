test: ../target/release/libpyaerso.so
	ln -fs $< pyaerso.so
	python3 test.py

test-debug: ../target/debug/libpyaerso.so
	ln -fs $< pyaerso.so
	python3 test.py

../target/debug/libpyaerso.so: ../pyaerso/src/*
	cd .. && cargo build

../target/release/libpyaerso.so: ../pyaerso/src/*
	cd .. && cargo build --release

.PHONY: test test-debug