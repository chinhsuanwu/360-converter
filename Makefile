OUTPUT_DIR = out

all:
	@clang++ example/example.cpp -o example/example

run:
	@if [ ! -d $(OUTPUT_DIR) ]; then mkdir $(OUTPUT_DIR); fi
	@./example/example assets/front.png assets/right.png assets/back.png assets/left.png assets/top.png assets/down.png

.PHONY: clean
clean:
	@rm example/example
