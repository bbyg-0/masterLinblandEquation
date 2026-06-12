CC = gcc
CFLAGS = -Wall -O2 -std=c99 -lm
TARGET = lindblad_sim
SOURCES = main.c simple_tensor.c lindblad_config.c matrix_operations.c
HEADERS = simple_tensor.h lindblad_config.h matrix_operations.h

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES)

run: $(TARGET)
	./$(TARGET) config.json density_matrix.dat expectations.dat

run_txt: $(TARGET)
	./$(TARGET) config.txt density_matrix.dat expectations.dat

clean:
	rm -f $(TARGET) density_matrix.dat expectations.dat

.PHONY: all clean run run_txt
