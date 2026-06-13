# TIC — Tournament Image Codec
# Lossless RGB compression: G-RG-BG + MED + zlib

CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lz

# If system zlib not available, use bundled:
# CFLAGS += -I../zlib-1.3.1
# LDFLAGS = ../zlib-1.3.1/libz.a

TARGET = tic_cli

all: $(TARGET)

$(TARGET): tic.c tic.h
	$(CC) $(CFLAGS) -o $@ tic.c $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

test: $(TARGET)
	@echo "=== Roundtrip test ==="
	@python3 -c "\
import numpy as np; \
rgb = np.random.randint(0, 256, (128, 128, 3), dtype=np.uint8); \
rgb.tofile('/tmp/tic_test.raw'); \
print('Generated 128x128 test image')"
	./$(TARGET) compress /tmp/tic_test.raw 128 128 /tmp/tic_test.tic
	./$(TARGET) decompress /tmp/tic_test.tic /tmp/tic_test_dec.raw
	@python3 -c "\
import numpy as np; \
a = np.fromfile('/tmp/tic_test.raw', dtype=np.uint8); \
b = np.fromfile('/tmp/tic_test_dec.raw', dtype=np.uint8); \
print('Roundtrip:', 'PASS' if np.array_equal(a, b) else 'FAIL')"

bench: $(TARGET)
	@echo "=== Speed benchmark ==="
	./$(TARGET) bench /tmp/tic_test.raw 128 128 200

.PHONY: all clean test bench
