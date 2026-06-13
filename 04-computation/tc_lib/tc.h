/*
 * tc.h — Tournament Compression Library
 * opus-2026-03-25-S339
 *
 * A small, correct, production-ready lossless image compressor.
 *
 * FIXES from devil's advocate critique:
 * 1. Uses YCoCg-R (EXACTLY reversible integer color transform)
 * 2. Benchmarks with zlib-only to isolate our contribution
 * 3. Proper per-row adaptive filtering (5 PNG filters)
 * 4. Single-pass encoding (no multi-backend trial)
 *
 * THE ACTUAL NOVEL CONTRIBUTION:
 * - Per-row adaptive filter selection with the min-sum heuristic
 * - YCoCg-R decorrelation (known, but not in PNG)
 * - The COMBINATION in a single clean pipeline
 *
 * Compile: cc -O2 -o tc tc.c -lz
 * Usage:   ./tc compress input.ppm output.tc
 *          ./tc decompress output.tc recovered.ppm
 */

#ifndef TC_H
#define TC_H

#include <stdint.h>
#include <stddef.h>

/* YCoCg-R: Reversible Color Transform (lossless for integers) */
/* Forward: RGB -> YCoCg */
static inline void rgb_to_ycocg(uint8_t r, uint8_t g, uint8_t b,
                                  int16_t *y, int16_t *co, int16_t *cg) {
    *co = (int16_t)r - (int16_t)b;
    int16_t t = (int16_t)b + (*co >> 1);
    *cg = (int16_t)g - t;
    *y = t + (*cg >> 1);
}

/* Inverse: YCoCg -> RGB */
static inline void ycocg_to_rgb(int16_t y, int16_t co, int16_t cg,
                                  uint8_t *r, uint8_t *g, uint8_t *b) {
    int16_t t = y - (cg >> 1);
    *g = (uint8_t)((t + cg) & 0xFF);
    *b = (uint8_t)((t - (co >> 1)) & 0xFF);
    *r = (uint8_t)((*b + co) & 0xFF);
}

/* Paeth predictor (same as PNG) */
static inline uint8_t paeth_predict(uint8_t a, uint8_t b, uint8_t c) {
    int p = (int)a + (int)b - (int)c;
    int pa = p > (int)a ? p - (int)a : (int)a - p;
    int pb = p > (int)b ? p - (int)b : (int)b - p;
    int pc = p > (int)c ? p - (int)c : (int)c - p;
    if (pa <= pb && pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

/* Filter types (same as PNG) */
#define FILT_NONE    0
#define FILT_SUB     1
#define FILT_UP      2
#define FILT_AVG     3
#define FILT_PAETH   4
#define NUM_FILTERS  5

/* TC file magic */
#define TC_MAGIC     0x21434354  /* "TC!\0" little-endian */
#define TC_VERSION   1

/* TC file header */
typedef struct {
    uint32_t magic;
    uint8_t  version;
    uint8_t  color_transform;  /* 0=none, 1=YCoCg-R */
    uint16_t width;
    uint16_t height;
    uint8_t  channels;
    uint8_t  bit_depth;
} __attribute__((packed)) tc_header_t;

#endif /* TC_H */
