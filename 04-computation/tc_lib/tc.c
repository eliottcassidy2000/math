/*
 * tc.c — Tournament Compression Library: core implementation
 * opus-2026-03-25-S339
 *
 * A correct, fast lossless image compressor in C.
 *
 * Pipeline: YCoCg-R decorrelation → per-row adaptive filtering → zlib
 *
 * Compile: cc -O2 -o tc tc.c -lz
 * Usage:   ./tc compress input.ppm output.tc
 *          ./tc decompress output.tc recovered.ppm
 *          ./tc benchmark input.ppm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "tc.h"

/* ================================================================
 * FILTERING: per-row adaptive (5 PNG filters, pick best per row)
 * ================================================================ */

static void apply_filter(uint8_t filt, const uint8_t *row, const uint8_t *above,
                          uint8_t *out, int w, int bpp) {
    int j;
    for (j = 0; j < w; j++) {
        uint8_t a = (j >= bpp) ? row[j - bpp] : 0;
        uint8_t b = above[j];
        uint8_t c = (j >= bpp) ? above[j - bpp] : 0;
        switch (filt) {
            case FILT_NONE:  out[j] = row[j]; break;
            case FILT_SUB:   out[j] = row[j] - a; break;
            case FILT_UP:    out[j] = row[j] - b; break;
            case FILT_AVG:   out[j] = row[j] - ((a + b) >> 1); break;
            case FILT_PAETH: out[j] = row[j] - paeth_predict(a, b, c); break;
        }
    }
}

static void undo_filter(uint8_t filt, const uint8_t *filtered, const uint8_t *above,
                          uint8_t *out, int w, int bpp) {
    int j;
    for (j = 0; j < w; j++) {
        uint8_t a = (j >= bpp) ? out[j - bpp] : 0;
        uint8_t b = above[j];
        uint8_t c = (j >= bpp) ? above[j - bpp] : 0;
        switch (filt) {
            case FILT_NONE:  out[j] = filtered[j]; break;
            case FILT_SUB:   out[j] = filtered[j] + a; break;
            case FILT_UP:    out[j] = filtered[j] + b; break;
            case FILT_AVG:   out[j] = filtered[j] + ((a + b) >> 1); break;
            case FILT_PAETH: out[j] = filtered[j] + paeth_predict(a, b, c); break;
        }
    }
}

/* Pick the filter that minimizes sum of min(byte, 256-byte) */
static uint8_t best_filter(const uint8_t *row, const uint8_t *above, int w, int bpp) {
    uint8_t temp[w];
    uint8_t best = 0;
    long best_score = (long)w * 256;

    for (uint8_t f = 0; f < NUM_FILTERS; f++) {
        apply_filter(f, row, above, temp, w, bpp);
        long score = 0;
        for (int j = 0; j < w; j++) {
            int v = temp[j];
            score += (v < 128) ? v : (256 - v);
        }
        if (score < best_score) {
            best_score = score;
            best = f;
        }
    }
    return best;
}

/* ================================================================
 * COMPRESS
 * ================================================================ */

int tc_compress(const uint8_t *pixels, int width, int height, int channels,
                uint8_t **out, size_t *out_len) {
    int row_bytes = width * channels;
    int bpp = channels;

    /* Allocate workspace for color transform + filtering */
    /* For RGB: apply YCoCg-R, expanding to 3 int16 channels */
    /* Then clamp to uint8 for filtering */
    int use_ycocg = (channels == 3);

    uint8_t *transformed = NULL;
    int trans_channels = channels;

    if (use_ycocg) {
        /* YCoCg-R produces signed values. We store as uint8 with offset 128 for Co,Cg */
        transformed = (uint8_t *)malloc(width * height * 3);
        if (!transformed) return -1;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int idx = (i * width + j) * 3;
                uint8_t r = pixels[idx], g = pixels[idx+1], b = pixels[idx+2];
                int16_t y, co, cg;
                rgb_to_ycocg(r, g, b, &y, &co, &cg);
                transformed[idx]   = (uint8_t)(y & 0xFF);
                transformed[idx+1] = (uint8_t)((co + 128) & 0xFF);  /* offset to unsigned */
                transformed[idx+2] = (uint8_t)((cg + 128) & 0xFF);
            }
        }
    } else {
        transformed = (uint8_t *)malloc(width * height * channels);
        if (!transformed) return -1;
        memcpy(transformed, pixels, width * height * channels);
    }

    /* Per-row adaptive filtering */
    size_t filtered_size = (size_t)height * (1 + row_bytes);  /* 1 byte filter type + row */
    uint8_t *filtered = (uint8_t *)malloc(filtered_size);
    if (!filtered) { free(transformed); return -1; }

    uint8_t *above = (uint8_t *)calloc(row_bytes, 1);
    uint8_t *temp = (uint8_t *)malloc(row_bytes);

    for (int i = 0; i < height; i++) {
        const uint8_t *row = transformed + i * row_bytes;
        uint8_t filt = best_filter(row, above, row_bytes, bpp);

        filtered[i * (1 + row_bytes)] = filt;
        apply_filter(filt, row, above, filtered + i * (1 + row_bytes) + 1, row_bytes, bpp);

        memcpy(above, row, row_bytes);
    }

    free(temp);
    free(above);
    free(transformed);

    /* Compress with zlib */
    uLongf compressed_len = compressBound(filtered_size);
    *out = (uint8_t *)malloc(sizeof(tc_header_t) + compressed_len);
    if (!*out) { free(filtered); return -1; }

    /* Write header */
    tc_header_t hdr;
    hdr.magic = TC_MAGIC;
    hdr.version = TC_VERSION;
    hdr.color_transform = use_ycocg ? 1 : 0;
    hdr.width = (uint16_t)width;
    hdr.height = (uint16_t)height;
    hdr.channels = (uint8_t)channels;
    hdr.bit_depth = 8;
    memcpy(*out, &hdr, sizeof(hdr));

    /* Compress */
    int ret = compress2(*out + sizeof(hdr), &compressed_len,
                         filtered, filtered_size, 9);
    free(filtered);

    if (ret != Z_OK) { free(*out); *out = NULL; return -1; }

    *out_len = sizeof(hdr) + compressed_len;
    return 0;
}

/* ================================================================
 * DECOMPRESS
 * ================================================================ */

int tc_decompress(const uint8_t *compressed, size_t comp_len,
                   uint8_t **pixels, int *width, int *height, int *channels) {
    if (comp_len < sizeof(tc_header_t)) return -1;

    tc_header_t hdr;
    memcpy(&hdr, compressed, sizeof(hdr));
    if (hdr.magic != TC_MAGIC || hdr.version != TC_VERSION) return -1;

    *width = hdr.width;
    *height = hdr.height;
    *channels = hdr.channels;
    int row_bytes = hdr.width * hdr.channels;
    int bpp = hdr.channels;

    /* Decompress */
    size_t filtered_size = (size_t)hdr.height * (1 + row_bytes);
    uint8_t *filtered = (uint8_t *)malloc(filtered_size);
    if (!filtered) return -1;

    uLongf filt_len = filtered_size;
    int ret = uncompress(filtered, &filt_len,
                          compressed + sizeof(hdr), comp_len - sizeof(hdr));
    if (ret != Z_OK) { free(filtered); return -1; }

    /* Undo filtering */
    uint8_t *transformed = (uint8_t *)malloc((size_t)hdr.width * hdr.height * hdr.channels);
    if (!transformed) { free(filtered); return -1; }

    uint8_t *above = (uint8_t *)calloc(row_bytes, 1);

    for (int i = 0; i < hdr.height; i++) {
        uint8_t filt = filtered[i * (1 + row_bytes)];
        uint8_t *row_out = transformed + i * row_bytes;
        undo_filter(filt, filtered + i * (1 + row_bytes) + 1, above, row_out, row_bytes, bpp);
        memcpy(above, row_out, row_bytes);
    }

    free(above);
    free(filtered);

    /* Undo YCoCg-R if applicable */
    *pixels = (uint8_t *)malloc((size_t)hdr.width * hdr.height * hdr.channels);
    if (!*pixels) { free(transformed); return -1; }

    if (hdr.color_transform == 1 && hdr.channels == 3) {
        for (int i = 0; i < hdr.height; i++) {
            for (int j = 0; j < hdr.width; j++) {
                int idx = (i * hdr.width + j) * 3;
                int16_t y  = (int16_t)transformed[idx];
                int16_t co = (int16_t)transformed[idx+1] - 128;
                int16_t cg = (int16_t)transformed[idx+2] - 128;
                uint8_t r, g, b;
                ycocg_to_rgb(y, co, cg, &r, &g, &b);
                (*pixels)[idx] = r;
                (*pixels)[idx+1] = g;
                (*pixels)[idx+2] = b;
            }
        }
        free(transformed);
    } else {
        free(*pixels);
        *pixels = transformed;
    }

    return 0;
}

/* ================================================================
 * PPM I/O (simple image format for testing)
 * ================================================================ */

static uint8_t *read_ppm(const char *path, int *w, int *h, int *ch) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    char magic[3]; int maxval;
    if (fscanf(f, "%2s", magic) != 1) { fclose(f); return NULL; }

    if (strcmp(magic, "P6") == 0) *ch = 3;
    else if (strcmp(magic, "P5") == 0) *ch = 1;
    else { fclose(f); return NULL; }

    if (fscanf(f, "%d %d %d", w, h, &maxval) != 3) { fclose(f); return NULL; }
    fgetc(f); /* consume newline */

    size_t sz = (size_t)(*w) * (*h) * (*ch);
    uint8_t *data = (uint8_t *)malloc(sz);
    if (!data || fread(data, 1, sz, f) != sz) { free(data); fclose(f); return NULL; }
    fclose(f);
    return data;
}

static int write_ppm(const char *path, const uint8_t *data, int w, int h, int ch) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    fprintf(f, "%s\n%d %d\n255\n", ch == 3 ? "P6" : "P5", w, h);
    fwrite(data, 1, (size_t)w * h * ch, f);
    fclose(f);
    return 0;
}

/* ================================================================
 * MAIN
 * ================================================================ */

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "TC v%d — Tournament Compression Library\n", TC_VERSION);
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s compress   input.ppm output.tc\n", argv[0]);
        fprintf(stderr, "  %s decompress output.tc recovered.ppm\n", argv[0]);
        fprintf(stderr, "  %s benchmark  input.ppm\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "compress") == 0 && argc >= 4) {
        int w, h, ch;
        uint8_t *pixels = read_ppm(argv[2], &w, &h, &ch);
        if (!pixels) { fprintf(stderr, "Failed to read %s\n", argv[2]); return 1; }

        uint8_t *compressed;
        size_t comp_len;
        if (tc_compress(pixels, w, h, ch, &compressed, &comp_len) != 0) {
            fprintf(stderr, "Compression failed\n"); free(pixels); return 1;
        }

        FILE *f = fopen(argv[3], "wb");
        fwrite(compressed, 1, comp_len, f);
        fclose(f);

        size_t raw = (size_t)w * h * ch;
        printf("TC: %zu -> %zu bytes (%.2fx)\n", raw, comp_len, (double)raw / comp_len);
        printf("  YCoCg-R: %s, Adaptive filters, zlib-9\n", ch == 3 ? "yes" : "no");

        free(pixels);
        free(compressed);

    } else if (strcmp(argv[1], "decompress") == 0 && argc >= 4) {
        FILE *f = fopen(argv[2], "rb");
        if (!f) { fprintf(stderr, "Failed to read %s\n", argv[2]); return 1; }
        fseek(f, 0, SEEK_END); size_t sz = ftell(f); fseek(f, 0, SEEK_SET);
        uint8_t *compressed = (uint8_t *)malloc(sz);
        fread(compressed, 1, sz, f); fclose(f);

        int w, h, ch;
        uint8_t *pixels;
        if (tc_decompress(compressed, sz, &pixels, &w, &h, &ch) != 0) {
            fprintf(stderr, "Decompression failed\n"); free(compressed); return 1;
        }

        write_ppm(argv[3], pixels, w, h, ch);
        printf("TC: decompressed %dx%d (%d channels)\n", w, h, ch);

        free(compressed);
        free(pixels);

    } else if (strcmp(argv[1], "benchmark") == 0 && argc >= 3) {
        int w, h, ch;
        uint8_t *pixels = read_ppm(argv[2], &w, &h, &ch);
        if (!pixels) { fprintf(stderr, "Failed to read %s\n", argv[2]); return 1; }

        size_t raw = (size_t)w * h * ch;

        /* TC compress */
        uint8_t *compressed;
        size_t comp_len;
        tc_compress(pixels, w, h, ch, &compressed, &comp_len);

        /* Verify roundtrip */
        int w2, h2, ch2;
        uint8_t *recovered;
        tc_decompress(compressed, comp_len, &recovered, &w2, &h2, &ch2);

        int match = (w == w2 && h == h2 && ch == ch2 &&
                     memcmp(pixels, recovered, raw) == 0);

        /* Plain zlib (no filtering, no decorrelation) for comparison */
        uLongf plain_len = compressBound(raw);
        uint8_t *plain = (uint8_t *)malloc(plain_len);
        compress2(plain, &plain_len, pixels, raw, 9);

        printf("TC benchmark: %s (%dx%d, %d ch)\n", argv[2], w, h, ch);
        printf("  Raw:      %10zu bytes\n", raw);
        printf("  zlib-9:   %10lu bytes (%.2fx)\n", (unsigned long)plain_len, (double)raw / plain_len);
        printf("  TC:       %10zu bytes (%.2fx)\n", comp_len, (double)raw / comp_len);
        printf("  TC/zlib:  %.4fx\n", (double)comp_len / plain_len);
        printf("  Roundtrip: %s\n", match ? "PASS" : "FAIL");

        free(plain);
        free(recovered);
        free(compressed);
        free(pixels);
    }

    return 0;
}
