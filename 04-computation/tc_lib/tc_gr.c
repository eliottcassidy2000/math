/*
 * tc_gr.c — Tournament Compression with Golomb-Rice (no zlib dependency)
 * opus-2026-03-26-S339
 *
 * Based on kind-pasteur's tic_gr2.c architecture.
 * Pipeline: YCoCg-R → MED prediction → 64-context adaptive Golomb-Rice
 *
 * Compile: cc -O2 -o tc_gr tc_gr.c -lm
 * Usage:   ./tc_gr compress   input.ppm output.tcg
 *          ./tc_gr decompress output.tcg recovered.ppm
 *          ./tc_gr benchmark  input.ppm
 *
 * NO EXTERNAL DEPENDENCIES beyond libc and libm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "tc.h"  /* for YCoCg-R and Paeth */

/* ================================================================
 * BITSTREAM
 * ================================================================ */

typedef struct { uint8_t *data; size_t cap, pos; int bit; } bitstream_t;

static void bs_init(bitstream_t *b, uint8_t *d, size_t cap) {
    b->data = d; b->cap = cap; b->pos = 0; b->bit = 7;
    memset(d, 0, cap);
}

static void bs_write_bit(bitstream_t *b, int v) {
    if (b->pos >= b->cap) return;
    if (v) b->data[b->pos] |= (1 << b->bit);
    if (--b->bit < 0) { b->bit = 7; b->pos++; }
}

static void bs_write_bits(bitstream_t *b, uint32_t v, int n) {
    for (int i = n - 1; i >= 0; i--) bs_write_bit(b, (v >> i) & 1);
}

static size_t bs_size(bitstream_t *b) {
    return b->pos + (b->bit < 7 ? 1 : 0);
}

static int bs_read_bit(bitstream_t *b) {
    if (b->pos >= b->cap) return 0;
    int v = (b->data[b->pos] >> b->bit) & 1;
    if (--b->bit < 0) { b->bit = 7; b->pos++; }
    return v;
}

static uint32_t bs_read_bits(bitstream_t *b, int n) {
    uint32_t v = 0;
    for (int i = 0; i < n; i++) v = (v << 1) | bs_read_bit(b);
    return v;
}

/* ================================================================
 * GOLOMB-RICE CODING
 * ================================================================ */

static int zigzag_enc(int v) { return v >= 0 ? 2 * v : 2 * (-v) - 1; }
static int zigzag_dec(int v) { return (v & 1) ? -((v + 1) / 2) : v / 2; }

static void gr_encode(bitstream_t *b, int val, int k) {
    int v = zigzag_enc(val);
    int q = v >> k, r = v & ((1 << k) - 1);
    if (q > 30) {
        for (int i = 0; i < 31; i++) bs_write_bit(b, 1);
        bs_write_bit(b, 0);
        bs_write_bits(b, v, 16);
    } else {
        for (int i = 0; i < q; i++) bs_write_bit(b, 1);
        bs_write_bit(b, 0);
        bs_write_bits(b, r, k);
    }
}

static int gr_decode(bitstream_t *b, int k) {
    int q = 0;
    while (bs_read_bit(b) == 1 && q < 31) q++;
    int v;
    if (q >= 31) v = bs_read_bits(b, 16);
    else v = (q << k) | bs_read_bits(b, k);
    return zigzag_dec(v);
}

/* ================================================================
 * MED PREDICTOR + CONTEXT
 * ================================================================ */

static inline int med_pred(int a, int b, int c) {
    int mn = a < b ? a : b, mx = a > b ? a : b;
    if (c >= mx) return mn;
    if (c <= mn) return mx;
    return a + b - c;
}

#define NCTX 64

static int get_context(int a, int b, int c, int d) {
    int g = abs(a - c) + abs(b - c) + abs(d - b);
    if (g < 3) return 0;
    if (g < 8) return 1;
    if (g < 15) return 2;
    if (g < 25) return 3;
    if (g < 40) return 4 + (g - 25) / 5;
    if (g < 100) return 7 + (g - 40) / 10;
    int r = 13 + (g - 100) / 20;
    return r < NCTX ? r : NCTX - 1;
}

static int k_from_avg(double avg) {
    if (avg < 0.5) return 0;
    int k = (int)(log2(avg / 0.693) + 0.5);
    return k < 0 ? 0 : (k > 12 ? 12 : k);
}

/* ================================================================
 * PLANE ENCODE / DECODE
 * ================================================================ */

static size_t encode_plane(const uint8_t *plane, int h, int w,
                            uint8_t *out, size_t cap) {
    double ctx_avg[NCTX];
    for (int i = 0; i < NCTX; i++) ctx_avg[i] = 4.0;

    bitstream_t bs;
    bs_init(&bs, out, cap);

    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int cv = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            int d = (r > 0 && c + 1 < w) ? plane[idx - w + 1] : b;

            int pred = med_pred(a, b, cv);
            int res = (int)plane[idx] - pred;
            if (res > 128) res -= 256;
            if (res < -128) res += 256;

            int cx = get_context(a, b, cv, d);
            int k = k_from_avg(ctx_avg[cx]);
            gr_encode(&bs, res, k);
            ctx_avg[cx] = 0.9 * ctx_avg[cx] + 0.1 * abs(res);
        }
    }
    return bs_size(&bs);
}

static void decode_plane(const uint8_t *data, size_t len,
                          uint8_t *plane, int h, int w) {
    double ctx_avg[NCTX];
    for (int i = 0; i < NCTX; i++) ctx_avg[i] = 4.0;

    bitstream_t bs;
    bs.data = (uint8_t *)data;
    bs.cap = len;
    bs.pos = 0;
    bs.bit = 7;

    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int cv = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            int d = (r > 0 && c + 1 < w) ? plane[idx - w + 1] : b;

            int cx = get_context(a, b, cv, d);
            int k = k_from_avg(ctx_avg[cx]);
            int res = gr_decode(&bs, k);
            plane[idx] = (uint8_t)((med_pred(a, b, cv) + res) & 0xFF);
            ctx_avg[cx] = 0.9 * ctx_avg[cx] + 0.1 * abs(res);
        }
    }
}

/* ================================================================
 * RGB COMPRESS / DECOMPRESS
 * ================================================================ */

#define TCG_MAGIC 0x47435421  /* "!TCG" */

int tcg_compress(const uint8_t *pixels, int w, int h, int ch,
                  uint8_t **out, size_t *out_len) {
    int n = w * h;

    /* Color transform: YCoCg-R for RGB, identity for gray */
    uint8_t *p1 = malloc(n), *p2 = malloc(n), *p3 = malloc(n);
    if (!p1 || !p2 || !p3) return -1;

    if (ch == 3) {
        for (int i = 0; i < n; i++) {
            p1[i] = pixels[3*i + 1];                         /* G */
            p2[i] = (pixels[3*i] - pixels[3*i + 1]) & 0xFF;  /* R-G */
            p3[i] = (pixels[3*i + 2] - pixels[3*i + 1]) & 0xFF; /* B-G */
        }
    } else {
        memcpy(p1, pixels, n);
    }

    /* Allocate output (worst case: 2× raw) */
    size_t max_out = 16 + 3 * n * 2;
    *out = malloc(max_out);
    if (!*out) { free(p1); free(p2); free(p3); return -1; }

    /* Header */
    uint8_t *ptr = *out;
    *(uint32_t *)ptr = TCG_MAGIC; ptr += 4;
    *(uint16_t *)ptr = w; ptr += 2;
    *(uint16_t *)ptr = h; ptr += 2;
    *ptr++ = ch;
    *ptr++ = 0; /* reserved */

    /* Encode each plane */
    size_t s1 = encode_plane(p1, h, w, ptr + 4, max_out - (ptr - *out) - 4);
    *(uint32_t *)ptr = (uint32_t)s1; ptr += 4;
    ptr += s1;

    if (ch == 3) {
        size_t s2 = encode_plane(p2, h, w, ptr + 4, max_out - (ptr - *out) - 4);
        *(uint32_t *)ptr = (uint32_t)s2; ptr += 4;
        ptr += s2;

        size_t s3 = encode_plane(p3, h, w, ptr + 4, max_out - (ptr - *out) - 4);
        *(uint32_t *)ptr = (uint32_t)s3; ptr += 4;
        ptr += s3;
    }

    *out_len = ptr - *out;
    free(p1); free(p2); free(p3);
    return 0;
}

int tcg_decompress(const uint8_t *comp, size_t comp_len,
                    uint8_t **pixels, int *w, int *h, int *ch) {
    if (comp_len < 10) return -1;
    if (*(uint32_t *)comp != TCG_MAGIC) return -1;

    *w = *(uint16_t *)(comp + 4);
    *h = *(uint16_t *)(comp + 6);
    *ch = comp[8];
    int n = (*w) * (*h);

    uint8_t *p1 = calloc(n, 1), *p2 = calloc(n, 1), *p3 = calloc(n, 1);
    if (!p1 || !p2 || !p3) return -1;

    const uint8_t *ptr = comp + 10;
    uint32_t s1 = *(uint32_t *)ptr; ptr += 4;
    decode_plane(ptr, s1, p1, *h, *w); ptr += s1;

    if (*ch == 3) {
        uint32_t s2 = *(uint32_t *)ptr; ptr += 4;
        decode_plane(ptr, s2, p2, *h, *w); ptr += s2;

        uint32_t s3 = *(uint32_t *)ptr; ptr += 4;
        decode_plane(ptr, s3, p3, *h, *w); ptr += s3;
    }

    /* Undo color transform */
    *pixels = malloc(n * (*ch));
    if (!*pixels) { free(p1); free(p2); free(p3); return -1; }

    if (*ch == 3) {
        for (int i = 0; i < n; i++) {
            (*pixels)[3*i + 1] = p1[i];
            (*pixels)[3*i] = (p2[i] + p1[i]) & 0xFF;
            (*pixels)[3*i + 2] = (p3[i] + p1[i]) & 0xFF;
        }
    } else {
        memcpy(*pixels, p1, n);
    }

    free(p1); free(p2); free(p3);
    return 0;
}

/* ================================================================
 * PPM I/O + MAIN (same as tc.c)
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
    fgetc(f);
    size_t sz = (size_t)(*w) * (*h) * (*ch);
    uint8_t *data = malloc(sz);
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

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "TC-GR — Tournament Compression with Golomb-Rice\n");
        fprintf(stderr, "No external dependencies. Pure C + libm.\n\n");
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s compress   input.ppm output.tcg\n", argv[0]);
        fprintf(stderr, "  %s decompress output.tcg recovered.ppm\n", argv[0]);
        fprintf(stderr, "  %s benchmark  input.ppm\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "compress") == 0 && argc >= 4) {
        int w, h, ch;
        uint8_t *pixels = read_ppm(argv[2], &w, &h, &ch);
        if (!pixels) { fprintf(stderr, "Failed to read %s\n", argv[2]); return 1; }

        uint8_t *compressed; size_t comp_len;
        clock_t t0 = clock();
        tcg_compress(pixels, w, h, ch, &compressed, &comp_len);
        double ms = (double)(clock() - t0) / CLOCKS_PER_SEC * 1000;

        FILE *f = fopen(argv[3], "wb");
        fwrite(compressed, 1, comp_len, f);
        fclose(f);

        size_t raw = (size_t)w * h * ch;
        printf("TC-GR: %zu -> %zu bytes (%.2fx) in %.1fms\n", raw, comp_len, (double)raw / comp_len, ms);

        free(pixels); free(compressed);

    } else if (strcmp(argv[1], "decompress") == 0 && argc >= 4) {
        FILE *f = fopen(argv[2], "rb");
        fseek(f, 0, SEEK_END); size_t sz = ftell(f); fseek(f, 0, SEEK_SET);
        uint8_t *comp = malloc(sz); fread(comp, 1, sz, f); fclose(f);

        int w, h, ch; uint8_t *pixels;
        tcg_decompress(comp, sz, &pixels, &w, &h, &ch);
        write_ppm(argv[3], pixels, w, h, ch);
        printf("TC-GR: decompressed %dx%d (%d ch)\n", w, h, ch);
        free(comp); free(pixels);

    } else if (strcmp(argv[1], "benchmark") == 0 && argc >= 3) {
        int w, h, ch;
        uint8_t *pixels = read_ppm(argv[2], &w, &h, &ch);
        if (!pixels) { fprintf(stderr, "Failed to read %s\n", argv[2]); return 1; }
        size_t raw = (size_t)w * h * ch;

        uint8_t *compressed; size_t comp_len;
        clock_t t0 = clock();
        tcg_compress(pixels, w, h, ch, &compressed, &comp_len);
        double enc_ms = (double)(clock() - t0) / CLOCKS_PER_SEC * 1000;

        /* Verify roundtrip */
        int w2, h2, ch2; uint8_t *recovered;
        tcg_decompress(compressed, comp_len, &recovered, &w2, &h2, &ch2);
        int match = (w == w2 && h == h2 && ch == ch2 && memcmp(pixels, recovered, raw) == 0);

        printf("TC-GR benchmark: %s (%dx%d, %d ch)\n", argv[2], w, h, ch);
        printf("  Raw:    %10zu bytes\n", raw);
        printf("  TC-GR:  %10zu bytes (%.2fx) in %.1fms\n", comp_len, (double)raw / comp_len, enc_ms);
        printf("  Roundtrip: %s\n", match ? "PASS" : "FAIL");

        free(pixels); free(compressed); free(recovered);
    }

    return 0;
}
