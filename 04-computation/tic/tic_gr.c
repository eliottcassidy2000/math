/*
 * tic_gr.c — TIC with Golomb-Rice entropy coding
 *
 * Replaces zlib with context-adaptive Golomb-Rice.
 * This is the JPEG-LS-like approach: MED prediction + adaptive GR.
 *
 * Pipeline: G-RG-BG decorrelation → MED prediction → context-adaptive GR
 *
 * Context: quantized gradient magnitude (64 contexts).
 * Each context maintains a running average of |residual|.
 * The GR parameter k is derived from the running average.
 *
 * Run-length mode: when residual = 0, encode the RUN LENGTH instead
 * of individual zeros. This dramatically helps flat regions.
 *
 * Compile: gcc -O2 -o tic_gr tic_gr.c -lm
 * (No zlib dependency!)
 *
 * kind-pasteur-2026-03-26-S32
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

/* ================================================================
 * BIT WRITER / READER
 * ================================================================ */

typedef struct {
    uint8_t *data;
    size_t capacity;
    size_t byte_pos;
    int bit_pos;  /* 7 = MSB, 0 = LSB */
} bitstream_t;

static void bs_init(bitstream_t *bs, uint8_t *buf, size_t cap) {
    bs->data = buf; bs->capacity = cap; bs->byte_pos = 0; bs->bit_pos = 7;
}

static void bs_write_bit(bitstream_t *bs, int bit) {
    if (bs->byte_pos >= bs->capacity) return;
    if (bit) bs->data[bs->byte_pos] |= (1 << bs->bit_pos);
    if (--bs->bit_pos < 0) { bs->bit_pos = 7; bs->byte_pos++; if (bs->byte_pos < bs->capacity) bs->data[bs->byte_pos] = 0; }
}

static void bs_write_bits(bitstream_t *bs, uint32_t val, int n) {
    for (int i = n - 1; i >= 0; i--)
        bs_write_bit(bs, (val >> i) & 1);
}

static size_t bs_bytes_written(bitstream_t *bs) {
    return bs->byte_pos + (bs->bit_pos < 7 ? 1 : 0);
}

static int bs_read_bit(bitstream_t *bs) {
    if (bs->byte_pos >= bs->capacity) return 0;
    int bit = (bs->data[bs->byte_pos] >> bs->bit_pos) & 1;
    if (--bs->bit_pos < 0) { bs->bit_pos = 7; bs->byte_pos++; }
    return bit;
}

static uint32_t bs_read_bits(bitstream_t *bs, int n) {
    uint32_t val = 0;
    for (int i = 0; i < n; i++)
        val = (val << 1) | bs_read_bit(bs);
    return val;
}

/* ================================================================
 * GOLOMB-RICE CODING
 * ================================================================ */

static int zigzag_enc(int v) { return v >= 0 ? 2*v : 2*(-v)-1; }
static int zigzag_dec(int v) { return (v & 1) ? -((v+1)/2) : v/2; }

static void gr_encode(bitstream_t *bs, int value, int k) {
    int v = zigzag_enc(value);
    int q = v >> k;
    int r = v & ((1 << k) - 1);
    /* Unary: q ones + 1 zero. Cap q to prevent runaway. */
    int capped_q = q < 32 ? q : 31;
    for (int i = 0; i < capped_q; i++) bs_write_bit(bs, 1);
    bs_write_bit(bs, 0);
    if (q >= 32) bs_write_bits(bs, v, 16); /* escape: raw 16-bit */
    else bs_write_bits(bs, r, k);
}

static int gr_decode(bitstream_t *bs, int k) {
    int q = 0;
    while (bs_read_bit(bs) == 1 && q < 32) q++;
    int v;
    if (q >= 32) { v = bs_read_bits(bs, 16); }
    else { int r = bs_read_bits(bs, k); v = (q << k) | r; }
    return zigzag_dec(v);
}

/* ================================================================
 * MED PREDICTOR
 * ================================================================ */

static inline int med_pred(int a, int b, int c) {
    int mn = a < b ? a : b, mx = a > b ? a : b;
    if (c >= mx) return mn;
    if (c <= mn) return mx;
    return a + b - c;
}

/* ================================================================
 * CONTEXT MODEL: gradient-based, 64 contexts
 * ================================================================ */

#define N_CTX 64
#define RUN_THRESHOLD 0  /* Enter run mode when residual == 0 */

static int get_context(int a, int b, int c, int d) {
    /* a=left, b=above, c=above-left, d=above-right */
    int gh = abs(a - c) + abs(d - b);  /* horizontal gradient */
    int gv = abs(b - c) + abs(a - (a > 0 ? a : 0));  /* simplified vertical */
    int g = gh + gv;
    /* Quantize to 64 levels */
    if (g < 2) return 0;
    if (g < 5) return 1;
    if (g < 10) return 2;
    if (g < 20) return 3;
    if (g < 40) return 4 + (g - 20) / 5;  /* 4-7 */
    if (g < 100) return 8 + (g - 40) / 10; /* 8-13 */
    return 14 + (g - 100) / 20;  /* 14-63, capped */
}

static int k_from_average(double avg) {
    if (avg < 0.5) return 0;
    int k = (int)(log2(avg / 0.6931) + 0.5);  /* ln(2) ≈ 0.6931 */
    return k < 0 ? 0 : (k > 12 ? 12 : k);
}

/* ================================================================
 * ENCODE ONE PLANE: MED + context-adaptive GR + run-length
 * ================================================================ */

static size_t encode_plane_gr(const uint8_t *plane, int h, int w,
                                uint8_t *out, size_t out_cap) {
    double ctx_avg[N_CTX];
    for (int i = 0; i < N_CTX; i++) ctx_avg[i] = 4.0;

    bitstream_t bs;
    memset(out, 0, out_cap);
    bs_init(&bs, out, out_cap);

    int run_count = 0;
    int in_run = 0;

    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx-1] : 0;
            int b = (r > 0) ? plane[idx-w] : 0;
            int cv = (r > 0 && c > 0) ? plane[idx-w-1] : 0;
            int d = (r > 0 && c+1 < w) ? plane[idx-w+1] : b;

            int pred = med_pred(a, b, cv);
            int residual = (int)plane[idx] - pred;
            if (residual > 128) residual -= 256;
            if (residual < -128) residual += 256;

            /* Run-length mode: encode consecutive zeros as a run */
            if (residual == 0) {
                run_count++;
                continue;
            }

            /* Flush run if we had one */
            if (run_count > 0) {
                /* Encode run length using GR with k=2 */
                gr_encode(&bs, run_count, 2);
                run_count = 0;
                /* Signal end of run with a marker bit */
                bs_write_bit(&bs, 1);
            } else {
                bs_write_bit(&bs, 0); /* no run before this residual */
            }

            /* Context-adaptive GR */
            int ctx = get_context(a, b, cv, d);
            if (ctx >= N_CTX) ctx = N_CTX - 1;
            int k = k_from_average(ctx_avg[ctx]);
            gr_encode(&bs, residual, k);

            /* Update context */
            ctx_avg[ctx] = 0.9 * ctx_avg[ctx] + 0.1 * abs(residual);
        }
    }

    /* Final run flush */
    if (run_count > 0) {
        gr_encode(&bs, run_count, 2);
        run_count = 0;
        bs_write_bit(&bs, 1);
    }

    return bs_bytes_written(&bs);
}

/* ================================================================
 * DECODE ONE PLANE
 * ================================================================ */

static int decode_plane_gr(const uint8_t *data, size_t data_len,
                             uint8_t *plane, int h, int w) {
    double ctx_avg[N_CTX];
    for (int i = 0; i < N_CTX; i++) ctx_avg[i] = 4.0;

    bitstream_t bs;
    bs_init((bitstream_t*)&bs, (uint8_t*)data, data_len);
    bs.byte_pos = 0; bs.bit_pos = 7;

    int pos = 0;
    while (pos < h * w) {
        /* Check for run */
        int run_bit = bs_read_bit((bitstream_t*)&bs);
        if (run_bit) {
            /* Decode run */
            int run_len = gr_decode((bitstream_t*)&bs, 2);
            /* Fill with MED predictions (residual = 0) */
            for (int i = 0; i < run_len && pos < h * w; i++, pos++) {
                int r = pos / w, c = pos % w;
                int a = (c > 0) ? plane[pos-1] : 0;
                int b = (r > 0) ? plane[pos-w] : 0;
                int cv = (r > 0 && c > 0) ? plane[pos-w-1] : 0;
                plane[pos] = (uint8_t)(med_pred(a, b, cv) & 0xFF);
            }
            continue;
        }

        /* Non-zero residual */
        int r = pos / w, c = pos % w;
        int a = (c > 0) ? plane[pos-1] : 0;
        int b = (r > 0) ? plane[pos-w] : 0;
        int cv = (r > 0 && c > 0) ? plane[pos-w-1] : 0;
        int d = (r > 0 && c+1 < w) ? plane[pos-w+1] : b;

        int ctx = get_context(a, b, cv, d);
        if (ctx >= N_CTX) ctx = N_CTX - 1;
        int k = k_from_average(ctx_avg[ctx]);
        int residual = gr_decode((bitstream_t*)&bs, k);

        int pred = med_pred(a, b, cv);
        plane[pos] = (uint8_t)((pred + residual) & 0xFF);

        ctx_avg[ctx] = 0.9 * ctx_avg[ctx] + 0.1 * abs(residual);
        pos++;
    }

    return 0;
}

/* ================================================================
 * COLOR TRANSFORM
 * ================================================================ */

static void rgb_to_grd(const uint8_t *rgb, uint8_t *g, uint8_t *rg,
                         uint8_t *bg, int n) {
    for (int i = 0; i < n; i++) {
        g[i] = rgb[3*i+1];
        rg[i] = (rgb[3*i] - rgb[3*i+1]) & 0xFF;
        bg[i] = (rgb[3*i+2] - rgb[3*i+1]) & 0xFF;
    }
}

static void grd_to_rgb(const uint8_t *g, const uint8_t *rg,
                         const uint8_t *bg, uint8_t *rgb, int n) {
    for (int i = 0; i < n; i++) {
        rgb[3*i+1] = g[i];
        rgb[3*i] = (rg[i] + g[i]) & 0xFF;
        rgb[3*i+2] = (bg[i] + g[i]) & 0xFF;
    }
}

/* ================================================================
 * FULL ENCODE / DECODE
 * ================================================================ */

static const uint8_t MAGIC[4] = {'T','G','R','1'};

long tic_gr_encode(const uint8_t *rgb, int w, int h,
                    uint8_t *out, size_t out_cap) {
    int n = w * h;
    uint8_t *g = malloc(n), *rg = malloc(n), *bg = malloc(n);
    size_t buf_sz = n + n/2;  /* GR can exceed raw for pathological input */
    uint8_t *buf = malloc(buf_sz);

    rgb_to_grd(rgb, g, rg, bg, n);

    size_t pos = 0;
    memcpy(out + pos, MAGIC, 4); pos += 4;
    out[pos++] = w & 0xFF; out[pos++] = (w >> 8) & 0xFF;
    out[pos++] = h & 0xFF; out[pos++] = (h >> 8) & 0xFF;
    out[pos++] = 0; /* flags */

    uint8_t *planes[3] = {g, rg, bg};
    for (int p = 0; p < 3; p++) {
        size_t plane_sz = encode_plane_gr(planes[p], h, w, buf, buf_sz);
        /* Write plane length + data */
        out[pos++] = plane_sz & 0xFF;
        out[pos++] = (plane_sz >> 8) & 0xFF;
        out[pos++] = (plane_sz >> 16) & 0xFF;
        out[pos++] = (plane_sz >> 24) & 0xFF;
        memcpy(out + pos, buf, plane_sz);
        pos += plane_sz;
    }

    free(g); free(rg); free(bg); free(buf);
    return (long)pos;
}

int tic_gr_decode(const uint8_t *data, size_t data_len,
                   int w, int h, uint8_t *rgb) {
    int n = w * h;
    uint8_t *g = malloc(n), *rg = malloc(n), *bg = malloc(n);

    size_t pos = 9;
    uint8_t *planes[3] = {g, rg, bg};
    for (int p = 0; p < 3; p++) {
        uint32_t plen = data[pos] | (data[pos+1]<<8) | (data[pos+2]<<16) | (data[pos+3]<<24);
        pos += 4;
        decode_plane_gr(data + pos, plen, planes[p], h, w);
        pos += plen;
    }

    grd_to_rgb(g, rg, bg, rgb, n);
    free(g); free(rg); free(bg);
    return 0;
}

/* ================================================================
 * MAIN: benchmark
 * ================================================================ */

int main(int argc, char **argv) {
    if (argc < 5) {
        printf("Usage: %s bench input.raw W H [iters]\n", argv[0]);
        return 1;
    }

    int w = atoi(argv[3]), h = atoi(argv[4]);
    int iters = (argc >= 6) ? atoi(argv[5]) : 50;
    int n = w * h * 3;

    uint8_t *rgb = malloc(n);
    FILE *f = fopen(argv[2], "rb");
    if (!f) { printf("Cannot open %s\n", argv[2]); return 1; }
    if (fread(rgb, 1, n, f) != (size_t)n) { printf("Short read\n"); fclose(f); return 1; }
    fclose(f);

    size_t out_cap = n * 2;
    uint8_t *out = malloc(out_cap);

    /* Warmup */
    long sz = tic_gr_encode(rgb, w, h, out, out_cap);

    /* Encode benchmark */
    clock_t t0 = clock();
    for (int i = 0; i < iters; i++)
        sz = tic_gr_encode(rgb, w, h, out, out_cap);
    double enc_t = (double)(clock() - t0) / CLOCKS_PER_SEC / iters;

    /* Decode + verify */
    uint8_t *dec = malloc(n);
    t0 = clock();
    for (int i = 0; i < iters; i++)
        tic_gr_decode(out, sz, w, h, dec);
    double dec_t = (double)(clock() - t0) / CLOCKS_PER_SEC / iters;

    int ok = (memcmp(rgb, dec, n) == 0);

    printf("TIC-GR: %dx%d, %d iters\n", w, h, iters);
    printf("  Raw:        %d\n", n);
    printf("  Compressed: %ld (%.2f:1)\n", sz, (double)n / sz);
    printf("  Encode:     %.3f ms (%.1f Mpix/s)\n", enc_t*1000, (double)(w*h)/enc_t/1e6);
    printf("  Decode:     %.3f ms (%.1f Mpix/s)\n", dec_t*1000, (double)(w*h)/dec_t/1e6);
    printf("  Roundtrip:  %s\n", ok ? "PASS" : "FAIL");

    free(rgb); free(out); free(dec);
    return ok ? 0 : 1;
}
