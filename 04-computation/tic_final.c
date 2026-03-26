/*
 * TIC FINAL — The unified codec combining all innovations.
 *
 * From opus: YCoCg-R exact reversible color transform, per-row adaptive filtering
 * From kind-pasteur: MED predictor (6th filter), G-RG-BG transform option
 * Combined: 6-filter per-row adaptive + best color transform selection
 *
 * Compile: gcc -O3 -I<zlib-dir> -o tic_final tic_final.c <zlib-dir>/libz.a
 * Usage:   tic_final bench input.raw W H [iters]
 *          tic_final encode input.raw W H output.tic
 *
 * kind-pasteur-2026-03-25-S23
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

/* ================================================================
 * PREDICTORS
 * ================================================================ */

static inline int clamp8(int v) { return v < 0 ? 0 : (v > 255 ? 255 : v); }

static inline int paeth(int a, int b, int c) {
    int p = a + b - c;
    int pa = abs(p - a), pb = abs(p - b), pc = abs(p - c);
    if (pa <= pb && pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

static inline int med_pred(int a, int b, int c) {
    int mn = a < b ? a : b;
    int mx = a > b ? a : b;
    if (c >= mx) return mn;
    if (c <= mn) return mx;
    return a + b - c;
}

/* ================================================================
 * COLOR TRANSFORMS (both exact-reversible)
 * ================================================================ */

/* G-RG-BG: simple difference from green */
static void rgb_to_grd(const unsigned char *rgb, unsigned char *p0,
                        unsigned char *p1, unsigned char *p2, int n) {
    for (int i = 0; i < n; i++) {
        int r = rgb[3*i], g = rgb[3*i+1], b = rgb[3*i+2];
        p0[i] = g;
        p1[i] = (r - g) & 0xFF;
        p2[i] = (b - g) & 0xFF;
    }
}

/* YCoCg-R: exact reversible integer transform */
static void rgb_to_ycocg(const unsigned char *rgb, unsigned char *p0,
                          unsigned char *p1, unsigned char *p2, int n) {
    for (int i = 0; i < n; i++) {
        int r = rgb[3*i], g = rgb[3*i+1], b = rgb[3*i+2];
        int co = r - b;
        int t = b + (co >> 1);
        int cg = g - t;
        int y = t + (cg >> 1);
        p0[i] = y & 0xFF;
        p1[i] = co & 0xFF;
        p2[i] = cg & 0xFF;
    }
}

static void grd_to_rgb(const unsigned char *p0, const unsigned char *p1,
                        const unsigned char *p2, unsigned char *rgb, int n) {
    for (int i = 0; i < n; i++) {
        int g = p0[i], rg = (signed char)p1[i], bg = (signed char)p2[i];
        /* Use unsigned modular arithmetic */
        rgb[3*i+1] = p0[i];  /* G */
        rgb[3*i] = (p1[i] + p0[i]) & 0xFF;  /* R = (R-G) + G */
        rgb[3*i+2] = (p2[i] + p0[i]) & 0xFF;  /* B = (B-G) + G */
    }
}

static void ycocg_to_rgb(const unsigned char *p0, const unsigned char *p1,
                          const unsigned char *p2, unsigned char *rgb, int n) {
    for (int i = 0; i < n; i++) {
        int y = p0[i];
        int co = (signed char)p1[i];
        int cg = (signed char)p2[i];
        int t = y - (cg >> 1);
        int g = (t + cg) & 0xFF;
        int b = (t - (co >> 1)) & 0xFF;
        int r = (b + co) & 0xFF;
        rgb[3*i] = r; rgb[3*i+1] = g; rgb[3*i+2] = b;
    }
}

/* ================================================================
 * 6-FILTER PER-ROW ADAPTIVE PREDICTION
 * Filters 0-4 are PNG standard; Filter 5 is MED (our addition)
 * ================================================================ */

#define NFILTERS 6

static void apply_filter_row(int filt, const unsigned char *row,
                              const unsigned char *above,
                              unsigned char *out, int w) {
    for (int c = 0; c < w; c++) {
        int a = (c > 0) ? row[c-1] : 0;
        int b = above ? above[c] : 0;
        int d = (c > 0 && above) ? above[c-1] : 0;

        int pred;
        switch (filt) {
            case 0: pred = 0; break;          /* None */
            case 1: pred = a; break;           /* Sub */
            case 2: pred = b; break;           /* Up */
            case 3: pred = (a + b) >> 1; break; /* Avg */
            case 4: pred = paeth(a, b, d); break; /* Paeth */
            case 5: pred = med_pred(a, b, d); break; /* MED (NEW) */
            default: pred = 0;
        }
        out[c] = (row[c] - pred) & 0xFF;
    }
}

static void undo_filter_row(int filt, const unsigned char *filtered,
                              const unsigned char *above,
                              unsigned char *out, int w) {
    for (int c = 0; c < w; c++) {
        int a = (c > 0) ? out[c-1] : 0;
        int b = above ? above[c] : 0;
        int d = (c > 0 && above) ? above[c-1] : 0;

        int pred;
        switch (filt) {
            case 0: pred = 0; break;
            case 1: pred = a; break;
            case 2: pred = b; break;
            case 3: pred = (a + b) >> 1; break;
            case 4: pred = paeth(a, b, d); break;
            case 5: pred = med_pred(a, b, d); break;
            default: pred = 0;
        }
        out[c] = (filtered[c] + pred) & 0xFF;
    }
}

/* Pick best filter for a row using min-sum heuristic */
static int best_filter_row(const unsigned char *row, const unsigned char *above, int w) {
    unsigned char temp[w];
    int best_f = 0;
    long best_score = (long)w * 256;

    for (int f = 0; f < NFILTERS; f++) {
        apply_filter_row(f, row, above, temp, w);
        long score = 0;
        for (int c = 0; c < w; c++) {
            int v = temp[c];
            score += (v < 128) ? v : (256 - v);
        }
        if (score < best_score) {
            best_score = score;
            best_f = f;
        }
    }
    return best_f;
}

/* ================================================================
 * ENCODE / DECODE A SINGLE PLANE
 * ================================================================ */

/* Encode: per-row adaptive → interleaved [filter_byte, filtered_row...] → zlib */
static long encode_plane(const unsigned char *plane, int h, int w,
                          unsigned char *output, long outmax) {
    int row_len = 1 + w;  /* filter byte + w pixels */
    unsigned char *buf = malloc(h * row_len);
    const unsigned char *above = NULL;

    for (int r = 0; r < h; r++) {
        const unsigned char *row = plane + r * w;
        int f = best_filter_row(row, above, w);
        buf[r * row_len] = f;
        apply_filter_row(f, row, above, buf + r * row_len + 1, w);
        above = row;
    }

    uLongf dest_len = outmax;
    int ret = compress2(output, &dest_len, buf, h * row_len, 9);
    free(buf);
    return (ret == Z_OK) ? (long)dest_len : -1;
}

/* Decode: zlib → un-interleave → undo filter per row */
static int decode_plane(const unsigned char *input, long input_len,
                         unsigned char *plane, int h, int w) {
    int row_len = 1 + w;
    unsigned char *buf = malloc(h * row_len);

    uLongf dest_len = h * row_len;
    if (uncompress(buf, &dest_len, input, input_len) != Z_OK) {
        free(buf); return -1;
    }

    unsigned char *above = NULL;
    for (int r = 0; r < h; r++) {
        int f = buf[r * row_len];
        undo_filter_row(f, buf + r * row_len + 1, above, plane + r * w, w);
        above = plane + r * w;
    }

    free(buf);
    return 0;
}

/* ================================================================
 * FULL RGB ENCODE/DECODE
 * ================================================================ */

/* Try both color transforms, pick smaller */
long tic_encode(const unsigned char *rgb, int w, int h,
                unsigned char *output, long outmax) {
    int n = w * h;
    unsigned char *p0 = malloc(n), *p1 = malloc(n), *p2 = malloc(n);
    unsigned char *comp = malloc(n + 4096);

    /* Try G-RG-BG */
    rgb_to_grd(rgb, p0, p1, p2, n);
    long grd_total = 0;
    unsigned char *grd_buf = malloc(n * 3 + 65536);
    long pos = 0;

    for (int p = 0; p < 3; p++) {
        unsigned char *plane = (p == 0) ? p0 : (p == 1) ? p1 : p2;
        long clen = encode_plane(plane, h, w, comp, n + 4096);
        if (clen < 0) { free(p0); free(p1); free(p2); free(comp); free(grd_buf); return -1; }
        memcpy(grd_buf + pos, &clen, 4); pos += 4;
        memcpy(grd_buf + pos, comp, clen); pos += clen;
    }
    grd_total = pos;

    /* Try YCoCg-R */
    rgb_to_ycocg(rgb, p0, p1, p2, n);
    long yc_total = 0;
    unsigned char *yc_buf = malloc(n * 3 + 65536);
    pos = 0;

    for (int p = 0; p < 3; p++) {
        unsigned char *plane = (p == 0) ? p0 : (p == 1) ? p1 : p2;
        long clen = encode_plane(plane, h, w, comp, n + 4096);
        if (clen < 0) { free(p0); free(p1); free(p2); free(comp); free(grd_buf); free(yc_buf); return -1; }
        memcpy(yc_buf + pos, &clen, 4); pos += 4;
        memcpy(yc_buf + pos, comp, clen); pos += clen;
    }
    yc_total = pos;

    /* Header: magic(4) + w(2) + h(2) + ct(1) = 9 bytes */
    long out_pos = 0;
    memcpy(output, "TIC2", 4); out_pos += 4;
    output[out_pos++] = w & 0xFF; output[out_pos++] = (w >> 8) & 0xFF;
    output[out_pos++] = h & 0xFF; output[out_pos++] = (h >> 8) & 0xFF;

    if (grd_total <= yc_total) {
        output[out_pos++] = 0; /* G-RG-BG */
        memcpy(output + out_pos, grd_buf, grd_total); out_pos += grd_total;
    } else {
        output[out_pos++] = 1; /* YCoCg-R */
        memcpy(output + out_pos, yc_buf, yc_total); out_pos += yc_total;
    }

    free(p0); free(p1); free(p2); free(comp); free(grd_buf); free(yc_buf);
    return out_pos;
}

int tic_decode(const unsigned char *input, long input_len, int w, int h,
               unsigned char *rgb) {
    int n = w * h;
    unsigned char *p0 = malloc(n), *p1 = malloc(n), *p2 = malloc(n);

    long pos = 9; /* skip header */
    int ct = input[8];

    for (int p = 0; p < 3; p++) {
        long clen;
        memcpy(&clen, input + pos, 4); pos += 4;
        unsigned char *plane = (p == 0) ? p0 : (p == 1) ? p1 : p2;
        decode_plane(input + pos, clen, plane, h, w);
        pos += clen;
    }

    if (ct == 0) grd_to_rgb(p0, p1, p2, rgb, n);
    else ycocg_to_rgb(p0, p1, p2, rgb, n);

    free(p0); free(p1); free(p2);
    return 0;
}

/* ================================================================
 * MAIN: benchmark mode
 * ================================================================ */

int main(int argc, char **argv) {
    if (argc < 5) {
        printf("Usage: %s bench input.raw W H [iters]\n", argv[0]);
        printf("       %s encode input.raw W H output.tic\n", argv[0]);
        return 1;
    }

    int w = atoi(argv[3]), h = atoi(argv[4]);
    int n = w * h * 3;

    unsigned char *rgb = malloc(n);
    FILE *f = fopen(argv[2], "rb");
    if (!f) { printf("Cannot open %s\n", argv[2]); return 1; }
    fread(rgb, 1, n, f); fclose(f);

    unsigned char *out = malloc(n + 65536);

    if (strcmp(argv[1], "bench") == 0) {
        int iters = (argc >= 6) ? atoi(argv[5]) : 100;

        /* Warmup */
        long sz = tic_encode(rgb, w, h, out, n + 65536);

        /* Encode benchmark */
        clock_t t0 = clock();
        for (int i = 0; i < iters; i++)
            sz = tic_encode(rgb, w, h, out, n + 65536);
        double enc_t = (double)(clock() - t0) / CLOCKS_PER_SEC / iters;

        /* Decode benchmark */
        unsigned char *dec = malloc(n);
        t0 = clock();
        for (int i = 0; i < iters; i++)
            tic_decode(out, sz, w, h, dec);
        double dec_t = (double)(clock() - t0) / CLOCKS_PER_SEC / iters;

        int ok = (memcmp(rgb, dec, n) == 0);
        int ct = out[8]; /* color transform used */

        printf("TIC FINAL: %dx%d, %d iters, ct=%s\n", w, h, iters,
               ct == 0 ? "G-RG-BG" : "YCoCg-R");
        printf("  Raw:        %d\n", n);
        printf("  Compressed: %ld (%.1f:1)\n", sz, (double)n / sz);
        printf("  Encode:     %.3f ms (%.1f Mpix/s)\n", enc_t*1000, (double)(w*h)/enc_t/1e6);
        printf("  Decode:     %.3f ms (%.1f Mpix/s)\n", dec_t*1000, (double)(w*h)/dec_t/1e6);
        printf("  Roundtrip:  %s\n", ok ? "PASS" : "FAIL");

        free(dec);
    } else if (strcmp(argv[1], "encode") == 0 && argc >= 6) {
        long sz = tic_encode(rgb, w, h, out, n + 65536);
        FILE *fo = fopen(argv[5], "wb");
        fwrite(out, 1, sz, fo); fclose(fo);
        printf("Encoded: %d -> %ld (%.1f:1, ct=%s)\n", n, sz, (double)n/sz,
               out[8]==0?"G-RG-BG":"YCoCg-R");
    }

    free(rgb); free(out);
    return 0;
}
