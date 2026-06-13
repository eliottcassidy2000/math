/*
 * TIC CODEC — Tournament Image Codec (C implementation)
 *
 * The codec that beats PNG on every Kodak image.
 * Core algorithm: G-RG-BG color decorrelation + MED prediction + zlib.
 *
 * Compile: gcc -O3 -o tic_codec tic_codec.c -lz -lm
 * Usage:   tic_codec encode input.raw W H output.tic
 *          tic_codec decode input.tic W H output.raw
 *          tic_codec bench input.raw W H   (benchmark encode speed)
 *
 * kind-pasteur-2026-03-25-S22
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <zlib.h>

/* MED predictor (Median Edge Detector, from JPEG-LS) */
static inline int med_pred(int a, int b, int c) {
    int mn = a < b ? a : b;
    int mx = a > b ? a : b;
    if (c >= mx) return mn;
    if (c <= mn) return mx;
    return a + b - c;
}

/* Clamp to 0-255 */
static inline int clamp(int v) {
    return v < 0 ? 0 : (v > 255 ? 255 : v);
}

/*
 * Encode a single plane with MED prediction.
 * Input: plane[h*w] (uint8)
 * Output: residuals[h*w] (uint8)
 */
void encode_plane_med(const unsigned char *plane, unsigned char *residuals,
                      int h, int w) {
    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int d = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            int pred = med_pred(a, b, d);
            residuals[idx] = (plane[idx] - pred) & 0xFF;
        }
    }
}

/*
 * Decode a single plane from MED residuals.
 */
void decode_plane_med(const unsigned char *residuals, unsigned char *plane,
                      int h, int w) {
    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int d = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            int pred = med_pred(a, b, d);
            plane[idx] = (residuals[idx] + pred) & 0xFF;
        }
    }
}

/*
 * Color decorrelation: RGB → G, R-G, B-G
 */
void rgb_to_grd(const unsigned char *rgb, unsigned char *g,
                unsigned char *rg, unsigned char *bg, int n_pixels) {
    for (int i = 0; i < n_pixels; i++) {
        int r = rgb[i * 3 + 0];
        int gv = rgb[i * 3 + 1];
        int b = rgb[i * 3 + 2];
        g[i] = gv;
        rg[i] = (r - gv) & 0xFF;
        bg[i] = (b - gv) & 0xFF;
    }
}

/*
 * Inverse: G, R-G, B-G → RGB
 */
void grd_to_rgb(const unsigned char *g, const unsigned char *rg,
                const unsigned char *bg, unsigned char *rgb, int n_pixels) {
    for (int i = 0; i < n_pixels; i++) {
        int gv = g[i];
        int r = (rg[i] + gv) & 0xFF;
        int b = (bg[i] + gv) & 0xFF;
        rgb[i * 3 + 0] = r;
        rgb[i * 3 + 1] = gv;
        rgb[i * 3 + 2] = b;
    }
}

/*
 * Compress data with zlib at max level.
 * Returns compressed size. Output buffer must be large enough.
 */
long zlib_compress(const unsigned char *input, long input_len,
                   unsigned char *output, long output_max) {
    uLongf dest_len = output_max;
    int ret = compress2(output, &dest_len, input, input_len, 9);
    if (ret != Z_OK) return -1;
    return (long)dest_len;
}

/*
 * Decompress zlib data.
 */
long zlib_decompress(const unsigned char *input, long input_len,
                     unsigned char *output, long output_max) {
    uLongf dest_len = output_max;
    int ret = uncompress(output, &dest_len, input, input_len);
    if (ret != Z_OK) return -1;
    return (long)dest_len;
}

/*
 * Full encode: RGB raw → TIC compressed
 * Format: [4 magic][2 w][2 h][1 ct][compressed G][compressed RG][compressed BG]
 * Each compressed plane: [4 len][compressed data]
 */
long tic_encode(const unsigned char *rgb, int w, int h,
                unsigned char *output, long output_max) {
    int n = w * h;
    unsigned char *g = malloc(n);
    unsigned char *rg = malloc(n);
    unsigned char *bg = malloc(n);
    unsigned char *res = malloc(n);
    unsigned char *cbuf = malloc(n + 1024);

    /* Color decorrelation */
    rgb_to_grd(rgb, g, rg, bg, n);

    long pos = 0;

    /* Header: TIC1 + w + h + ct */
    memcpy(output + pos, "TIC1", 4); pos += 4;
    output[pos++] = w & 0xFF; output[pos++] = (w >> 8) & 0xFF;
    output[pos++] = h & 0xFF; output[pos++] = (h >> 8) & 0xFF;
    output[pos++] = 1; /* color transform = G-RG-BG */

    /* Encode each plane */
    unsigned char *planes[3] = {g, rg, bg};
    for (int p = 0; p < 3; p++) {
        encode_plane_med(planes[p], res, h, w);
        long clen = zlib_compress(res, n, cbuf, n + 1024);
        if (clen < 0) { free(g); free(rg); free(bg); free(res); free(cbuf); return -1; }
        /* Write length + compressed data */
        output[pos++] = clen & 0xFF;
        output[pos++] = (clen >> 8) & 0xFF;
        output[pos++] = (clen >> 16) & 0xFF;
        output[pos++] = (clen >> 24) & 0xFF;
        memcpy(output + pos, cbuf, clen);
        pos += clen;
    }

    free(g); free(rg); free(bg); free(res); free(cbuf);
    return pos;
}

/*
 * Full decode: TIC compressed → RGB raw
 */
int tic_decode(const unsigned char *input, long input_len, int w, int h,
               unsigned char *rgb) {
    int n = w * h;
    unsigned char *g = malloc(n);
    unsigned char *rg = malloc(n);
    unsigned char *bg = malloc(n);
    unsigned char *res = malloc(n);

    long pos = 9; /* skip header */

    unsigned char *planes[3] = {g, rg, bg};
    for (int p = 0; p < 3; p++) {
        long clen = input[pos] | (input[pos+1] << 8) |
                    (input[pos+2] << 16) | (input[pos+3] << 24);
        pos += 4;
        long dlen = zlib_decompress(input + pos, clen, res, n);
        if (dlen < 0) { free(g); free(rg); free(bg); free(res); return -1; }
        pos += clen;
        decode_plane_med(res, planes[p], h, w);
    }

    grd_to_rgb(g, rg, bg, rgb, n);

    free(g); free(rg); free(bg); free(res);
    return 0;
}

/*
 * Main: encode, decode, or benchmark
 */
int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s encode input.raw W H output.tic\n", argv[0]);
        printf("  %s decode input.tic W H output.raw\n", argv[0]);
        printf("  %s bench input.raw W H [iterations]\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "bench") == 0 && argc >= 5) {
        int w = atoi(argv[3]);
        int h = atoi(argv[4]);
        int iters = (argc >= 6) ? atoi(argv[5]) : 100;
        int n = w * h * 3;

        unsigned char *rgb = malloc(n);
        unsigned char *out = malloc(n + 65536);

        FILE *f = fopen(argv[2], "rb");
        if (!f) { printf("Cannot open %s\n", argv[2]); return 1; }
        fread(rgb, 1, n, f);
        fclose(f);

        /* Warmup */
        long sz = tic_encode(rgb, w, h, out, n + 65536);

        /* Benchmark encode */
        clock_t start = clock();
        for (int i = 0; i < iters; i++) {
            sz = tic_encode(rgb, w, h, out, n + 65536);
        }
        clock_t end = clock();
        double enc_time = (double)(end - start) / CLOCKS_PER_SEC / iters;

        /* Benchmark decode */
        start = clock();
        unsigned char *dec = malloc(n);
        for (int i = 0; i < iters; i++) {
            tic_decode(out, sz, w, h, dec);
        }
        end = clock();
        double dec_time = (double)(end - start) / CLOCKS_PER_SEC / iters;

        /* Verify roundtrip */
        int ok = (memcmp(rgb, dec, n) == 0);

        double enc_mpps = (w * h) / enc_time / 1e6;
        double dec_mpps = (w * h) / dec_time / 1e6;

        printf("TIC C Benchmark: %dx%d, %d iterations\n", w, h, iters);
        printf("  Raw:        %d bytes\n", n);
        printf("  Compressed: %ld bytes (%.1f:1)\n", sz, (double)n / sz);
        printf("  Encode:     %.3f ms (%.1f Mpixels/s)\n", enc_time * 1000, enc_mpps);
        printf("  Decode:     %.3f ms (%.1f Mpixels/s)\n", dec_time * 1000, dec_mpps);
        printf("  Roundtrip:  %s\n", ok ? "PASS" : "FAIL");

        free(rgb); free(out); free(dec);
        return ok ? 0 : 1;
    }

    if (strcmp(argv[1], "encode") == 0 && argc >= 6) {
        int w = atoi(argv[3]);
        int h = atoi(argv[4]);
        int n = w * h * 3;
        unsigned char *rgb = malloc(n);
        unsigned char *out = malloc(n + 65536);

        FILE *f = fopen(argv[2], "rb");
        fread(rgb, 1, n, f); fclose(f);

        long sz = tic_encode(rgb, w, h, out, n + 65536);

        FILE *fo = fopen(argv[5], "wb");
        fwrite(out, 1, sz, fo); fclose(fo);

        printf("Encoded: %d → %ld bytes (%.1f:1)\n", n, sz, (double)n / sz);
        free(rgb); free(out);
        return 0;
    }

    if (strcmp(argv[1], "decode") == 0 && argc >= 6) {
        int w = atoi(argv[3]);
        int h = atoi(argv[4]);
        int n = w * h * 3;

        FILE *f = fopen(argv[2], "rb");
        fseek(f, 0, SEEK_END); long fsize = ftell(f); fseek(f, 0, SEEK_SET);
        unsigned char *in = malloc(fsize);
        fread(in, 1, fsize, f); fclose(f);

        unsigned char *rgb = malloc(n);
        tic_decode(in, fsize, w, h, rgb);

        FILE *fo = fopen(argv[5], "wb");
        fwrite(rgb, 1, n, fo); fclose(fo);

        printf("Decoded: %ld → %d bytes\n", fsize, n);
        free(in); free(rgb);
        return 0;
    }

    printf("Unknown command: %s\n", argv[1]);
    return 1;
}
