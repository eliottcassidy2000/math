/*
 * TIC FAST — MED prediction kernel in C, called from Python via ctypes.
 * This handles the hot loop; Python handles zlib + file I/O.
 *
 * Compile: gcc -O3 -shared -o tic_fast.dll tic_fast.c  (Windows)
 *          gcc -O3 -shared -fPIC -o tic_fast.so tic_fast.c  (Linux)
 *
 * kind-pasteur-2026-03-25-S22
 */

#include <string.h>
#include <stdlib.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

static inline int med_pred(int a, int b, int c) {
    int mn = a < b ? a : b;
    int mx = a > b ? a : b;
    if (c >= mx) return mn;
    if (c <= mn) return mx;
    return a + b - c;
}

/*
 * Encode a plane with MED prediction.
 * plane: input pixels [h*w], uint8
 * residuals: output [h*w], uint8
 */
EXPORT void encode_med(const unsigned char *plane, unsigned char *residuals,
                       int h, int w) {
    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int d = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            residuals[idx] = (plane[idx] - med_pred(a, b, d)) & 0xFF;
        }
    }
}

/*
 * Decode MED residuals back to plane.
 */
EXPORT void decode_med(const unsigned char *residuals, unsigned char *plane,
                       int h, int w) {
    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int idx = r * w + c;
            int a = (c > 0) ? plane[idx - 1] : 0;
            int b = (r > 0) ? plane[idx - w] : 0;
            int d = (r > 0 && c > 0) ? plane[idx - w - 1] : 0;
            plane[idx] = (residuals[idx] + med_pred(a, b, d)) & 0xFF;
        }
    }
}

/*
 * Color decorrelation: RGB → G, R-G, B-G
 */
EXPORT void rgb_to_grd(const unsigned char *rgb, unsigned char *g,
                       unsigned char *rg, unsigned char *bg, int n) {
    for (int i = 0; i < n; i++) {
        int rv = rgb[i * 3];
        int gv = rgb[i * 3 + 1];
        int bv = rgb[i * 3 + 2];
        g[i] = gv;
        rg[i] = (rv - gv) & 0xFF;
        bg[i] = (bv - gv) & 0xFF;
    }
}

/*
 * Inverse color decorrelation
 */
EXPORT void grd_to_rgb(const unsigned char *g, const unsigned char *rg,
                       const unsigned char *bg, unsigned char *rgb, int n) {
    for (int i = 0; i < n; i++) {
        int gv = g[i];
        rgb[i * 3] = (rg[i] + gv) & 0xFF;
        rgb[i * 3 + 1] = gv;
        rgb[i * 3 + 2] = (bg[i] + gv) & 0xFF;
    }
}

/*
 * Full encode pipeline: RGB → G-RG-BG → MED residuals for each plane.
 * Returns 3 residual planes concatenated.
 */
EXPORT void encode_rgb_full(const unsigned char *rgb, unsigned char *out,
                            int h, int w) {
    int n = h * w;
    unsigned char *g = out;
    unsigned char *rg_out = out + n;
    unsigned char *bg_out = out + 2 * n;

    /* Temp buffers for color transform */
    unsigned char *g_raw = (unsigned char*)malloc(n);
    unsigned char *rg_raw = (unsigned char*)malloc(n);
    unsigned char *bg_raw = (unsigned char*)malloc(n);

    rgb_to_grd(rgb, g_raw, rg_raw, bg_raw, n);
    encode_med(g_raw, g, h, w);
    encode_med(rg_raw, rg_out, h, w);
    encode_med(bg_raw, bg_out, h, w);

    free(g_raw); free(rg_raw); free(bg_raw);
}

/*
 * Full decode pipeline
 */
EXPORT void decode_rgb_full(const unsigned char *in, unsigned char *rgb,
                            int h, int w) {
    int n = h * w;
    const unsigned char *g_res = in;
    const unsigned char *rg_res = in + n;
    const unsigned char *bg_res = in + 2 * n;

    unsigned char *g = (unsigned char*)malloc(n);
    unsigned char *rg = (unsigned char*)malloc(n);
    unsigned char *bg = (unsigned char*)malloc(n);

    decode_med(g_res, g, h, w);
    decode_med(rg_res, rg, h, w);
    decode_med(bg_res, bg, h, w);
    grd_to_rgb(g, rg, bg, rgb, n);

    free(g); free(rg); free(bg);
}
