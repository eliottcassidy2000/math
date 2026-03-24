#!/usr/bin/env python3
"""
hic_v3_s336.py — Hilbert Image Codec v3: Per-Row Adaptive + Extended Filters

Like PNG's per-row filter selection, but with 8 filters instead of 5:
  0: None (raw)
  1: Sub (left neighbor)
  2: Up (above neighbor)
  3: Average (mean of left + above)
  4: Paeth (PNG predictor)
  5: Diagonal (upper-left neighbor only)
  6: Gradient (a + b - c, clamped — like Paeth but without abs comparison)
  7: Median (median of left, above, upper-left)

For each row: try all 8 filters, pick the one that gives smallest
zlib output for that row. Store: 1 byte per row + filtered data.

This is essentially PNG with MORE filter options.
The extra filters (5-7) capture diagonal and median correlations
that PNG misses.

Author: opus-2026-03-24-S336
"""

import sys
import numpy as np
import zlib
import struct
import time
sys.stdout.reconfigure(line_buffering=True)

def paeth_val(a, b, c):
    a, b, c = int(a), int(b), int(c)
    p = a + b - c
    pa, pb, pc = abs(p-a), abs(p-b), abs(p-c)
    if pa <= pb and pa <= pc: return a
    elif pb <= pc: return b
    return c

def apply_filter(row, prev_row, filt):
    """Apply filter to a row given the previous row. Returns filtered bytes."""
    W = len(row)
    out = bytearray(W)

    for c in range(W):
        a = int(row[c-1]) if c > 0 else 0      # left
        b = int(prev_row[c]) if prev_row is not None else 0  # above
        d = int(prev_row[c-1]) if prev_row is not None and c > 0 else 0  # upper-left

        x = int(row[c])

        if filt == 0:    # None
            out[c] = x
        elif filt == 1:  # Sub
            out[c] = (x - a) % 256
        elif filt == 2:  # Up
            out[c] = (x - b) % 256
        elif filt == 3:  # Average
            out[c] = (x - (a + b) // 2) % 256
        elif filt == 4:  # Paeth
            out[c] = (x - paeth_val(a, b, d)) % 256
        elif filt == 5:  # Diagonal (upper-left only)
            out[c] = (x - d) % 256
        elif filt == 6:  # Gradient (clamped linear prediction)
            pred = max(0, min(255, a + b - d))
            out[c] = (x - pred) % 256
        elif filt == 7:  # Median
            pred = int(np.median([a, b, d]))
            out[c] = (x - pred) % 256

    return bytes(out)

def reverse_filter(filtered, prev_row, filt):
    """Reverse a filter to get original row."""
    W = len(filtered)
    row = bytearray(W)

    for c in range(W):
        a = int(row[c-1]) if c > 0 else 0
        b = int(prev_row[c]) if prev_row is not None else 0
        d = int(prev_row[c-1]) if prev_row is not None and c > 0 else 0

        f = filtered[c]

        if filt == 0:
            row[c] = f
        elif filt == 1:
            row[c] = (f + a) % 256
        elif filt == 2:
            row[c] = (f + b) % 256
        elif filt == 3:
            row[c] = (f + (a + b) // 2) % 256
        elif filt == 4:
            row[c] = (f + paeth_val(a, b, d)) % 256
        elif filt == 5:
            row[c] = (f + d) % 256
        elif filt == 6:
            pred = max(0, min(255, a + b - d))
            row[c] = (f + pred) % 256
        elif filt == 7:
            pred = int(np.median([a, b, d]))
            row[c] = (f + pred) % 256

    return bytes(row)

N_FILTERS = 8

def compress_image(image):
    """Compress using per-row adaptive filtering."""
    H, W = image.shape

    # For each row: try all filters, pick best
    filter_choices = bytearray(H)
    all_filtered = bytearray()

    prev_row = None
    for r in range(H):
        row = image[r]
        best_filt = 0
        best_filtered = apply_filter(row, prev_row, 0)
        best_score = sum(min(b, 256-b) for b in best_filtered)

        for filt in range(N_FILTERS):
            filtered = apply_filter(row, prev_row, filt)
            # Score: sum of abs(filtered) — proxy for compressibility
            score = sum(min(b, 256-b) for b in filtered)
            if score < best_score:
                best_score = score
                best_filt = filt
                best_filtered = filtered

        filter_choices[r] = best_filt
        all_filtered.extend(best_filtered)
        prev_row = row

    # Pack: header + filter choices + filtered data, all compressed together
    header = struct.pack('<HH', H, W)
    payload = bytes(filter_choices) + bytes(all_filtered)
    return header + zlib.compress(payload, 9)

def decompress_image(data):
    """Decompress per-row adaptive filtered image."""
    H, W = struct.unpack_from('<HH', data, 0)
    payload = zlib.decompress(data[4:])

    filter_choices = list(payload[:H])
    filtered_data = payload[H:]

    image = np.zeros((H, W), dtype=np.uint8)
    prev_row = None

    for r in range(H):
        filt = filter_choices[r]
        row_data = filtered_data[r*W:(r+1)*W]
        row = reverse_filter(row_data, prev_row, filt)
        image[r] = np.frombuffer(row, dtype=np.uint8)
        prev_row = image[r]

    return image

# ============================================================
# BENCHMARK
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 70)
    print("  HIC v3 — Per-Row Adaptive + 8 Filters")
    print("=" * 70)
    print()

    tests = {}
    for N in [64, 128, 256]:
        tests[f'gradient_{N}'] = np.array(
            [[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f'circle_{N}'] = np.array(
            [[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8))))
              for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f'smooth_{N}'] = np.array(
            [[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)], dtype=np.uint8).clip(0,255)

    tests['noise_128'] = np.random.randint(0, 256, (128, 128), dtype=np.uint8)
    tests['edges_128'] = np.array(
        [[255 if abs(r-c)<3 or abs(r+c-128)<3 else 0 for c in range(128)] for r in range(128)], dtype=np.uint8)
    tests['photo_128'] = np.array(
        [[int(128+60*np.sin(r/5)*np.cos(c/4)+np.random.normal(0,15))
          for c in range(128)] for r in range(128)], dtype=np.uint8).clip(0,255)

    print("  %-20s %-8s %-8s %-8s %-6s %-6s %-20s" %
          ("Image", "Raw", "HICv3", "zlib", "v3/z", "ratio", "filters_used"))
    print("  " + "-"*75)

    wins = 0; total = 0
    for name, img in sorted(tests.items()):
        raw = img.size
        total += 1

        t0 = time.time()
        v3 = compress_image(img)
        t_c = time.time() - t0

        recovered = decompress_image(v3)
        ok = "✓" if np.array_equal(img, recovered) else "✗"

        z = len(zlib.compress(img.tobytes(), 9))
        vs = z / len(v3) if len(v3) > 0 else 0
        ratio = raw / len(v3) if len(v3) > 0 else 0

        if len(v3) < z: wins += 1

        # Count filter usage
        payload = zlib.decompress(v3[4:])
        H = img.shape[0]
        filts = list(payload[:H])
        filt_counts = {}
        for f in filts:
            fname = ['none','sub','up','avg','paeth','diag','grad','median'][f]
            filt_counts[fname] = filt_counts.get(fname, 0) + 1
        top = sorted(filt_counts.items(), key=lambda x: -x[1])[:3]
        filt_str = " ".join(f"{n}:{c}" for n, c in top)

        print("  %-20s %-8d %-7d%s %-8d %-6.2f %-6.1f %-20s" %
              (name, raw, len(v3), ok, z, vs, ratio, filt_str))

    print(f"\n  Wins vs zlib: {wins}/{total}")
    print(f"  All lossless: {'✓' if all(np.array_equal(t, decompress_image(compress_image(t))) for t in tests.values()) else '✗'}")
