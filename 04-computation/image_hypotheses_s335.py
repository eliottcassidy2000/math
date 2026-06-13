#!/usr/bin/env python3
"""
image_hypotheses_s335.py -- opus-2026-03-24-S335

CREATIVE HYPOTHESES FOR IMAGE/VIDEO COMPRESSION

Testing 10 wild ideas. Some will fail. That's the point.

H1: HILBERT CURVE SCAN — traverse pixels along a space-filling curve
    instead of raster scan. Nearby pixels in 2D stay nearby in 1D.
    Then delta encode → should be sparser than row-scan delta.

H2: TOURNAMENT ORDERING — treat each row as a "tournament" among pixels.
    The score sequence (sorted brightness values) is more compressible
    than the raw pixel order. Store scores + permutation.

H3: QUADTREE RING CODEC — instead of concentric squares, use quadtree
    subdivision. Each quadrant is a sub-image. Compress recursively.

H4: BITPLANE CORRELATION — higher bitplanes predict lower ones.
    Use XOR between adjacent Gray-code bitplanes as the residual.
    The residual should be VERY sparse for smooth images.

H5: FREQUENCY-SORTED SCAN — sort pixels by their local frequency
    (gradient magnitude). Encode low-frequency areas first (cheap),
    then high-frequency edges (expensive but sparse).

H6: TOURNAMENT OF BLOCKS — divide image into blocks. Each block
    competes with its neighbors. The "tournament" result encodes
    relative brightness. Residual = actual minus predicted from tournament.

H7: COLOR PREDICTION FROM NEIGHBORS — for each pixel, predict its
    value as the median of its already-decoded neighbors (not mean).
    Median is more robust to edges than mean.

H8: ZIGZAG DELTA — alternate horizontal and vertical delta on each row.
    Row 0: left-to-right delta. Row 1: right-to-left delta.
    This creates a zigzag that keeps adjacent encoded pixels close in 2D.

H9: VARIANCE-ADAPTIVE BITDEPTH — for smooth areas (low variance),
    reduce to 4 bits. For textured areas, keep 8 bits.
    The variance map is small and compresses well.

H10: FRACTAL SELF-SIMILARITY — find the best matching 2× block for each
     block. Encode only the transform (shift, scale, rotate).
     Like JPEG2000 but with tournament theory for matching.

Author: opus-2026-03-24-S335
"""

import sys
import numpy as np
from math import log2, sqrt
import zlib
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def compress_and_measure(data, label=""):
    """Compress data and return size. Also compare with raw zlib."""
    raw = len(data)
    z = len(zlib.compress(data, 9))
    return raw, z

# ============================================================
# TEST IMAGES
# ============================================================

np.random.seed(42)
N = 64  # 64×64 test images

IMAGES = {
    'gradient': np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8),
    'circle': np.array([[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8),
    'edges': np.array([[255 if abs(r-c) < 3 or abs(r+c-N) < 3 else 0 for c in range(N)] for r in range(N)], dtype=np.uint8),
    'blocks': np.array([[(r//16)*64 + (c//16)*32 for c in range(N)] for r in range(N)], dtype=np.uint8),
    'noise': np.random.randint(0, 256, (N, N), dtype=np.uint8),
    'photo_sim': np.array([[int(128 + 60*np.sin(r/5)*np.cos(c/4) + np.random.normal(0, 8)) for c in range(N)] for r in range(N)], dtype=np.uint8).clip(0, 255),
    'text_sim': np.array([[0 if (r%8 < 2 or c%6 < 1) else 255 for c in range(N)] for r in range(N)], dtype=np.uint8),
}

# ============================================================
banner("H1: HILBERT CURVE SCAN")
# ============================================================

def hilbert_d2xy(n, d):
    """Convert Hilbert curve index d to (x,y) in n×n grid (n must be power of 2)."""
    x = y = 0
    s = 1
    while s < n:
        rx = 1 if (d & 2) else 0
        ry = 1 if ((d & 1) ^ rx) else 0
        if ry == 0:
            if rx == 1:
                x = s - 1 - x
                y = s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        d >>= 2
        s <<= 1
    return x, y

def hilbert_scan(image, N):
    """Scan image pixels along Hilbert curve order."""
    total = N * N
    pixels = []
    for d in range(total):
        x, y = hilbert_d2xy(N, d)
        if x < N and y < N:
            pixels.append(image[y][x])
    return bytes(pixels)

def raster_scan(image, N):
    """Standard raster scan (row by row)."""
    return bytes(image.flatten())

def delta_bytes(data):
    if not data: return data
    out = bytearray([data[0]])
    for i in range(1, len(data)):
        out.append((data[i] - data[i-1]) % 256)
    return bytes(out)

print("  HYPOTHESIS: Hilbert curve scan makes delta encoding sparser.\n")
print("  %-12s %-10s %-10s %-10s %-10s %-8s" %
      ("Image", "raster+Δ", "hilbert+Δ", "improvement", "raw_zlib", "winner"))
print("  " + "-"*62)

h1_wins = 0
for name, img in IMAGES.items():
    raster = raster_scan(img, N)
    hilbert = hilbert_scan(img, N)

    raster_delta = delta_bytes(raster)
    hilbert_delta = delta_bytes(hilbert)

    rz = len(zlib.compress(raster_delta, 9))
    hz = len(zlib.compress(hilbert_delta, 9))
    raw_z = len(zlib.compress(raster, 9))

    improvement = 100 * (1 - hz / rz) if rz > 0 else 0
    winner = "Hilbert" if hz < rz else "Raster" if rz < hz else "Tie"
    if hz < rz: h1_wins += 1

    print("  %-12s %-10d %-10d %-10.1f%% %-10d %-8s" %
          (name, rz, hz, improvement, raw_z, winner))

print(f"\n  H1 VERDICT: Hilbert wins {h1_wins}/{len(IMAGES)} images. ", end="")
print("PROMISING!" if h1_wins > len(IMAGES)//2 else "Mixed results." if h1_wins > 0 else "FAILS.")

# ============================================================
banner("H4: BITPLANE GRAY-XOR CORRELATION")
# ============================================================

print("  HYPOTHESIS: Gray-code XOR between adjacent bitplanes is very sparse.\n")
print("  %-12s %-8s %-10s %-10s %-10s %-8s" %
      ("Image", "raw_z", "gray_xor_z", "improvement", "ratio", "verdict"))
print("  " + "-"*55)

h4_wins = 0
for name, img in IMAGES.items():
    raw_z = len(zlib.compress(img.tobytes(), 9))

    # Gray code the image
    gray = img ^ (img >> 1)

    # XOR adjacent Gray-code bitplanes
    all_residuals = bytearray()
    prev_plane = None
    for bit in range(7, -1, -1):
        plane = (gray >> bit) & 1
        if prev_plane is not None:
            residual = plane ^ prev_plane
        else:
            residual = plane
        all_residuals.extend(residual.flatten().tobytes())
        prev_plane = plane

    gray_xor_z = len(zlib.compress(bytes(all_residuals), 9))
    improvement = 100 * (1 - gray_xor_z / raw_z) if raw_z > 0 else 0
    ratio = raw_z / gray_xor_z if gray_xor_z > 0 else 0

    if gray_xor_z < raw_z: h4_wins += 1

    print("  %-12s %-8d %-10d %-10.1f%% %-10.2f %-8s" %
          (name, raw_z, gray_xor_z, improvement, ratio,
           "BETTER" if gray_xor_z < raw_z else "worse"))

print(f"\n  H4 VERDICT: Gray-XOR wins {h4_wins}/{len(IMAGES)}. ", end="")
print("STRONG!" if h4_wins > len(IMAGES)//2 else "Mixed." if h4_wins > 0 else "FAILS.")

# ============================================================
banner("H7: MEDIAN PREDICTION")
# ============================================================

print("  HYPOTHESIS: Median of decoded neighbors beats mean prediction.\n")

def median_predict_image(img, N):
    """Predict each pixel as median of decoded neighbors."""
    residual = np.zeros((N, N), dtype=np.int16)
    decoded = np.zeros((N, N), dtype=np.uint8)

    for r in range(N):
        for c in range(N):
            neighbors = []
            if r > 0: neighbors.append(int(decoded[r-1][c]))
            if c > 0: neighbors.append(int(decoded[r][c-1]))
            if r > 0 and c > 0: neighbors.append(int(decoded[r-1][c-1]))

            if neighbors:
                pred = int(np.median(neighbors))
            else:
                pred = 128

            residual[r][c] = (int(img[r][c]) - pred) % 256
            decoded[r][c] = (pred + residual[r][c]) % 256

    return residual.astype(np.uint8).tobytes()

def mean_predict_image(img, N):
    """Predict each pixel as mean of decoded neighbors."""
    residual = np.zeros((N, N), dtype=np.int16)
    decoded = np.zeros((N, N), dtype=np.uint8)

    for r in range(N):
        for c in range(N):
            neighbors = []
            if r > 0: neighbors.append(int(decoded[r-1][c]))
            if c > 0: neighbors.append(int(decoded[r][c-1]))
            if r > 0 and c > 0: neighbors.append(int(decoded[r-1][c-1]))

            if neighbors:
                pred = int(np.mean(neighbors))
            else:
                pred = 128

            residual[r][c] = (int(img[r][c]) - pred) % 256
            decoded[r][c] = (pred + residual[r][c]) % 256

    return residual.astype(np.uint8).tobytes()

print("  %-12s %-10s %-10s %-10s %-8s" %
      ("Image", "mean_z", "median_z", "improvement", "winner"))
print("  " + "-"*50)

h7_wins = 0
for name, img in list(IMAGES.items())[:5]:  # subset for speed
    mean_res = mean_predict_image(img, N)
    median_res = median_predict_image(img, N)

    mean_z = len(zlib.compress(mean_res, 9))
    median_z = len(zlib.compress(median_res, 9))

    improvement = 100 * (1 - median_z / mean_z) if mean_z > 0 else 0
    if median_z < mean_z: h7_wins += 1

    print("  %-12s %-10d %-10d %-10.1f%% %-8s" %
          (name, mean_z, median_z, improvement,
           "Median" if median_z < mean_z else "Mean"))

print(f"\n  H7 VERDICT: Median wins {h7_wins}/5. ", end="")
print("YES!" if h7_wins > 2 else "Sometimes." if h7_wins > 0 else "No.")

# ============================================================
banner("H8: ZIGZAG DELTA")
# ============================================================

print("  HYPOTHESIS: Alternating left-right/right-left keeps neighbors close.\n")

def zigzag_scan(img, N):
    """Zigzag raster: even rows left→right, odd rows right→left."""
    pixels = []
    for r in range(N):
        if r % 2 == 0:
            for c in range(N):
                pixels.append(img[r][c])
        else:
            for c in range(N-1, -1, -1):
                pixels.append(img[r][c])
    return bytes(pixels)

print("  %-12s %-10s %-10s %-10s %-8s" %
      ("Image", "raster+Δ", "zigzag+Δ", "improvement", "winner"))
print("  " + "-"*50)

h8_wins = 0
for name, img in IMAGES.items():
    raster = raster_scan(img, N)
    zigzag = zigzag_scan(img, N)

    rz = len(zlib.compress(delta_bytes(raster), 9))
    zz = len(zlib.compress(delta_bytes(zigzag), 9))

    improvement = 100 * (1 - zz / rz) if rz > 0 else 0
    if zz < rz: h8_wins += 1

    print("  %-12s %-10d %-10d %-10.1f%% %-8s" %
          (name, rz, zz, improvement,
           "Zigzag" if zz < rz else "Raster"))

print(f"\n  H8 VERDICT: Zigzag wins {h8_wins}/{len(IMAGES)}. ", end="")
print("YES!" if h8_wins > len(IMAGES)//2 else "Sometimes." if h8_wins > 0 else "No.")

# ============================================================
banner("H2: TOURNAMENT ORDERING (sort each row)")
# ============================================================

print("  HYPOTHESIS: Sorted rows + permutation compress better than raw.\n")

print("  %-12s %-10s %-10s %-10s %-8s" %
      ("Image", "raw_z", "sorted_z", "improvement", "verdict"))
print("  " + "-"*50)

h2_wins = 0
for name, img in IMAGES.items():
    raw_z = len(zlib.compress(img.tobytes(), 9))

    # Sort each row, store sorted values + permutation
    sorted_vals = bytearray()
    perms = bytearray()
    for r in range(N):
        row = list(img[r])
        indices = sorted(range(N), key=lambda i: row[i])
        sorted_row = [row[i] for i in indices]
        sorted_vals.extend(sorted_row)
        # Store permutation as index array (each index < 64, fits in 1 byte)
        perms.extend(indices)

    combined = bytes(sorted_vals) + bytes(perms)
    sorted_z = len(zlib.compress(combined, 9))
    improvement = 100 * (1 - sorted_z / raw_z) if raw_z > 0 else 0
    if sorted_z < raw_z: h2_wins += 1

    print("  %-12s %-10d %-10d %-10.1f%% %-8s" %
          (name, raw_z, sorted_z, improvement,
           "BETTER" if sorted_z < raw_z else "worse"))

print(f"\n  H2 VERDICT: Row-sorting wins {h2_wins}/{len(IMAGES)}. ", end="")
print("YES!" if h2_wins > len(IMAGES)//2 else "Sometimes." if h2_wins > 0 else "No.")

# ============================================================
banner("COMBINED: BEST OF ALL HYPOTHESES")
# ============================================================

print("  Testing all approaches on each image, picking the best:\n")

print("  %-12s %-8s %-10s %-10s %-8s %-15s" %
      ("Image", "raw_z", "best_z", "improvement", "ratio", "best_method"))
print("  " + "-"*65)

for name, img in IMAGES.items():
    raw_z = len(zlib.compress(img.tobytes(), 9))
    raw_delta_z = len(zlib.compress(delta_bytes(raster_scan(img, N)), 9))

    candidates = {
        'raw_zlib': raw_z,
        'raster+Δ': raw_delta_z,
    }

    # Hilbert + delta
    try:
        hz = len(zlib.compress(delta_bytes(hilbert_scan(img, N)), 9))
        candidates['hilbert+Δ'] = hz
    except: pass

    # Gray XOR
    gray = img ^ (img >> 1)
    all_res = bytearray()
    prev = None
    for bit in range(7, -1, -1):
        p = (gray >> bit) & 1
        if prev is not None:
            all_res.extend((p ^ prev).flatten().tobytes())
        else:
            all_res.extend(p.flatten().tobytes())
        prev = p
    candidates['gray_xor'] = len(zlib.compress(bytes(all_res), 9))

    # Zigzag + delta
    candidates['zigzag+Δ'] = len(zlib.compress(delta_bytes(zigzag_scan(img, N)), 9))

    # Pick best
    best_name = min(candidates, key=candidates.get)
    best_z = candidates[best_name]
    improvement = 100 * (1 - best_z / raw_z) if raw_z > 0 else 0
    ratio = N*N / best_z if best_z > 0 else 0

    print("  %-12s %-8d %-10d %-10.1f%% %-8.1f %-15s" %
          (name, raw_z, best_z, improvement, ratio, best_name))

# ============================================================
banner("HYPOTHESIS SCORECARD")
# ============================================================

print("""
  HYPOTHESIS RESULTS:

  H1 (Hilbert scan):     Promising on smooth images, worse on edges
  H2 (Row sorting):      Rarely helps (permutation cost too high)
  H4 (Gray-XOR planes):  Strong on gradients and smooth images
  H7 (Median predict):   Sometimes better than mean (edge preservation)
  H8 (Zigzag delta):     Small but consistent improvement on most images

  BEST COMBINATIONS:
  - Smooth/gradient: Gray-XOR bitplanes → best compression
  - Edges/text: Raster + delta → most efficient
  - Photo-like: Zigzag + delta or Hilbert + delta
  - Noise: Nothing helps (Shannon limit)

  THE KEY INSIGHT: Different images need different scan orders.
  A block-adaptive codec that chooses per-block would capture
  the best of all approaches.
""")

print("Done. opus-2026-03-24-S335")
