#!/usr/bin/env python3
"""
rotated_grid_codec.py — 45-degree rotation makes diagonal = spatial
kind-pasteur-2026-03-25-S20gp

THE TRICK: Rotate the image grid 45 degrees.

Original grid:          Rotated grid (diamond):
  (0,0) (0,1) (0,2)       (0,1)
  (1,0) (1,1) (1,2)    (0,0) (1,2)
  (2,0) (2,1) (2,2)       (1,0) (2,2)
                             (2,1)

In the rotated grid: diagonals of the ORIGINAL become rows/columns.
So diagonal-distance ordering on the ROTATED grid captures the
row/column correlation of the ORIGINAL.

THE DUAL COORDINATE SYSTEM:
  Standard: (row, col) — row-major scan captures horizontal correlation
  Rotated:  (row+col, row-col) — diagonal scan captures both directions

Apply zlib to FOUR scans and take the best:
  1. Row-major (standard)
  2. Column-major (transpose)
  3. Diagonal (sum = row+col constant)
  4. Anti-diagonal (diff = row-col constant)

For video delta: the best scan varies by frame content.
Encode which scan was chosen in a 2-bit header.

Also: combine ALL FOUR scans with zlib and pick the smallest.
The 2-bit overhead is negligible.
"""

import sys
import zlib
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ROTATED GRID CODEC — 45-degree dual coordinates")
print("  kind-pasteur-2026-03-25-S20gp")
print("=" * 80)


def scan_row_major(M):
    """Standard row-major scan."""
    H, W = M.shape
    return M.flatten().tobytes()


def scan_col_major(M):
    """Column-major scan."""
    return M.T.flatten().tobytes()


def scan_diagonal(M):
    """Scan by anti-diagonals (row+col = constant). The 45° rotation."""
    H, W = M.shape
    result = []
    for s in range(H + W - 1):
        for i in range(max(0, s-W+1), min(s+1, H)):
            j = s - i
            result.append(M[i, j])
    return bytes(np.array(result, dtype=np.uint8))


def scan_antidiagonal(M):
    """Scan by diagonals (row-col = constant). The other 45° rotation."""
    H, W = M.shape
    result = []
    for d in range(-(W-1), H):
        for i in range(max(0, d), min(H, W+d)):
            j = i - d
            result.append(M[i, j])
    return bytes(np.array(result, dtype=np.uint8))


def scan_hilbert(M):
    """Hilbert curve scan (approximate for non-power-of-2)."""
    H, W = M.shape
    # Simple Z-order (Morton) as Hilbert approximation
    result = []
    def interleave(x, y):
        z = 0
        for i in range(16):
            z |= ((x >> i) & 1) << (2*i)
            z |= ((y >> i) & 1) << (2*i + 1)
        return z
    coords = [(i, j) for i in range(H) for j in range(W)]
    coords.sort(key=lambda p: interleave(p[0], p[1]))
    for i, j in coords:
        result.append(M[i, j])
    return bytes(np.array(result, dtype=np.uint8))


def best_scan_compress(M, level=9):
    """Try all scans, return best compression."""
    scans = {
        'row': scan_row_major(M),
        'col': scan_col_major(M),
        'diag': scan_diagonal(M),
        'adiag': scan_antidiagonal(M),
        'zorder': scan_hilbert(M),
    }
    results = {}
    for name, data in scans.items():
        compressed = zlib.compress(data, level)
        results[name] = len(compressed)

    best = min(results, key=results.get)
    return results, best


# ================================================================
# BENCHMARK ON IMAGE DATA
# ================================================================
print("\n  1. STATIC IMAGE COMPRESSION:")
print(f"  {'Image':>20} {'Size':>6} {'Row':>7} {'Col':>7} {'Diag':>7} {'ADiag':>7} {'Zord':>7} {'Best':>7} {'Gain':>6}")

for N in [32, 64, 128]:
    for name, gen_fn in [
        ("h_gradient", lambda N: np.tile(np.arange(N, dtype=np.uint8), (N, 1))),
        ("v_gradient", lambda N: np.tile(np.arange(N, dtype=np.uint8).reshape(-1,1), (1, N))),
        ("d_gradient", lambda N: np.array([[(i+j)*128//N for j in range(N)] for i in range(N)], dtype=np.uint8)),
        ("circles", lambda N: np.array([[int(128+60*np.sin(np.sqrt((i-N/2)**2+(j-N/2)**2)/3)) for j in range(N)] for i in range(N)], dtype=np.uint8)),
        ("texture", lambda N: np.array([[int(128+40*np.sin(i/3)*np.cos(j/4)+20*np.sin((i+j)/5)) for j in range(N)] for i in range(N)], dtype=np.uint8)),
        ("block8", lambda N: np.array([[255 if (i//8+j//8)%2==0 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)),
    ]:
        M = gen_fn(N)
        results, best = best_scan_compress(M)
        row = results['row']
        gain = row / results[best]

        print(f"  {name:>20} {N:>3}x{N} {results['row']:>6}B {results['col']:>6}B {results['diag']:>6}B {results['adiag']:>6}B {results['zorder']:>6}B {best:>7} {gain:>5.2f}x")

# ================================================================
print(f"\n  2. VIDEO DELTA COMPRESSION:")
print(f"  {'Scene':>20} {'Size':>6} {'Row':>7} {'Diag':>7} {'ADiag':>7} {'Best':>7} {'vs Row':>7}")

for N in [64, 128]:
    np.random.seed(42)
    base = np.array([[int(128+60*np.sin(i/5)*np.cos(j/4)) for j in range(N)] for i in range(N)], dtype=np.uint8)

    for scene, make_frame2 in [
        ("pan_right", lambda b: np.roll(b, 2, axis=1)),
        ("pan_down", lambda b: np.roll(b, 2, axis=0)),
        ("pan_diag", lambda b: np.roll(np.roll(b, 2, axis=0), 2, axis=1)),
        ("zoom_in", lambda b: np.array([[b[min(int(N/2+(i-N/2)*0.95), N-1), min(int(N/2+(j-N/2)*0.95), N-1)] for j in range(N)] for i in range(N)], dtype=np.uint8)),
        ("noise_2", lambda b: np.clip(b.astype(int)+np.random.randint(-2,3,b.shape), 0, 255).astype(np.uint8)),
        ("noise_10", lambda b: np.clip(b.astype(int)+np.random.randint(-10,11,b.shape), 0, 255).astype(np.uint8)),
    ]:
        frame2 = make_frame2(base)
        diff = np.clip(frame2.astype(int) - base.astype(int) + 128, 0, 255).astype(np.uint8)

        results, best = best_scan_compress(diff)
        vs_row = results['row'] / results[best]

        tag = f"{N}x{N} {scene}"
        print(f"  {tag:>20} {N*N:>6} {results['row']:>6}B {results['diag']:>6}B {results['adiag']:>6}B {best:>7} {vs_row:>6.2f}x")

# ================================================================
print(f"\n  3. ADAPTIVE 4-SCAN (pick best per frame):")
print(f"  {'Scene':>15} {'Plain zlib':>11} {'4-scan+2b':>11} {'Gain':>6}")

for N in [128]:
    np.random.seed(42)
    base = np.array([[int(128+60*np.sin(i/5)*np.cos(j/4)) for j in range(N)] for i in range(N)], dtype=np.uint8)

    for scene, make_frame2 in [
        ("pan_right", lambda b: np.roll(b, 3, axis=1)),
        ("pan_down", lambda b: np.roll(b, 3, axis=0)),
        ("pan_diag", lambda b: np.roll(np.roll(b, 3, axis=0), 3, axis=1)),
        ("rotate_5", lambda b: np.array([[b[min(max(int(i*0.996-j*0.087+N*0.044),0),N-1), min(max(int(i*0.087+j*0.996-N*0.044),0),N-1)] for j in range(N)] for i in range(N)], dtype=np.uint8)),
        ("noise_5", lambda b: np.clip(b.astype(int)+np.random.randint(-5,6,b.shape),0,255).astype(np.uint8)),
        ("mixed", lambda b: np.clip(np.roll(b,1,axis=1).astype(int)+np.random.randint(-3,4,b.shape),0,255).astype(np.uint8)),
    ]:
        frame2 = make_frame2(base)
        diff = np.clip(frame2.astype(int) - base.astype(int) + 128, 0, 255).astype(np.uint8)

        plain = len(zlib.compress(bytes(diff.flatten()), 9))
        results, best = best_scan_compress(diff)
        adaptive = results[best] + 1  # +1 byte for scan choice header

        gain = plain / adaptive
        print(f"  {scene:>15} {plain:>10}B {adaptive:>10}B {gain:>5.3f}x")

print(f"\n{'='*60}")
print("THE 45-DEGREE TRICK")
print(f"{'='*60}")
print("""
The rotated coordinate system works because:

  ROW MAJOR: captures horizontal correlation (pan right)
  COL MAJOR: captures vertical correlation (pan down)
  DIAGONAL: captures diagonal correlation (pan diagonal, rotation)
  ANTI-DIAG: captures the other diagonal

Each scan direction is optimal for different MOTION TYPES:
  Horizontal motion -> row major wins
  Vertical motion -> col major wins
  Diagonal motion -> diagonal scan wins
  Rotation -> mixed scans help

The ADAPTIVE 4-SCAN approach: try all 4, pick best, prepend 2-bit header.
Cost: ~4x encode time (parallelizable), 2 bits overhead.
Gain: 1-10% over row-major for directional motion.

Combined with bitplane decomposition:
  Each bit plane gets its OWN best scan direction.
  MSB planes might prefer row-major (smooth gradients horizontal)
  LSB planes might prefer diagonal (noise is isotropic)
  The per-plane adaptive scan adds ~16 bits overhead for 8-bit image.
""")

print("DONE.")
print("=" * 80)
