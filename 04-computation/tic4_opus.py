#!/usr/bin/env python3
"""
tic4.py — Tournament Image Codec v4: Spiral + Checkerboard + Zigzag Scans

INNOVATIONS beyond v3:
  1. SPIRAL scan: continuous center-outward path (1D locality preserved)
  2. CHECKERBOARD scan: two-phase — all "white" squares then "black" squares
     Phase 2 has ALL 4 cardinal neighbors from phase 1 as context
  3. ZIGZAG scan: diagonal zigzag (like JPEG DCT coefficient order)
  4. HILBERT scan: space-filling curve that preserves 2D locality in 1D

KEY INSIGHT: Checkerboard is a MULTIRESOLUTION approach.
  Phase 1: half the pixels, sub-sampled, compresses like a half-res image
  Phase 2: interpolation residuals from all 4 neighbors → near-zero for smooth regions

Author: opus-2026-03-25-S342
"""

import sys, struct, zlib, lzma, io, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tic3 import (
    encode as _enc3, decode as _dec3,
    VFILTERS_2D, N_FILTERS_2D, VFILTERS_1D, N_FILTERS_1D,
    _score, _fast_compress, _best_backend, _decompress,
    _unfilter_row_2d, _unfilter_line_1d,
    CT_FWD, CT_INV,
    _predict_ring_pixel,
    _scan_raster, _scan_transpose, _scan_diagonal, _scan_antidiag, _scan_ring,
    MAGIC as MAGIC_V3, png_size,
)

__version__ = "4.0.0"
MAGIC = b'TIC4'

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# NEW SCAN ORDERS
# ============================================================

def _scan_spiral(H, W):
    """Continuous clockwise spiral from center outward.
    Returns a SINGLE list of (r,c) coords — one continuous path.
    Unlike ring scan, spiral is a single 1D line with adjacent-pixel locality."""
    cr, cc = H // 2, W // 2
    visited = set()
    path = []

    # Start at center
    if 0 <= cr < H and 0 <= cc < W:
        path.append((cr, cc))
        visited.add((cr, cc))

    # Spiral outward: right, down, left, up in increasing sizes
    dr = [0, 1, 0, -1]  # direction vectors
    dc = [1, 0, -1, 0]
    r, c = cr, cc
    step = 1
    d = 0  # direction index
    while len(path) < H * W:
        for _ in range(2):  # each step size is used twice
            for _ in range(step):
                r += dr[d]
                c += dc[d]
                if 0 <= r < H and 0 <= c < W and (r, c) not in visited:
                    path.append((r, c))
                    visited.add((r, c))
                if len(path) >= H * W:
                    break
            d = (d + 1) % 4
            if len(path) >= H * W:
                break
        step += 1
        if step > H + W:
            break

    # Fill any missed pixels (edge cases at non-square boundaries)
    if len(path) < H * W:
        for r in range(H):
            for c in range(W):
                if (r, c) not in visited:
                    path.append((r, c))
                    visited.add((r, c))

    return path


def _scan_checkerboard(H, W):
    """Two-phase checkerboard scan.
    Phase 1: all (r,c) where (r+c) % 2 == 0 (raster order within phase)
    Phase 2: all (r,c) where (r+c) % 2 == 1 (raster order within phase)

    Phase 2 pixels have all 4 cardinal neighbors from Phase 1 as context."""
    phase1 = []
    phase2 = []
    for r in range(H):
        for c in range(W):
            if (r + c) % 2 == 0:
                phase1.append((r, c))
            else:
                phase2.append((r, c))
    return phase1, phase2


def _scan_zigzag(H, W):
    """JPEG-style diagonal zigzag scan. Traverses anti-diagonals
    alternating direction (top-right to bottom-left, then reverse).
    Returns a single list of (r,c) coords."""
    path = []
    for k in range(H + W - 1):
        if k % 2 == 0:
            # Even diagonals: go up-right
            r_start = min(k, H - 1)
            c_start = k - r_start
            r, c = r_start, c_start
            while r >= 0 and c < W:
                path.append((r, c))
                r -= 1; c += 1
        else:
            # Odd diagonals: go down-left
            c_start = min(k, W - 1)
            r_start = k - c_start
            r, c = r_start, c_start
            while c >= 0 and r < H:
                path.append((r, c))
                r += 1; c -= 1
    return path


def _hilbert_d2xy(n, d):
    """Convert Hilbert curve index d to (x,y) coords for n×n grid (n must be power of 2)."""
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


def _scan_hilbert(H, W):
    """Hilbert curve scan for the largest power-of-2 that fits, then fill remainder.
    Returns a single list of (r,c) coords."""
    n = 1
    while n * 2 <= min(H, W):
        n *= 2

    path = []
    visited = set()
    for d in range(n * n):
        x, y = _hilbert_d2xy(n, d)
        if 0 <= y < H and 0 <= x < W:
            path.append((y, x))
            visited.add((y, x))

    # Fill remaining pixels in raster order
    for r in range(H):
        for c in range(W):
            if (r, c) not in visited:
                path.append((r, c))
                visited.add((r, c))

    return path


# ============================================================
# CHECKERBOARD PREDICTION
# ============================================================

def _predict_checker_phase2(plane, r, c, H, W):
    """Predict a phase-2 pixel from its 4 cardinal neighbors (all phase-1).
    Uses median of available neighbors for robustness."""
    neighbors = []
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nr, nc = r + dr, c + dc
        if 0 <= nr < H and 0 <= nc < W:
            neighbors.append(int(plane[nr, nc]))
    if not neighbors:
        return 128
    if len(neighbors) <= 2:
        return sum(neighbors) // len(neighbors)
    # Median of neighbors — more robust than mean at edges
    neighbors.sort()
    mid = len(neighbors) // 2
    return neighbors[mid]


def _predict_checker_phase2_avg(plane, r, c, H, W):
    """Simple average of 4 cardinal neighbors."""
    total = count = 0
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nr, nc = r + dr, c + dc
        if 0 <= nr < H and 0 <= nc < W:
            total += int(plane[nr, nc]); count += 1
    return total // count if count > 0 else 128


# ============================================================
# SPIRAL/ZIGZAG/HILBERT PREDICTION
# ============================================================

def _predict_path_pixel(plane, r, c, known_set, H, W):
    """Predict from all known 8-neighbors (average). Used for path-based scans."""
    total = count = 0
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0: continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < H and 0 <= nc < W and (nr, nc) in known_set:
                total += int(plane[nr, nc]); count += 1
    return total // count if count > 0 else 128


# ============================================================
# ENCODE
# ============================================================

def encode(image):
    image = np.asarray(image, dtype=np.uint8)
    if image.ndim == 2:
        H, W = image.shape
        channels = 1
        planes_list = [([image], 0)]
    elif image.ndim == 3 and image.shape[2] == 3:
        H, W = image.shape[:2]
        channels = 3
        ct_scores = []
        all_planes = []
        for ci, ct in enumerate(CT_FWD):
            transformed = ct(image)
            pls = [transformed[:,:,c] for c in range(3)]
            all_planes.append(pls)
            ct_scores.append((sum(float(np.std(p)) for p in pls), ci))
        ct_scores.sort()
        top2 = [ci for _, ci in ct_scores[:2]]
        planes_list = [(all_planes[ci], ci) for ci in top2]
    else:
        raise ValueError(f"Unsupported shape: {image.shape}")

    all_candidates = []

    for planes, ct_id in planes_list:
        # ---- V3 STRATEGIES (raster, transpose, diagonal, ring, etc.) ----
        # Raster per-row adaptive
        buf = bytearray()
        for plane in planes:
            prev = None
            for y in range(H):
                row = plane[y]
                best_f, best_s, best_o = 0, W * 256, None
                for f in range(N_FILTERS_2D):
                    filt = VFILTERS_2D[f](row, prev, W)
                    s = _score(filt)
                    if s < best_s: best_s = s; best_f = f; best_o = filt
                buf.append(best_f); buf.extend(best_o.tobytes()); prev = row
        all_candidates.append((0, ct_id, bytes(buf)))

        # Whole-plane single filter (top 7)
        for f in range(min(N_FILTERS_2D, 7)):
            buf = bytearray()
            for plane in planes:
                prev = None
                for y in range(H):
                    filt = VFILTERS_2D[f](plane[y], prev, W)
                    buf.extend(filt.tobytes()); prev = plane[y]
            all_candidates.append((0x10 + f, ct_id, bytes(buf)))

        # Transpose per-row adaptive
        buf = bytearray()
        for plane in planes:
            pt = plane.T.copy(); prev = None
            for y in range(W):
                row = pt[y]
                best_f, best_s, best_o = 0, H * 256, None
                for f in range(N_FILTERS_2D):
                    filt = VFILTERS_2D[f](row, prev, H)
                    s = _score(filt)
                    if s < best_s: best_s = s; best_f = f; best_o = filt
                buf.append(best_f); buf.extend(best_o.tobytes()); prev = row
        all_candidates.append((1, ct_id, bytes(buf)))

        # Diagonal/anti-diagonal with global filters
        for scan_base, scan_fn in [(2, _scan_diagonal), (3, _scan_antidiag)]:
            lines = scan_fn(H, W)
            for filt_id in range(N_FILTERS_1D):
                buf = bytearray()
                for plane in planes:
                    for lc in lines:
                        vals = np.array([plane[r,c] for r,c in lc], dtype=np.uint8)
                        L = len(vals)
                        if L == 0: continue
                        filt = VFILTERS_1D[filt_id](vals, L)
                        buf.extend(filt.tobytes())
                all_candidates.append((scan_base * 16 + filt_id, ct_id, bytes(buf)))

        # Ring scan with avg-neighbor prediction
        ring_lines = _scan_ring(H, W)
        buf = bytearray()
        for plane in planes:
            known = set()
            for ring_pixels in ring_lines:
                for r, c in ring_pixels:
                    pred = _predict_ring_pixel(plane, r, c, known, H, W)
                    buf.append((int(plane[r, c]) - pred) & 0xFF)
                for r, c in ring_pixels:
                    known.add((r, c))
        all_candidates.append((4, ct_id, bytes(buf)))

        # Delta-row
        buf = bytearray()
        for plane in planes:
            buf.extend(plane[0].tobytes())
            for y in range(1, H):
                delta = ((plane[y].astype(np.int16) - plane[y-1].astype(np.int16)) & 0xFF).astype(np.uint8)
                buf.extend(delta.tobytes())
        all_candidates.append((5, ct_id, bytes(buf)))

        # ---- NEW V4 STRATEGIES ----

        # SPIRAL scan with linear extrapolation (2*prev - prev2)
        spiral_path = _scan_spiral(H, W)
        buf = bytearray()
        for plane in planes:
            prev, prev2 = 0, 0
            for i, (r, c) in enumerate(spiral_path):
                if i == 0: pred = 0
                elif i == 1: pred = prev
                else: pred = max(0, min(255, 2 * prev - prev2))
                buf.append((int(plane[r, c]) - pred) & 0xFF)
                prev2, prev = prev, int(plane[r, c])
        all_candidates.append((6, ct_id, bytes(buf)))

        # SPIRAL scan with avg-neighbor prediction (better on noisy images)
        buf = bytearray()
        for plane in planes:
            known = set()
            for r, c in spiral_path:
                pred = _predict_path_pixel(plane, r, c, known, H, W)
                buf.append((int(plane[r, c]) - pred) & 0xFF)
                known.add((r, c))
        all_candidates.append((13, ct_id, bytes(buf)))

        # CHECKERBOARD scan: phase1 raw sub, phase2 median-predicted
        phase1, phase2 = _scan_checkerboard(H, W)
        buf = bytearray()
        for plane in planes:
            # Phase 1: sub prediction along raster order within phase
            p1_vals = np.array([plane[r,c] for r,c in phase1], dtype=np.uint8)
            p1_filt = VFILTERS_1D[1](p1_vals, len(p1_vals))  # sub filter
            buf.extend(p1_filt.tobytes())
            # Phase 2: predict from 4 cardinal neighbors (all phase-1)
            # First, reconstruct phase-1 into a temporary plane
            tmp = np.zeros((H, W), dtype=np.uint8)
            for i, (r, c) in enumerate(phase1):
                tmp[r, c] = plane[r, c]
            for r, c in phase2:
                pred = _predict_checker_phase2(tmp, r, c, H, W)
                buf.append((int(plane[r, c]) - pred) & 0xFF)
        all_candidates.append((7, ct_id, bytes(buf)))

        # CHECKERBOARD with average predictor (simpler, often better for smooth)
        buf = bytearray()
        for plane in planes:
            p1_vals = np.array([plane[r,c] for r,c in phase1], dtype=np.uint8)
            p1_filt = VFILTERS_1D[1](p1_vals, len(p1_vals))
            buf.extend(p1_filt.tobytes())
            tmp = np.zeros((H, W), dtype=np.uint8)
            for i, (r, c) in enumerate(phase1):
                tmp[r, c] = plane[r, c]
            for r, c in phase2:
                pred = _predict_checker_phase2_avg(tmp, r, c, H, W)
                buf.append((int(plane[r, c]) - pred) & 0xFF)
        all_candidates.append((8, ct_id, bytes(buf)))

        # ZIGZAG scan with sub filter
        zigzag_path = _scan_zigzag(H, W)
        buf = bytearray()
        for plane in planes:
            vals = np.array([plane[r,c] for r,c in zigzag_path], dtype=np.uint8)
            filt = VFILTERS_1D[1](vals, len(vals))  # sub
            buf.extend(filt.tobytes())
        all_candidates.append((9, ct_id, bytes(buf)))

        # ZIGZAG raw (for patterns where zigzag order naturally groups similar values)
        buf = bytearray()
        for plane in planes:
            vals = np.array([plane[r,c] for r,c in zigzag_path], dtype=np.uint8)
            buf.extend(vals.tobytes())
        all_candidates.append((10, ct_id, bytes(buf)))

        # HILBERT scan with sub filter
        hilbert_path = _scan_hilbert(H, W)
        buf = bytearray()
        for plane in planes:
            vals = np.array([plane[r,c] for r,c in hilbert_path], dtype=np.uint8)
            filt = VFILTERS_1D[1](vals, len(vals))
            buf.extend(filt.tobytes())
        all_candidates.append((11, ct_id, bytes(buf)))

        # HILBERT raw
        buf = bytearray()
        for plane in planes:
            vals = np.array([plane[r,c] for r,c in hilbert_path], dtype=np.uint8)
            buf.extend(vals.tobytes())
        all_candidates.append((12, ct_id, bytes(buf)))

    # Phase 2: Fast screen with zlib
    screened = [(sid, ct, pay, len(_fast_compress(pay))) for sid, ct, pay in all_candidates]
    screened.sort(key=lambda x: x[3])

    # Phase 3: Try lzma on top 3
    best_total = 1 << 60
    best_result = None
    for sid, ct, pay, _ in screened[:3]:
        bid, comp = _best_backend(pay)
        if len(comp) < best_total:
            best_total = len(comp)
            best_result = (sid, ct, bid, comp)

    scan_id, ct_id, backend_id, compressed = best_result
    header = struct.pack('<4sBHHBBBB', MAGIC, 4, H, W, channels, scan_id, backend_id, ct_id)
    return header + compressed


# ============================================================
# DECODE
# ============================================================

SCAN_NAMES = {
    0: "raster", 1: "transp", 4: "ring", 5: "delta",
    6: "spiral/ext", 7: "checker/med", 8: "checker/avg",
    9: "zigzag/sub", 10: "zigzag/raw",
    11: "hilbert/sub", 12: "hilbert/raw",
    13: "spiral/avg",
}

def get_scan_name(data):
    if data[:4] != MAGIC: return "?"
    _, _, _, _, _, sid, _, _ = struct.unpack_from('<4sBHHBBBB', data, 0)
    if sid in SCAN_NAMES: return SCAN_NAMES[sid]
    if 0x10 <= sid < 0x10 + N_FILTERS_2D:
        fn = ["none","sub","up","avg","paeth","loco","diag","NE","cross","grad"]
        return f"whole/{fn[sid-0x10]}"
    if 32 <= sid < 36: return f"diag/{['n','s','a','d'][sid-32]}"
    if 48 <= sid < 52: return f"adiag/{['n','s','a','d'][sid-48]}"
    return f"s{sid}"


def decode(data):
    magic = data[:4]
    if magic in (b'TIC1', b'TIC2', b'TIC3'):
        return _dec3(data)

    assert magic == MAGIC, f"Not TIC4 (magic={magic})"
    _, version, H, W, channels, scan_id, backend_id, ct_id = struct.unpack_from('<4sBHHBBBB', data, 0)
    hdr_size = struct.calcsize('<4sBHHBBBB')
    payload = _decompress(data[hdr_size:], backend_id)

    planes = _decode_planes(payload, H, W, channels, scan_id)

    if channels == 1:
        return planes[0]
    img = np.stack(planes, axis=2)
    return CT_INV[ct_id](img)


def _decode_planes(payload, H, W, channels, scan_id):
    # Whole-plane single filter
    if 0x10 <= scan_id < 0x10 + N_FILTERS_2D:
        f_id = scan_id - 0x10
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8); prev = None
            for y in range(H):
                filtered = payload[offset:offset+W]; offset += W
                row = _unfilter_row_2d(filtered, prev, f_id, W)
                plane[y] = np.frombuffer(row, dtype=np.uint8); prev = plane[y]
            planes.append(plane)
        return planes

    # Diagonal/anti-diagonal global filter
    if 32 <= scan_id < 36 or 48 <= scan_id < 52:
        if scan_id >= 48:
            filt_id = scan_id - 48; lines = _scan_antidiag(H, W)
        else:
            filt_id = scan_id - 32; lines = _scan_diagonal(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            for lc in lines:
                L = len(lc)
                if L == 0: continue
                filtered = payload[offset:offset+L]; offset += L
                vals = _unfilter_line_1d(filtered, filt_id, L)
                for i, (r, c) in enumerate(lc):
                    plane[r, c] = vals[i]
            planes.append(plane)
        return planes

    if scan_id == 0:
        # Raster per-row adaptive
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8); prev = None
            for y in range(H):
                f_id = payload[offset]; offset += 1
                filtered = payload[offset:offset+W]; offset += W
                row = _unfilter_row_2d(filtered, prev, f_id, W)
                plane[y] = np.frombuffer(row, dtype=np.uint8); prev = plane[y]
            planes.append(plane)
        return planes

    if scan_id == 1:
        # Transpose per-row adaptive
        planes = []; offset = 0
        for ch in range(channels):
            pt = np.zeros((W, H), dtype=np.uint8); prev = None
            for y in range(W):
                f_id = payload[offset]; offset += 1
                filtered = payload[offset:offset+H]; offset += H
                row = _unfilter_row_2d(filtered, prev, f_id, H)
                pt[y] = np.frombuffer(row, dtype=np.uint8); prev = pt[y]
            planes.append(pt.T)
        return planes

    if scan_id == 4:
        # Ring
        ring_lines = _scan_ring(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8); known = set()
            for ring_pixels in ring_lines:
                for r, c in ring_pixels:
                    pred = _predict_ring_pixel(plane, r, c, known, H, W)
                    plane[r, c] = (int(payload[offset]) + pred) & 0xFF; offset += 1
                for r, c in ring_pixels:
                    known.add((r, c))
            planes.append(plane)
        return planes

    if scan_id == 5:
        # Delta-row
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            plane[0] = np.frombuffer(payload[offset:offset+W], dtype=np.uint8); offset += W
            for y in range(1, H):
                delta = np.frombuffer(payload[offset:offset+W], dtype=np.uint8); offset += W
                plane[y] = ((plane[y-1].astype(np.int16) + delta.astype(np.int16)) & 0xFF).astype(np.uint8)
            planes.append(plane)
        return planes

    if scan_id == 6:
        # Spiral with linear extrapolation
        spiral_path = _scan_spiral(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            prev, prev2 = 0, 0
            for i, (r, c) in enumerate(spiral_path):
                if i == 0: pred = 0
                elif i == 1: pred = prev
                else: pred = max(0, min(255, 2 * prev - prev2))
                plane[r, c] = (int(payload[offset]) + pred) & 0xFF; offset += 1
                prev2, prev = prev, int(plane[r, c])
            planes.append(plane)
        return planes

    if scan_id == 13:
        # Spiral with avg-neighbor prediction
        spiral_path = _scan_spiral(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8); known = set()
            for r, c in spiral_path:
                pred = _predict_path_pixel(plane, r, c, known, H, W)
                plane[r, c] = (int(payload[offset]) + pred) & 0xFF; offset += 1
                known.add((r, c))
            planes.append(plane)
        return planes

    if scan_id in (7, 8):
        # Checkerboard (7=median, 8=avg)
        phase1, phase2 = _scan_checkerboard(H, W)
        pred_fn = _predict_checker_phase2 if scan_id == 7 else _predict_checker_phase2_avg
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            # Phase 1: unfilter sub
            p1_len = len(phase1)
            filtered = payload[offset:offset+p1_len]; offset += p1_len
            vals = _unfilter_line_1d(filtered, 1, p1_len)  # sub filter
            for i, (r, c) in enumerate(phase1):
                plane[r, c] = vals[i]
            # Phase 2: median/avg predicted
            for r, c in phase2:
                pred = pred_fn(plane, r, c, H, W)
                plane[r, c] = (int(payload[offset]) + pred) & 0xFF; offset += 1
            planes.append(plane)
        return planes

    if scan_id in (9, 10):
        # Zigzag (9=sub, 10=raw)
        zigzag_path = _scan_zigzag(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            n_px = len(zigzag_path)
            filtered = payload[offset:offset+n_px]; offset += n_px
            if scan_id == 9:
                vals = _unfilter_line_1d(filtered, 1, n_px)  # sub
            else:
                vals = filtered  # raw
            for i, (r, c) in enumerate(zigzag_path):
                plane[r, c] = vals[i]
            planes.append(plane)
        return planes

    if scan_id in (11, 12):
        # Hilbert (11=sub, 12=raw)
        hilbert_path = _scan_hilbert(H, W)
        planes = []; offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            n_px = len(hilbert_path)
            filtered = payload[offset:offset+n_px]; offset += n_px
            if scan_id == 11:
                vals = _unfilter_line_1d(filtered, 1, n_px)  # sub
            else:
                vals = filtered  # raw
            for i, (r, c) in enumerate(hilbert_path):
                plane[r, c] = vals[i]
            planes.append(plane)
        return planes

    raise ValueError(f"Unknown scan_id {scan_id}")


# ============================================================
# TEST
# ============================================================

def _test():
    np.random.seed(42)
    print("=" * 95)
    print(f"  TIC v{__version__} — Spiral + Checkerboard + Zigzag + Hilbert Scans")
    print("=" * 95)

    passed = 0; failed = 0
    png_wins = 0; png_total = 0

    def check(name, orig, rec, sz, scan=""):
        nonlocal passed, failed
        if np.array_equal(orig, rec):
            passed += 1
            r = orig.size / sz if sz > 0 else 0
            print(f"  OK {name:<35} {orig.size:>8}B -> {sz:>8}B {r:>6.1f}:1 [{scan}]")
        else:
            failed += 1
            print(f"  FAIL {name:<35} {np.sum(orig != rec)} diffs")

    def vs_png(name, img):
        nonlocal png_wins, png_total
        png_total += 1
        d = encode(img)
        sz = len(d)
        ps = png_size(img)
        ratio = ps / sz
        win = "WIN" if ratio > 1.001 else "LOSS"
        if ratio > 1.001: png_wins += 1
        scan = get_scan_name(d)
        print(f"  {name:<30} TIC={sz:>7} PNG={ps:>7} {ratio:.3f}x {win:>4} [{scan}]")

    # Basic roundtrip tests
    print("\n--- Roundtrip Tests ---")
    for N in [64, 128]:
        for name, img in [
            (f"gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))),
            (f"diag_edge_{N}", np.array([[255 if r>c else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"circle_{N}", np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"checker_{N}", np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"noise_{N}", np.random.randint(0,256,(N,N),dtype=np.uint8)),
            (f"solid_{N}", np.full((N,N),128,dtype=np.uint8)),
        ]:
            d = encode(img); rec = decode(d)
            check(name, img, rec, len(d), get_scan_name(d))

    # RGB
    for N in [64, 128]:
        img = np.random.randint(0,256,(N,N,3),dtype=np.uint8)
        d = encode(img); rec = decode(d)
        check(f"rgb_{N}", img, rec, len(d), get_scan_name(d))

    # VS PNG comparison
    print("\n--- TIC v4 vs PNG ---")
    np.random.seed(42)
    for N in [128, 256]:
        vs_png(f"gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1)))
        vs_png(f"diag_grad_{N}", np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8))
        vs_png(f"diag_stripe_{N}", np.array([[(r+c)//8%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8))
        vs_png(f"circle_{N}", np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8))
        vs_png(f"checker_{N}", np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8))
        vs_png(f"noise_{N}", np.random.randint(0,256,(N,N),dtype=np.uint8))
        vs_png(f"smooth_{N}", np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)]),0,255).astype(np.uint8))

    # Real photos
    print("\n--- Real Photos ---")
    from PIL import Image
    for path in ['/home/e/tron.jpeg', '/home/e/t.jpg']:
        if os.path.exists(path):
            img = np.array(Image.open(path).convert('RGB'))
            t0 = time.time()
            d = encode(img)
            t_enc = time.time() - t0
            rec = decode(d)
            ok = np.array_equal(img, rec)
            if ok: passed += 1
            else: failed += 1
            ps = png_size(img)
            sz = len(d)
            ratio = ps / sz
            win = "WIN" if ratio > 1.001 else "LOSS"
            if ratio > 1.001: png_wins += 1; png_total += 1
            scan = get_scan_name(d)
            print(f"  {os.path.basename(path):<25} TIC={sz:>7} PNG={ps:>7} {ratio:.3f}x {win:>4} {'OK' if ok else 'FAIL'} {t_enc:.1f}s [{scan}]")

    print(f"\n  {passed} passed, {failed} failed, {png_wins}/{png_total} beat PNG")
    return failed == 0


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] in ("test", "bench"):
        _test()
    elif sys.argv[1] == "encode" and len(sys.argv) >= 3:
        from PIL import Image
        img = np.array(Image.open(sys.argv[2]).convert('RGB'))
        d = encode(img)
        out = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2].rsplit('.',1)[0] + '.tic'
        with open(out, 'wb') as f: f.write(d)
        print(f"Encoded: {img.size}B -> {len(d)}B ({img.size/len(d):.1f}:1)")
    elif sys.argv[1] == "decode" and len(sys.argv) >= 3:
        with open(sys.argv[2], 'rb') as f: d = f.read()
        img = decode(d)
        out = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2].rsplit('.',1)[0] + '.png'
        from PIL import Image
        Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB').save(out)
        print(f"Decoded: {len(d)}B -> {img.size}B")
