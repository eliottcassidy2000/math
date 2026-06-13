#!/usr/bin/env python3
"""
wavefront.py — Adaptive Wavefront Codec

THE AXIOM THIS BREAKS: "The scan order must be a predetermined 1D sequence."

Instead: pixels are processed in order of EXPANDABLE CONFIDENCE.
The "wavefront" is the boundary between decoded and undecoded regions.
At each step, the wavefront expands to the pixel that can be predicted
with HIGHEST CONFIDENCE from the already-decoded region.

Confidence = inverse of LOCAL VARIANCE among known neighbors.
Low variance = smooth region = confident prediction = process first.
High variance = edge/texture = uncertain prediction = defer.

THE KEY INSIGHT: this is NOT a tree walk (MST failed because tree has overhead).
This is a BFS-LIKE expansion where the wavefront is implicit — it's the
set of undecoded pixels that have at least one decoded neighbor.
The expansion order is deterministic from the decoded values alone.

STRUCTURE: It's a GRAPH (not a tree). Each pixel can be predicted from
MULTIPLE decoded neighbors (weighted by distance and confidence).
The prediction graph has cycles — pixel A's prediction influenced B,
and B's prediction influenced C, but C and A might be neighbors.

Second idea: ADAPTIVE LATTICE WAVELET
Instead of fixed wavelet basis, use the checkerboard lattice but with
ADAPTIVE lifting: the lifting coefficients depend on local content.

Author: opus-2026-03-25-S343
"""

import sys, struct, zlib, lzma, io, time, os
import numpy as np
import heapq

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# ADAPTIVE WAVEFRONT SCAN
# ============================================================

def _wavefront_encode(plane):
    """Encode using adaptive wavefront expansion.

    Process pixels in order of prediction confidence.
    Confidence = 1 / (local_variance + 1) among known 8-neighbors.
    Higher confidence → process first → smaller residuals.
    """
    H, W = plane.shape
    decoded = np.full((H, W), -1, dtype=np.int16)  # -1 = unknown
    in_queue = np.zeros((H, W), dtype=bool)
    residuals = bytearray()
    order = []  # NOT stored — recomputed at decode

    # Priority queue: (-confidence, tiebreaker, r, c)
    heap = []
    tiebreaker = 0

    # Seed: pixel with minimum local variance (smoothest region)
    # For speed, just use center pixel
    cr, cc = H // 2, W // 2
    heapq.heappush(heap, (0, tiebreaker, cr, cc))
    in_queue[cr, cc] = True
    tiebreaker += 1

    while heap:
        neg_conf, _, r, c = heapq.heappop(heap)

        if decoded[r, c] >= 0:
            continue  # already processed

        # Predict from known 8-neighbors
        total = 0.0
        weight = 0.0
        values = []
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and decoded[nr, nc] >= 0:
                    v = float(decoded[nr, nc])
                    w = 1.0 if abs(dr) + abs(dc) == 1 else 0.707
                    total += w * v
                    weight += w
                    values.append(v)

        if weight > 0:
            pred = int(round(total / weight))
        else:
            pred = 128

        residual = (int(plane[r, c]) - pred) & 0xFF
        residuals.append(residual)
        decoded[r, c] = int(plane[r, c])
        order.append((r, c))

        # Add unprocessed neighbors to queue with confidence based on
        # THEIR known neighbors' variance
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and decoded[nr, nc] < 0 and not in_queue[nr, nc]:
                    # Compute confidence for (nr, nc)
                    nvals = []
                    for dr2 in [-1, 0, 1]:
                        for dc2 in [-1, 0, 1]:
                            if dr2 == 0 and dc2 == 0: continue
                            nnr, nnc = nr + dr2, nc + dc2
                            if 0 <= nnr < H and 0 <= nnc < W and decoded[nnr, nnc] >= 0:
                                nvals.append(float(decoded[nnr, nnc]))

                    if nvals:
                        # Confidence = -variance (we want LOW variance first)
                        # Plus count of known neighbors (more = better)
                        variance = np.var(nvals) if len(nvals) > 1 else 0
                        n_known = len(nvals)
                        # Priority: lower is processed first
                        # Prefer: many known neighbors + low variance
                        priority = variance - n_known * 100
                    else:
                        priority = 1e9  # no known neighbors = very low confidence

                    heapq.heappush(heap, (priority, tiebreaker, nr, nc))
                    in_queue[nr, nc] = True
                    tiebreaker += 1

    return bytes(residuals)


def _wavefront_decode(data, H, W):
    """Decode wavefront — MUST reproduce exact same expansion order."""
    decoded = np.full((H, W), -1, dtype=np.int16)
    in_queue = np.zeros((H, W), dtype=bool)
    heap = []
    tiebreaker = 0

    cr, cc = H // 2, W // 2
    heapq.heappush(heap, (0, tiebreaker, cr, cc))
    in_queue[cr, cc] = True
    tiebreaker += 1

    idx = 0
    while heap:
        neg_conf, _, r, c = heapq.heappop(heap)
        if decoded[r, c] >= 0:
            continue

        # Same prediction as encoder
        total = 0.0
        weight = 0.0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and decoded[nr, nc] >= 0:
                    w = 1.0 if abs(dr) + abs(dc) == 1 else 0.707
                    total += w * float(decoded[nr, nc])
                    weight += w

        pred = int(round(total / weight)) if weight > 0 else 128
        decoded[r, c] = (int(data[idx]) + pred) & 0xFF
        idx += 1

        # Same expansion as encoder
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and decoded[nr, nc] < 0 and not in_queue[nr, nc]:
                    nvals = []
                    for dr2 in [-1, 0, 1]:
                        for dc2 in [-1, 0, 1]:
                            if dr2 == 0 and dc2 == 0: continue
                            nnr, nnc = nr + dr2, nc + dc2
                            if 0 <= nnr < H and 0 <= nnc < W and decoded[nnr, nnc] >= 0:
                                nvals.append(float(decoded[nnr, nnc]))
                    if nvals:
                        variance = np.var(nvals) if len(nvals) > 1 else 0
                        n_known = len(nvals)
                        priority = variance - n_known * 100
                    else:
                        priority = 1e9
                    heapq.heappush(heap, (priority, tiebreaker, nr, nc))
                    in_queue[nr, nc] = True
                    tiebreaker += 1

    return decoded.astype(np.uint8)


# ============================================================
# ADAPTIVE LIFTING ON CHECKERBOARD LATTICE
# ============================================================

def _adaptive_lift_encode(plane):
    """Checkerboard lifting with adaptive coefficients.

    Phase 1 (white squares): encode with sub filter along diag lines
    Phase 2 (black squares): predict from 4 cardinal (phase-1) neighbors
    using ADAPTIVE weights based on local gradient.

    If horizontal gradient > vertical: weight left/right more.
    If vertical gradient > horizontal: weight up/down more.
    """
    H, W = plane.shape
    buf = bytearray()

    # Phase 1: white squares (r+c even) in raster order
    white = []
    for r in range(H):
        for c in range(W):
            if (r + c) % 2 == 0:
                white.append((r, c))

    # Encode phase 1 with simple sub filter
    prev = 0
    for r, c in white:
        v = int(plane[r, c])
        buf.append((v - prev) & 0xFF)
        prev = v

    # Phase 2: black squares, adaptive prediction from 4 cardinal neighbors
    # Store phase-1 values in temp plane
    tmp = np.zeros((H, W), dtype=np.int16)
    for r, c in white:
        tmp[r, c] = plane[r, c]

    for r in range(H):
        for c in range(W):
            if (r + c) % 2 == 1:
                # Get cardinal neighbors (all phase-1)
                n = int(tmp[r-1, c]) if r > 0 else -1
                s = int(tmp[r+1, c]) if r < H-1 else -1
                w = int(tmp[r, c-1]) if c > 0 else -1
                e = int(tmp[r, c+1]) if c < W-1 else -1

                vals = []
                if n >= 0: vals.append(('v', n))
                if s >= 0: vals.append(('v', s))
                if w >= 0: vals.append(('h', w))
                if e >= 0: vals.append(('h', e))

                if len(vals) >= 2:
                    # Compute directional gradients
                    h_vals = [v for d, v in vals if d == 'h']
                    v_vals = [v for d, v in vals if d == 'v']
                    h_range = max(h_vals) - min(h_vals) if len(h_vals) >= 2 else 0
                    v_range = max(v_vals) - min(v_vals) if len(v_vals) >= 2 else 0

                    if h_range > v_range * 2:
                        # Strong horizontal gradient → predict from vertical neighbors
                        pred = sum(v for d, v in vals if d == 'v') // max(1, len(v_vals))
                    elif v_range > h_range * 2:
                        # Strong vertical gradient → predict from horizontal neighbors
                        pred = sum(v for d, v in vals if d == 'h') // max(1, len(h_vals))
                    else:
                        # Mixed → simple average
                        pred = sum(v for _, v in vals) // len(vals)
                else:
                    pred = vals[0][1] if vals else 128

                buf.append((int(plane[r, c]) - pred) & 0xFF)

    return bytes(buf)


def _adaptive_lift_decode(data, H, W):
    """Decode adaptive lifting."""
    plane = np.zeros((H, W), dtype=np.uint8)
    offset = 0

    # Phase 1: white squares
    white = []
    for r in range(H):
        for c in range(W):
            if (r + c) % 2 == 0:
                white.append((r, c))

    prev = 0
    for r, c in white:
        plane[r, c] = (int(data[offset]) + prev) & 0xFF
        prev = int(plane[r, c])
        offset += 1

    # Phase 2: black squares with adaptive prediction
    tmp = plane.astype(np.int16)
    for r in range(H):
        for c in range(W):
            if (r + c) % 2 == 1:
                n = int(tmp[r-1, c]) if r > 0 else -1
                s = int(tmp[r+1, c]) if r < H-1 else -1
                w = int(tmp[r, c-1]) if c > 0 else -1
                e = int(tmp[r, c+1]) if c < W-1 else -1

                vals = []
                if n >= 0: vals.append(('v', n))
                if s >= 0: vals.append(('v', s))
                if w >= 0: vals.append(('h', w))
                if e >= 0: vals.append(('h', e))

                if len(vals) >= 2:
                    h_vals = [v for d, v in vals if d == 'h']
                    v_vals = [v for d, v in vals if d == 'v']
                    h_range = max(h_vals) - min(h_vals) if len(h_vals) >= 2 else 0
                    v_range = max(v_vals) - min(v_vals) if len(v_vals) >= 2 else 0

                    if h_range > v_range * 2:
                        pred = sum(v for d, v in vals if d == 'v') // max(1, len(v_vals))
                    elif v_range > h_range * 2:
                        pred = sum(v for d, v in vals if d == 'h') // max(1, len(h_vals))
                    else:
                        pred = sum(v for _, v in vals) // len(vals)
                else:
                    pred = vals[0][1] if vals else 128

                plane[r, c] = (int(data[offset]) + pred) & 0xFF
                offset += 1

    return plane


# ============================================================
# TEST
# ============================================================

def png_size(img):
    from PIL import Image
    pil = Image.fromarray(img, 'L')
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def _loco(a, b, c):
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

def _encode_med(plane):
    H, W = plane.shape; buf = bytearray()
    for r in range(H):
        for c in range(W):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            buf.append((int(plane[r,c]) - _loco(a,b,d)) & 0xFF)
    return bytes(buf)


if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 85)
    print("  Wavefront + Adaptive Lifting Codec")
    print("  opus-2026-03-25-S343")
    print("=" * 85)

    tests = {}
    for N in [64, 128]:
        tests[f"circle_{N}"] = np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"gradient_{N}"] = np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))
        tests[f"diag_edge_{N}"] = np.array([[255 if r > c else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"checker_{N}"] = np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f"noise_{N}"] = np.random.randint(0, 256, (N, N), dtype=np.uint8)
        tests[f"smooth_{N}"] = np.clip(np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)]),0,255).astype(np.uint8)
        tests[f"blob_{N}"] = np.array([[min(255,int(200*np.exp(-((r-N//3)**2+(c-2*N//3)**2)/(N**2/16))+100*np.exp(-((r-2*N//3)**2+(c-N//3)**2)/(N**2/20)))) for c in range(N)] for r in range(N)], dtype=np.uint8)

    print(f"\n  {'Image':<25} {'PNG':>6} {'MED':>6} {'WF':>6} {'ALift':>6} {'Best':>8} {'OK':>4}")
    print("  " + "-" * 70)

    for name, img in sorted(tests.items()):
        H, W = img.shape
        ps = png_size(img)

        # MED
        med_data = _encode_med(img)
        c_med = len(zlib.compress(med_data, 9))

        # Wavefront
        t0 = time.time()
        wf_data = _wavefront_encode(img)
        t_wf = time.time() - t0
        wf_rec = _wavefront_decode(wf_data, H, W)
        wf_ok = np.array_equal(img, wf_rec)
        c_wf = len(zlib.compress(wf_data, 9))

        # Adaptive lifting
        al_data = _adaptive_lift_encode(img)
        al_rec = _adaptive_lift_decode(al_data, H, W)
        al_ok = np.array_equal(img, al_rec)
        c_al = len(zlib.compress(al_data, 9))

        best_name = "MED"
        best_val = c_med
        if c_wf < best_val: best_name = "WF"; best_val = c_wf
        if c_al < best_val: best_name = "ALift"; best_val = c_al

        ok = "OK" if (wf_ok and al_ok) else f"{'WF' if not wf_ok else ''}{'AL' if not al_ok else ''}"
        print(f"  {name:<25} {ps:>6} {c_med:>6} {c_wf:>6} {c_al:>6} {best_name:>5} {ps/best_val:.2f}x {ok:>4}")

    # Real photo
    print("\n--- Real Photo ---")
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg']:
        if os.path.exists(path):
            img = np.array(PILImage.open(path).convert('L'))[:64,:64]
            H, W = img.shape
            ps = png_size(img)
            c_med = len(zlib.compress(_encode_med(img), 9))
            t0 = time.time()
            wf_data = _wavefront_encode(img)
            t_wf = time.time() - t0
            wf_rec = _wavefront_decode(wf_data, H, W)
            c_wf = len(zlib.compress(wf_data, 9))
            al_data = _adaptive_lift_encode(img)
            al_rec = _adaptive_lift_decode(al_data, H, W)
            c_al = len(zlib.compress(al_data, 9))
            wf_ok = np.array_equal(img, wf_rec)
            al_ok = np.array_equal(img, al_rec)
            print(f"  tron_64x64  PNG={ps} MED={c_med} WF={c_wf} ALift={c_al} WF_ok={wf_ok} AL_ok={al_ok} ({t_wf*1000:.0f}ms)")
