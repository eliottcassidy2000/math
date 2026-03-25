#!/usr/bin/env python3
"""
CHECKERBOARD CODEC: The staircase paradox doesn't apply when you have
ALL 4 cardinal neighbors as context.

Phase 1: Encode EVEN pixels (r+c even) — these form a diamond lattice.
  Treat them as a half-resolution image, compress with MED.
Phase 2: Encode ODD pixels (r+c odd) — each has N,S,E,W from Phase 1.
  Perfect 2D context → predict from 4 neighbors → very low residual entropy.

This dodges the staircase paradox: Phase 2 doesn't scan in any direction.
Every pixel has symmetric context regardless of edge orientation.

Combine with structure-aligned scanning for Phase 1 (where direction matters).
Phase 2 is direction-free by construction.

Also try: for Phase 2, use MED(N,W,NW) which is how PNG would predict
if it could see the south and east neighbors too.

kind-pasteur-2026-03-25-S7
"""
import sys, io, struct, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def compress(data):
    """Compress and return (backend_id, compressed_data). 0=zlib, 1=brotli."""
    results = [(0, zlib.compress(data, 9))]
    if HAS_BROTLI:
        try: results.append((1, brotli.compress(data, quality=11)))
        except: pass
    return min(results, key=lambda x: len(x[1]))

def decomp(bid, data):
    if bid == 0: return zlib.decompress(data)
    return brotli.decompress(data)

def png_size(img):
    pil = Image.fromarray(img, 'L' if img.ndim == 2 else 'RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def med(a, b, c):
    mn, mx = min(a, b), max(a, b)
    if c >= mx: return mn
    if c <= mn: return mx
    return a + b - c

# ============================================================
# CHECKERBOARD CODEC
# ============================================================

def extract_phases(plane):
    """Split plane into even (r+c%2==0) and odd (r+c%2==1) pixels."""
    h, w = plane.shape
    even_coords = [(r,c) for r in range(h) for c in range(w) if (r+c)%2==0]
    odd_coords = [(r,c) for r in range(h) for c in range(w) if (r+c)%2==1]
    return even_coords, odd_coords

def encode_phase1_med(plane, even_coords):
    """Encode even pixels with MED using diagonal neighbors (2-step)."""
    h, w = plane.shape
    known = np.full((h,w), -1, dtype=np.int16)
    residuals = bytearray()

    for r, c in even_coords:
        # Diagonal neighbors at distance sqrt(2) (also even pixels)
        nw = int(known[r-1,c-1]) if r>0 and c>0 and known[r-1,c-1]>=0 else -1
        ne = int(known[r-1,c+1]) if r>0 and c+1<w and known[r-1,c+1]>=0 else -1
        sw = int(known[r+1,c-1]) if r+1<h and c>0 and known[r+1,c-1]>=0 else -1

        # Also check 2-step cardinal (other even pixels)
        up2 = int(known[r-2,c]) if r>=2 and known[r-2,c]>=0 else -1
        left2 = int(known[r,c-2]) if c>=2 and known[r,c-2]>=0 else -1

        # Build prediction from available context
        vals = [v for v in [nw, ne, up2, left2] if v >= 0]
        if len(vals) >= 3:
            # MED-like: use 3 nearest
            a = left2 if left2 >= 0 else (nw if nw >= 0 else vals[0])
            b = up2 if up2 >= 0 else (ne if ne >= 0 else vals[0])
            c2 = nw if nw >= 0 else vals[0]
            pred = med(a, b, c2)
        elif vals:
            pred = sum(vals) // len(vals)
        else:
            pred = 128

        residuals.append((int(plane[r,c]) - pred) & 0xFF)
        known[r,c] = plane[r,c]

    return bytes(residuals), known

def encode_phase2_cardinal(plane, known):
    """Encode odd pixels. Each has N,S,E,W neighbors all known (even pixels).
    Try multiple prediction strategies, pick best per pixel group."""
    h, w = plane.shape
    odd_coords = [(r,c) for r in range(h) for c in range(w) if (r+c)%2==1]

    # Strategy A: simple average of 4 neighbors
    res_avg = bytearray()
    for r, c in odd_coords:
        n = int(known[r-1,c]) if r>0 else -1
        s = int(known[r+1,c]) if r+1<h else -1
        e = int(known[r,c+1]) if c+1<w else -1
        w2 = int(known[r,c-1]) if c>0 else -1
        vals = [v for v in [n,s,e,w2] if v >= 0]
        pred = sum(vals) // len(vals) if vals else 128
        res_avg.append((int(plane[r,c]) - pred) & 0xFF)

    # Strategy B: MED(N,W,NW) — but NW is an odd pixel, not known!
    # Use MED(N,W,avg(E,S)) as a substitute
    res_med = bytearray()
    for r, c in odd_coords:
        n = int(known[r-1,c]) if r>0 else 128
        w2 = int(known[r,c-1]) if c>0 else 128
        s = int(known[r+1,c]) if r+1<h else n
        e = int(known[r,c+1]) if c+1<w else w2
        c2 = (s + e) >> 1  # diagonal estimate
        pred = med(n, w2, c2)
        res_med.append((int(plane[r,c]) - pred) & 0xFF)

    # Strategy C: gradient-adaptive
    res_grad = bytearray()
    for r, c in odd_coords:
        n = int(known[r-1,c]) if r>0 else 128
        s = int(known[r+1,c]) if r+1<h else 128
        e = int(known[r,c+1]) if c+1<w else 128
        w2 = int(known[r,c-1]) if c>0 else 128
        # Gradient detection
        dh = abs(e - w2)  # horizontal variation
        dv = abs(n - s)   # vertical variation
        if dh > 2*dv:     # horizontal edge → predict vertical avg
            pred = (n + s) >> 1
        elif dv > 2*dh:   # vertical edge → predict horizontal avg
            pred = (e + w2) >> 1
        else:
            pred = (n + s + e + w2) >> 2
        res_grad.append((int(plane[r,c]) - pred) & 0xFF)

    return bytes(res_avg), bytes(res_med), bytes(res_grad)

def encode_checker(plane):
    """Full checkerboard encoder: Phase 1 + Phase 2, pick best combo."""
    even_coords, odd_coords = extract_phases(plane)
    phase1_res, known = encode_phase1_med(plane, even_coords)
    phase1_comp = compress(phase1_res)

    p2_avg, p2_med, p2_grad = encode_phase2_cardinal(plane, known)
    p2_options = [compress(p2_avg), compress(p2_med), compress(p2_grad)]
    best_p2_idx = min(range(3), key=lambda i: len(p2_options[i]))
    phase2_comp = p2_options[best_p2_idx]

    total = len(phase1_comp) + len(phase2_comp) + 3  # 3 bytes header
    return total, phase1_comp, phase2_comp, best_p2_idx

def encode_checker_full(plane):
    """Encode plane with checkerboard, return full compressed data with decode info."""
    h, w = plane.shape
    even_coords, odd_coords = extract_phases(plane)
    phase1_res, known = encode_phase1_med(plane, even_coords)
    p1_bid, p1_comp = compress(phase1_res)

    p2_avg, p2_med, p2_grad = encode_phase2_cardinal(plane, known)
    p2_options = [compress(p2_avg), compress(p2_med), compress(p2_grad)]
    best_p2 = min(range(3), key=lambda i: len(p2_options[i][1]))
    p2_bid, p2_comp = p2_options[best_p2]

    # Pack: p1_bid(1) + p1_len(4) + p1_comp + p2_strat(1) + p2_bid(1) + p2_comp
    data = (bytes([p1_bid]) + struct.pack('<I', len(p1_comp)) + p1_comp +
            bytes([best_p2, p2_bid]) + p2_comp)
    return data

def decode_checker(data, h, w):
    """Decode checkerboard-encoded plane."""
    p1_bid = data[0]
    p1_len = struct.unpack_from('<I', data, 1)[0]
    p1_comp = data[5:5+p1_len]
    p2_strat = data[5+p1_len]
    p2_bid = data[6+p1_len]
    p2_comp = data[7+p1_len:]

    p1_raw = decomp(p1_bid, p1_comp)
    p2_raw = decomp(p2_bid, p2_comp)

    plane = np.zeros((h, w), dtype=np.uint8)
    known = np.full((h,w), -1, dtype=np.int16)
    even_coords = [(r,c) for r in range(h) for c in range(w) if (r+c)%2==0]
    odd_coords = [(r,c) for r in range(h) for c in range(w) if (r+c)%2==1]

    # Decode Phase 1
    for i, (r,c) in enumerate(even_coords):
        nw = int(known[r-1,c-1]) if r>0 and c>0 and known[r-1,c-1]>=0 else -1
        ne = int(known[r-1,c+1]) if r>0 and c+1<w and known[r-1,c+1]>=0 else -1
        up2 = int(known[r-2,c]) if r>=2 and known[r-2,c]>=0 else -1
        left2 = int(known[r,c-2]) if c>=2 and known[r,c-2]>=0 else -1

        vals = [v for v in [nw, ne, up2, left2] if v >= 0]
        if len(vals) >= 3:
            a = left2 if left2 >= 0 else (nw if nw >= 0 else vals[0])
            b = up2 if up2 >= 0 else (ne if ne >= 0 else vals[0])
            c2 = nw if nw >= 0 else vals[0]
            pred = med(a, b, c2)
        elif vals:
            pred = sum(vals) // len(vals)
        else:
            pred = 128

        plane[r,c] = (int(p1_raw[i]) + pred) & 0xFF
        known[r,c] = plane[r,c]

    # Decode Phase 2
    for i, (r,c) in enumerate(odd_coords):
        n = int(known[r-1,c]) if r>0 else 128
        s = int(known[r+1,c]) if r+1<h else 128
        e = int(known[r,c+1]) if c+1<w else 128
        w2 = int(known[r,c-1]) if c>0 else 128

        if p2_strat == 0:  # avg
            vals2 = [v for v in [n if r>0 else -1, s if r+1<h else -1,
                                 e if c+1<w else -1, w2 if c>0 else -1] if v >= 0]
            # Recompute properly
            vl = []
            if r>0: vl.append(int(known[r-1,c]))
            if r+1<h: vl.append(int(known[r+1,c]))
            if c+1<w: vl.append(int(known[r,c+1]))
            if c>0: vl.append(int(known[r,c-1]))
            pred = sum(vl)//len(vl) if vl else 128
        elif p2_strat == 1:  # med
            c3 = (s + e) >> 1
            pred = med(n, w2, c3)
        else:  # gradient
            dh = abs(e - w2)
            dv = abs(n - s)
            if dh > 2*dv: pred = (n + s) >> 1
            elif dv > 2*dh: pred = (e + w2) >> 1
            else: pred = (n + s + e + w2) >> 2

        plane[r,c] = (int(p2_raw[i]) + pred) & 0xFF

    return plane

# ============================================================
# ALSO: Standard MED raster for comparison
# ============================================================

def encode_med(plane):
    h, w = plane.shape
    res = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c]) - med(a,b,c2)) & 0xFF)
    _, cdata = compress(bytes(res))
    return cdata

# ============================================================
# BENCHMARK
# ============================================================

print("=" * 80)
print("  CHECKERBOARD CODEC: Direction-free prediction via phase decomposition")
print("  kind-pasteur-2026-03-25-S7")
print("=" * 80)

for SZ in [64, 128, 256]:
    np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2+(y-SZ/2)**2)

    tests = {
        "solid": np.full((SZ,SZ),128,dtype=np.uint8),
        "grad_h": np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1)),
        "grad_d": ((x+y)*255//(2*SZ-2)).astype(np.uint8),
        "stripes_45": ((((x+y)/8).astype(int))%2*255).astype(np.uint8),
        "circles": (np.sin(r/5)*127+128).astype(np.uint8),
        "blob": (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8),
        "random": np.random.randint(0,256,(SZ,SZ),dtype=np.uint8),
        "natural": np.clip(np.array(Image.fromarray(
            np.random.randint(0,256,(max(SZ//16,2),max(SZ//16,2)),dtype=np.uint8)
        ).resize((SZ,SZ),Image.BILINEAR)).astype(float)+np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8),
    }

    print(f"\n--- {SZ}x{SZ} ---")
    print(f"  {'Image':<12} {'PNG':>7} {'MED':>7} {'Checker':>8} {'Raw':>7} {'Best':>7} {'vs PNG':>7} {'RT':>4}")
    print("  " + "-" * 65)

    W, L = 0, 0
    for name, img in sorted(tests.items()):
        ps = png_size(img)
        med_sz = len(encode_med(img)) + 10
        _, raw_comp = compress(img.tobytes())
        raw_sz = len(raw_comp) + 10

        checker_data = encode_checker_full(img)
        checker_sz = len(checker_data) + 10

        # Roundtrip
        dec = decode_checker(checker_data, SZ, SZ)
        ok = np.array_equal(img, dec)

        best = min(med_sz, checker_sz, raw_sz)
        who = "MED" if best == med_sz else ("Checker" if best == checker_sz else "Raw")
        ratio = ps / best if best > 0 else 0
        if ratio > 1.001: W += 1
        elif ratio < 0.999: L += 1

        print(f"  {name:<12} {ps:>7} {med_sz:>7} {checker_sz:>8} {raw_sz:>7} {best:>7} {ratio:>6.2f}x {'OK' if ok else 'FAIL':>4}")

    n = len(tests)
    print(f"\n  {W}W / {n-W-L}T / {L}L out of {n}")
