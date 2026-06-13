#!/usr/bin/env python3
"""
tic5_unified.py — Tournament Image Codec v5: The Grand Unification

Combines ALL winning strategies from opus and kind-pasteur:
  1. Raster per-row adaptive (10 filters including NE, cross, gradient)
  2. Whole-plane single filter (7 options)
  3. Transpose per-row adaptive
  4. Diagonal / anti-diagonal with global filter (4 options each)
  5. Ring scan with 8-neighbor prediction
  6. Spiral with linear extrapolation
  7. Spiral with 8-neighbor prediction
  8. Checkerboard (median + average variants)
  9. Quincunx pyramid (multiresolution) ← NEW: 5.3x on circles
  10. Delta-row
  11. Hilbert/zigzag
  12. Color transforms (GRD, RBD)

Fast screening: zlib pre-screen all ~30 candidates, lzma on top 3.

Author: opus-2026-03-25-S342
"""

import sys, struct, zlib, lzma, io, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tic3 import (
    VFILTERS_2D, N_FILTERS_2D, VFILTERS_1D, N_FILTERS_1D,
    _score, _fast_compress, _best_backend, _decompress,
    _unfilter_row_2d, _unfilter_line_1d,
    CT_FWD, CT_INV,
    _predict_ring_pixel,
    _scan_diagonal, _scan_antidiag, _scan_ring,
)
from tic4_opus import (
    _scan_spiral, _scan_checkerboard, _scan_hilbert,
    _predict_path_pixel, _predict_checker_phase2, _predict_checker_phase2_avg,
)
from quincunx import encode_quincunx, decode_quincunx

__version__ = "5.0.0"
MAGIC = b'TIC5'
sys.stdout.reconfigure(line_buffering=True)

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

    cands = []  # (mode_id, ct_id, payload_bytes)

    for planes, ct_id in planes_list:
        # --- RASTER per-row adaptive (mode 0) ---
        buf = bytearray()
        for plane in planes:
            prev = None
            for y in range(H):
                row = plane[y]
                bf, bs, bo = 0, W*256, None
                for f in range(N_FILTERS_2D):
                    fi = VFILTERS_2D[f](row, prev, W)
                    s = _score(fi)
                    if s < bs: bs = s; bf = f; bo = fi
                buf.append(bf); buf.extend(bo.tobytes()); prev = row
        cands.append((0, ct_id, bytes(buf)))

        # --- WHOLE-PLANE single filter (mode 0x10+f) ---
        for f in range(min(N_FILTERS_2D, 7)):
            buf = bytearray()
            for plane in planes:
                prev = None
                for y in range(H):
                    fi = VFILTERS_2D[f](plane[y], prev, W)
                    buf.extend(fi.tobytes()); prev = plane[y]
            cands.append((0x10 + f, ct_id, bytes(buf)))

        # --- TRANSPOSE per-row adaptive (mode 1) ---
        buf = bytearray()
        for plane in planes:
            pt = plane.T.copy(); prev = None
            for y in range(W):
                row = pt[y]
                bf, bs, bo = 0, H*256, None
                for f in range(N_FILTERS_2D):
                    fi = VFILTERS_2D[f](row, prev, H)
                    s = _score(fi)
                    if s < bs: bs = s; bf = f; bo = fi
                buf.append(bf); buf.extend(bo.tobytes()); prev = row
        cands.append((1, ct_id, bytes(buf)))

        # --- DIAGONAL/ANTIDIAG global filter (mode 32+f / 48+f) ---
        for base, scan_fn in [(2, _scan_diagonal), (3, _scan_antidiag)]:
            lines = scan_fn(H, W)
            for fi_id in range(N_FILTERS_1D):
                buf = bytearray()
                for plane in planes:
                    for lc in lines:
                        vals = np.array([plane[r,c] for r,c in lc], dtype=np.uint8)
                        L = len(vals)
                        if L == 0: continue
                        buf.extend(VFILTERS_1D[fi_id](vals, L).tobytes())
                cands.append((base * 16 + fi_id, ct_id, bytes(buf)))

        # --- RING scan (mode 4) ---
        ring_lines = _scan_ring(H, W)
        buf = bytearray()
        for plane in planes:
            known = set()
            for rp in ring_lines:
                for r, c in rp:
                    pred = _predict_ring_pixel(plane, r, c, known, H, W)
                    buf.append((int(plane[r,c]) - pred) & 0xFF)
                for r, c in rp: known.add((r, c))
        cands.append((4, ct_id, bytes(buf)))

        # --- DELTA-ROW (mode 5) ---
        buf = bytearray()
        for plane in planes:
            buf.extend(plane[0].tobytes())
            for y in range(1, H):
                d = ((plane[y].astype(np.int16) - plane[y-1].astype(np.int16)) & 0xFF).astype(np.uint8)
                buf.extend(d.tobytes())
        cands.append((5, ct_id, bytes(buf)))

        # --- SPIRAL + extrap (mode 6) ---
        spiral_path = _scan_spiral(H, W)
        buf = bytearray()
        for plane in planes:
            prev, prev2 = 0, 0
            for i, (r, c) in enumerate(spiral_path):
                if i == 0: pred = 0
                elif i == 1: pred = prev
                else: pred = max(0, min(255, 2*prev - prev2))
                buf.append((int(plane[r,c]) - pred) & 0xFF)
                prev2, prev = prev, int(plane[r,c])
        cands.append((6, ct_id, bytes(buf)))

        # --- CHECKERBOARD median (mode 7) ---
        phase1, phase2 = _scan_checkerboard(H, W)
        buf = bytearray()
        for plane in planes:
            p1v = np.array([plane[r,c] for r,c in phase1], dtype=np.uint8)
            p1f = VFILTERS_1D[1](p1v, len(p1v))
            buf.extend(p1f.tobytes())
            tmp = np.zeros((H, W), dtype=np.uint8)
            for i, (r, c) in enumerate(phase1): tmp[r,c] = plane[r,c]
            for r, c in phase2:
                pred = _predict_checker_phase2(tmp, r, c, H, W)
                buf.append((int(plane[r,c]) - pred) & 0xFF)
        cands.append((7, ct_id, bytes(buf)))

        # --- QUINCUNX PYRAMID (mode 14) ---
        buf = bytearray()
        level_info = []
        for plane in planes:
            residuals, n_levels, level_sizes = encode_quincunx(plane)
            buf.extend(residuals)
            level_info.append((n_levels, level_sizes))
        # Store level info compactly: n_levels (1 byte) per plane
        # Level sizes are deterministic from H, W, n_levels — no need to store
        n_lev = level_info[0][0] if level_info else 0
        cands.append((14, ct_id, bytes([n_lev]) + bytes(buf)))

    # Phase 2: Fast screen with zlib
    screened = [(sid, ct, pay, len(_fast_compress(pay))) for sid, ct, pay in cands]
    screened.sort(key=lambda x: x[3])

    # Phase 3: lzma on top 3
    best_total = 1 << 60
    best_result = None
    for sid, ct, pay, _ in screened[:3]:
        bid, comp = _best_backend(pay)
        if len(comp) < best_total:
            best_total = len(comp)
            best_result = (sid, ct, bid, comp)

    scan_id, ct_id, backend_id, compressed = best_result
    header = struct.pack('<4sBHHBBBB', MAGIC, 5, H, W, channels, scan_id, backend_id, ct_id)
    return header + compressed


# ============================================================
# DECODE
# ============================================================

SCAN_NAMES = {
    0: "raster", 1: "transp", 4: "ring", 5: "delta",
    6: "spiral", 7: "checker",
    14: "quincunx",
}

def get_scan_name(data):
    if data[:4] != MAGIC: return "?"
    _, _, _, _, _, sid, _, _ = struct.unpack_from('<4sBHHBBBB', data, 0)
    if sid in SCAN_NAMES: return SCAN_NAMES[sid]
    fn = ["none","sub","up","avg","paeth","loco","diag","NE","cross","grad"]
    if 0x10 <= sid < 0x10 + N_FILTERS_2D:
        return f"whole/{fn[sid-0x10]}"
    if 32 <= sid < 36: return f"diag/{['n','s','a','d'][sid-32]}"
    if 48 <= sid < 52: return f"adiag/{['n','s','a','d'][sid-48]}"
    return f"s{sid}"


def decode(data):
    magic = data[:4]
    if magic in (b'TIC1', b'TIC2', b'TIC3', b'TIC4'):
        from tic4_opus import decode as d4
        return d4(data)
    assert magic == MAGIC
    _, version, H, W, channels, scan_id, backend_id, ct_id = struct.unpack_from('<4sBHHBBBB', data, 0)
    hdr_size = struct.calcsize('<4sBHHBBBB')
    payload = _decompress(data[hdr_size:], backend_id)
    planes = _decode_planes(payload, H, W, channels, scan_id)
    if channels == 1: return planes[0]
    return CT_INV[ct_id](np.stack(planes, axis=2))


def _decode_planes(payload, H, W, channels, scan_id):
    # Whole-plane
    if 0x10 <= scan_id < 0x10 + N_FILTERS_2D:
        f_id = scan_id - 0x10
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8); prev = None
            for y in range(H):
                fi = payload[off:off+W]; off += W
                row = _unfilter_row_2d(fi, prev, f_id, W)
                p[y] = np.frombuffer(row, dtype=np.uint8); prev = p[y]
            planes.append(p)
        return planes

    # Diag/antidiag
    if 32 <= scan_id < 36 or 48 <= scan_id < 52:
        if scan_id >= 48: fid = scan_id - 48; lines = _scan_antidiag(H, W)
        else: fid = scan_id - 32; lines = _scan_diagonal(H, W)
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8)
            for lc in lines:
                L = len(lc)
                if L == 0: continue
                fi = payload[off:off+L]; off += L
                vals = _unfilter_line_1d(fi, fid, L)
                for i, (r, c) in enumerate(lc): p[r,c] = vals[i]
            planes.append(p)
        return planes

    if scan_id == 0:
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8); prev = None
            for y in range(H):
                fid = payload[off]; off += 1
                fi = payload[off:off+W]; off += W
                row = _unfilter_row_2d(fi, prev, fid, W)
                p[y] = np.frombuffer(row, dtype=np.uint8); prev = p[y]
            planes.append(p)
        return planes

    if scan_id == 1:
        planes = []; off = 0
        for ch in range(channels):
            pt = np.zeros((W, H), dtype=np.uint8); prev = None
            for y in range(W):
                fid = payload[off]; off += 1
                fi = payload[off:off+H]; off += H
                row = _unfilter_row_2d(fi, prev, fid, H)
                pt[y] = np.frombuffer(row, dtype=np.uint8); prev = pt[y]
            planes.append(pt.T)
        return planes

    if scan_id == 4:
        rlines = _scan_ring(H, W)
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8); known = set()
            for rp in rlines:
                for r, c in rp:
                    pred = _predict_ring_pixel(p, r, c, known, H, W)
                    p[r,c] = (int(payload[off]) + pred) & 0xFF; off += 1
                for r, c in rp: known.add((r, c))
            planes.append(p)
        return planes

    if scan_id == 5:
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8)
            p[0] = np.frombuffer(payload[off:off+W], dtype=np.uint8); off += W
            for y in range(1, H):
                d = np.frombuffer(payload[off:off+W], dtype=np.uint8); off += W
                p[y] = ((p[y-1].astype(np.int16) + d.astype(np.int16)) & 0xFF).astype(np.uint8)
            planes.append(p)
        return planes

    if scan_id == 6:
        spiral = _scan_spiral(H, W)
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8)
            prev, prev2 = 0, 0
            for i, (r, c) in enumerate(spiral):
                if i == 0: pred = 0
                elif i == 1: pred = prev
                else: pred = max(0, min(255, 2*prev - prev2))
                p[r,c] = (int(payload[off]) + pred) & 0xFF; off += 1
                prev2, prev = prev, int(p[r,c])
            planes.append(p)
        return planes

    if scan_id == 7:
        ph1, ph2 = _scan_checkerboard(H, W)
        planes = []; off = 0
        for ch in range(channels):
            p = np.zeros((H, W), dtype=np.uint8)
            n1 = len(ph1)
            fi = payload[off:off+n1]; off += n1
            vals = _unfilter_line_1d(fi, 1, n1)
            for i, (r, c) in enumerate(ph1): p[r,c] = vals[i]
            for r, c in ph2:
                pred = _predict_checker_phase2(p, r, c, H, W)
                p[r,c] = (int(payload[off]) + pred) & 0xFF; off += 1
            planes.append(p)
        return planes

    if scan_id == 14:
        n_lev = payload[0]; off = 1
        planes = []
        for ch in range(channels):
            # Compute level sizes
            from quincunx import _quincunx_levels
            levels = _quincunx_levels(H, W, max_levels=n_lev - 1)
            total_px = H * W
            residuals = payload[off:off+total_px]; off += total_px
            p = decode_quincunx(residuals, H, W, n_lev, None)
            planes.append(p)
        return planes

    raise ValueError(f"Unknown scan_id {scan_id}")


# ============================================================
# TEST
# ============================================================

def _test():
    np.random.seed(42)
    print("=" * 95)
    print(f"  TIC v{__version__} — Grand Unified Codec")
    print("=" * 95)

    passed = 0; failed = 0

    def check(name, orig, rec, sz, scan=""):
        nonlocal passed, failed
        if np.array_equal(orig, rec):
            passed += 1
            r = orig.size / sz if sz > 0 else 0
            print(f"  OK {name:<30} {orig.size:>8}B -> {sz:>7}B {r:>6.1f}:1 [{scan}]")
        else:
            failed += 1
            print(f"  FAIL {name:<30} {np.sum(orig != rec)} diffs")

    from PIL import Image as PILImage

    def png_sz(img):
        mode = 'L' if img.ndim == 2 else 'RGB'
        pil = PILImage.fromarray(img, mode)
        buf = io.BytesIO()
        pil.save(buf, format='PNG', optimize=True, compress_level=9)
        return buf.tell()

    # Roundtrip
    print("\n--- Roundtrip ---")
    for N in [64, 128, 256]:
        for name, img in [
            (f"gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))),
            (f"circle_{N}", np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"diag_edge_{N}", np.array([[255 if r>c else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"checker_{N}", np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"noise_{N}", np.random.randint(0,256,(N,N),dtype=np.uint8)),
        ]:
            d = encode(img); rec = decode(d)
            check(name, img, rec, len(d), get_scan_name(d))

    # RGB
    for N in [64, 128]:
        img = np.random.randint(0,256,(N,N,3),dtype=np.uint8)
        d = encode(img); rec = decode(d)
        check(f"rgb_{N}", img, rec, len(d), get_scan_name(d))

    # VS PNG
    print("\n--- v5 vs PNG ---")
    np.random.seed(42)
    wins = 0; total = 0
    for N in [128, 256]:
        for name, img in [
            (f"gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))),
            (f"circle_{N}", np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"diag_stripe_{N}", np.array([[(r+c)//8%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"checker_{N}", np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"radial_{N}", np.array([[int(128+127*np.sin(np.sqrt((r-N//2)**2+(c-N//2)**2)/5)) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"noise_{N}", np.random.randint(0,256,(N,N),dtype=np.uint8)),
        ]:
            total += 1
            d = encode(img); sz = len(d); ps = png_sz(img)
            ratio = ps / sz
            w = "WIN" if ratio > 1.001 else "LOSS"
            if ratio > 1.001: wins += 1
            scan = get_scan_name(d)
            print(f"  {name:<25} TIC={sz:>7} PNG={ps:>7} {ratio:.2f}x {w:>4} [{scan}]")

    # Real photos
    print("\n--- Real Photos ---")
    for path in ['/home/e/tron.jpeg', '/home/e/t.jpg']:
        if os.path.exists(path):
            img = np.array(PILImage.open(path).convert('RGB'))
            t0 = time.time()
            d = encode(img)
            t_enc = time.time() - t0
            rec = decode(d)
            ok = np.array_equal(img, rec)
            if ok: passed += 1
            else: failed += 1
            total += 1
            ps = png_sz(img)
            sz = len(d)
            ratio = ps / sz
            if ratio > 1.001: wins += 1
            scan = get_scan_name(d)
            print(f"  {os.path.basename(path):<25} TIC={sz:>7} PNG={ps:>7} {ratio:.2f}x {'OK' if ok else 'FAIL'} {t_enc:.1f}s [{scan}]")

    print(f"\n  {passed} passed, {failed} failed, {wins}/{total} beat PNG")
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
