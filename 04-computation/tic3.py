#!/usr/bin/env python3
"""
tic3.py — Tournament Image Codec v3: Structure-Aligned Scanning

LANDMARK INNOVATION: Instead of always scanning left-to-right row-by-row,
detect the dominant structure direction and scan ALONG it.

Key insight from tournament theory: the staircase triangle has THREE natural
directions (horizontal leg, vertical leg, hypotenuse/diagonal). Standard PNG
only uses horizontal. We add diagonal and anti-diagonal scans.

SCAN MODES:
  0: Raster (standard row-by-row, left-to-right)
  1: Transpose (column-by-column, top-to-bottom)
  2: Diagonal scan (along r-c=const lines, SW-to-NE direction)
  3: Anti-diagonal scan (along r+c=const lines, NW-to-SE direction)
  4: Ring scan (center-outward concentric squares)

For each scan mode, per-row adaptive filtering (8 filters) is applied
along the 1D scan lines. The residuals are then compressed.

WHY THIS WORKS:
  Consider a diagonal edge at 45 degrees. In raster mode, the prediction
  crosses the edge every pixel → large residuals. In diagonal mode,
  the prediction follows the edge → near-zero residuals.

  The same logic applies at other angles: choose the scan direction that
  runs PARALLEL to the dominant structure, not perpendicular.

FORMAT: .tic (TIC3)
  Magic: TIC3 (4 bytes)
  Version: 3 (1 byte)
  Height: uint16
  Width: uint16
  Channels: uint8
  ScanMode: uint8
  Backend: uint8
  ColorTransform: uint8
  Compressed payload

Author: opus-2026-03-25-S341
"""

import sys, struct, zlib, lzma, io, time
import numpy as np

__version__ = "3.0.0"
MAGIC = b'TIC3'

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# SCAN ORDER GENERATORS
# ============================================================

def _scan_raster(H, W):
    """Standard row-by-row scan. Returns list of 1D line arrays of (r,c) coords."""
    return [[(r, c) for c in range(W)] for r in range(H)]

def _scan_transpose(H, W):
    """Column-by-column scan."""
    return [[(r, c) for r in range(H)] for c in range(W)]

def _scan_diagonal(H, W):
    """Scan along r-c=const lines (SW-to-NE diagonals).
    Each diagonal: all (r,c) where r-c = k, for k from -(W-1) to H-1.
    Within each diagonal, order by increasing r."""
    lines = []
    for k in range(-(W-1), H):
        line = []
        for r in range(max(0, k), min(H, k + W)):
            c = r - k
            line.append((r, c))
        if line:
            lines.append(line)
    return lines

def _scan_antidiag(H, W):
    """Scan along r+c=const lines (NW-to-SE anti-diagonals).
    Each anti-diagonal: all (r,c) where r+c = k, for k from 0 to H+W-2.
    Within each anti-diagonal, order by increasing r."""
    lines = []
    for k in range(H + W - 1):
        line = []
        for r in range(max(0, k - W + 1), min(H, k + 1)):
            c = k - r
            line.append((r, c))
        if line:
            lines.append(line)
    return lines

def _scan_ring(H, W):
    """Center-outward concentric square rings."""
    cr, cc = H // 2, W // 2
    max_k = max(cr, cc, H - 1 - cr, W - 1 - cc)
    visited = set()
    lines = []
    for k in range(max_k + 1):
        ring = []
        if k == 0:
            if 0 <= cr < H and 0 <= cc < W:
                ring.append((cr, cc))
        else:
            # Top edge
            for c in range(cc - k, cc + k + 1):
                r = cr - k
                if 0 <= r < H and 0 <= c < W and (r, c) not in visited:
                    ring.append((r, c))
            # Bottom edge
            for c in range(cc - k, cc + k + 1):
                r = cr + k
                if 0 <= r < H and 0 <= c < W and (r, c) not in visited:
                    ring.append((r, c))
            # Left edge (excluding corners)
            for r in range(cr - k + 1, cr + k):
                c = cc - k
                if 0 <= r < H and 0 <= c < W and (r, c) not in visited:
                    ring.append((r, c))
            # Right edge (excluding corners)
            for r in range(cr - k + 1, cr + k):
                c = cc + k
                if 0 <= r < H and 0 <= c < W and (r, c) not in visited:
                    ring.append((r, c))
        if ring:
            lines.append(ring)
            for p in ring:
                visited.add(p)
    return lines

SCAN_MODES = [
    ("raster", _scan_raster),
    ("transpose", _scan_transpose),
    ("diagonal", _scan_diagonal),
    ("antidiag", _scan_antidiag),
    ("ring", _scan_ring),
]
N_SCANS = len(SCAN_MODES)

# ============================================================
# VECTORIZED PREDICTION FILTERS (operate on 1D scan lines)
# ============================================================

def _vfilt_none(line_vals, W):
    return line_vals.copy()

def _vfilt_sub(line_vals, W):
    v = line_vals.astype(np.int16)
    pred = np.zeros(W, dtype=np.int16)
    pred[1:] = v[:-1]
    return ((v - pred) & 0xFF).astype(np.uint8)

def _vfilt_avg2(line_vals, W):
    """Average of previous two pixels in the scan line."""
    v = line_vals.astype(np.int16)
    pred = np.zeros(W, dtype=np.int16)
    if W > 1: pred[1] = v[0]
    if W > 2: pred[2:] = (v[:-2] + v[1:-1]) >> 1
    return ((v - pred) & 0xFF).astype(np.uint8)

def _vfilt_diff2(line_vals, W):
    """Linear extrapolation from previous two pixels: 2*prev - prev2."""
    v = line_vals.astype(np.int16)
    pred = np.zeros(W, dtype=np.int16)
    if W > 1: pred[1] = v[0]
    if W > 2: pred[2:] = np.clip(2 * v[1:-1] - v[:-2], 0, 255)
    return ((v - pred) & 0xFF).astype(np.uint8)

VFILTERS_1D = [_vfilt_none, _vfilt_sub, _vfilt_avg2, _vfilt_diff2]
N_FILTERS_1D = len(VFILTERS_1D)
FILTER_NAMES_1D = ["none", "sub", "avg2", "diff2"]

# ============================================================
# 2D PREDICTION FILTERS (for raster/transpose using above-row context)
# ============================================================

def _paeth(a, b, c):
    p = int(a) + int(b) - int(c)
    pa, pb, pc = abs(p - int(a)), abs(p - int(b)), abs(p - int(c))
    if pa <= pb and pa <= pc: return int(a)
    elif pb <= pc: return int(b)
    return int(c)

def _loco(a, b, c):
    a, b, c = int(a), int(b), int(c)
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

def _vfilt2d_none(row, prev, W):
    return row.copy()

def _vfilt2d_sub(row, prev, W):
    r = row.astype(np.int16)
    pred = np.zeros(W, dtype=np.int16); pred[1:] = r[:-1]
    return ((r - pred) & 0xFF).astype(np.uint8)

def _vfilt2d_up(row, prev, W):
    r = row.astype(np.int16)
    p = prev.astype(np.int16) if prev is not None else np.zeros(W, dtype=np.int16)
    return ((r - p) & 0xFF).astype(np.uint8)

def _vfilt2d_avg(row, prev, W):
    r = row.astype(np.int16)
    a = np.zeros(W, dtype=np.int16); a[1:] = r[:-1]
    b = prev.astype(np.int16) if prev is not None else np.zeros(W, dtype=np.int16)
    return ((r - ((a + b) >> 1)) & 0xFF).astype(np.uint8)

def _vfilt2d_paeth(row, prev, W):
    r = row.astype(np.int16)
    a = np.zeros(W, dtype=np.int16); a[1:] = r[:-1]
    if prev is not None:
        b = prev.astype(np.int16)
        c = np.zeros(W, dtype=np.int16); c[1:] = prev[:-1].astype(np.int16)
    else:
        b = np.zeros(W, dtype=np.int16)
        c = np.zeros(W, dtype=np.int16)
    p = a + b - c
    pa = np.abs(p - a); pb = np.abs(p - b); pc = np.abs(p - c)
    pred = np.where((pa <= pb) & (pa <= pc), a, np.where(pb <= pc, b, c))
    return ((r - pred) & 0xFF).astype(np.uint8)

def _vfilt2d_loco(row, prev, W):
    r = row.astype(np.int16)
    a = np.zeros(W, dtype=np.int16); a[1:] = r[:-1]
    if prev is not None:
        b = prev.astype(np.int16)
        c = np.zeros(W, dtype=np.int16); c[1:] = prev[:-1].astype(np.int16)
    else:
        b = np.zeros(W, dtype=np.int16)
        c = np.zeros(W, dtype=np.int16)
    mx = np.maximum(a, b); mn = np.minimum(a, b)
    pred = np.where(c >= mx, mn, np.where(c <= mn, mx, a + b - c))
    return ((r - pred) & 0xFF).astype(np.uint8)

def _vfilt2d_diag(row, prev, W):
    """Upper-left diagonal predictor."""
    r = row.astype(np.int16)
    if prev is not None:
        d = np.zeros(W, dtype=np.int16); d[1:] = prev[:-1].astype(np.int16)
    else:
        d = np.zeros(W, dtype=np.int16)
    return ((r - d) & 0xFF).astype(np.uint8)

def _vfilt2d_ne(row, prev, W):
    """Upper-right diagonal predictor — good for 135-degree edges."""
    r = row.astype(np.int16)
    if prev is not None:
        ne = np.zeros(W, dtype=np.int16)
        ne[:-1] = prev[1:].astype(np.int16)
        ne[-1] = prev[-1].astype(np.int16)
    else:
        ne = np.zeros(W, dtype=np.int16)
    return ((r - ne) & 0xFF).astype(np.uint8)

def _vfilt2d_cross(row, prev, W):
    """Cross-diagonal: average of upper-left and upper-right — good for smooth diagonal gradients."""
    r = row.astype(np.int16)
    if prev is not None:
        nw = np.zeros(W, dtype=np.int16); nw[1:] = prev[:-1].astype(np.int16)
        ne = np.zeros(W, dtype=np.int16); ne[:-1] = prev[1:].astype(np.int16); ne[-1] = prev[-1].astype(np.int16)
        pred = (nw + ne) >> 1
    else:
        pred = np.zeros(W, dtype=np.int16)
    return ((r - pred) & 0xFF).astype(np.uint8)

def _vfilt2d_gradient(row, prev, W):
    """Gradient predictor: a + b - c, clamped. Classic LOCO-I without branch."""
    r = row.astype(np.int16)
    a = np.zeros(W, dtype=np.int16); a[1:] = r[:-1]
    if prev is not None:
        b = prev.astype(np.int16)
        c = np.zeros(W, dtype=np.int16); c[1:] = prev[:-1].astype(np.int16)
    else:
        b = np.zeros(W, dtype=np.int16)
        c = np.zeros(W, dtype=np.int16)
    pred = np.clip(a + b - c, 0, 255)
    return ((r - pred) & 0xFF).astype(np.uint8)

VFILTERS_2D = [_vfilt2d_none, _vfilt2d_sub, _vfilt2d_up, _vfilt2d_avg,
               _vfilt2d_paeth, _vfilt2d_loco, _vfilt2d_diag,
               _vfilt2d_ne, _vfilt2d_cross, _vfilt2d_gradient]
N_FILTERS_2D = len(VFILTERS_2D)

# ============================================================
# RING PREDICTION
# ============================================================

def _predict_ring_pixel(plane, r, c, known_set, H, W):
    """Predict pixel from all known 8-neighbors (average)."""
    total = count = 0
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0: continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < H and 0 <= nc < W and (nr, nc) in known_set:
                total += int(plane[nr, nc]); count += 1
    return total // count if count > 0 else 128

# ============================================================
# SCALAR DECODE FILTERS (sequential, needed for per-row adaptive decode)
# ============================================================

def _scalar_ne(row, prev, x):
    if prev is not None and x + 1 < len(prev): return int(prev[x+1])
    if prev is not None: return int(prev[-1])
    return 0

def _scalar_cross(row, prev, x):
    nw = int(prev[x-1]) if prev is not None and x > 0 else 0
    ne = int(prev[x+1]) if prev is not None and x + 1 < len(prev) else (int(prev[-1]) if prev is not None else 0)
    return (nw + ne) >> 1

def _scalar_gradient(row, prev, x):
    a = int(row[x-1]) if x > 0 else 0
    b = int(prev[x]) if prev is not None else 0
    c = int(prev[x-1]) if prev is not None and x > 0 else 0
    return max(0, min(255, a + b - c))

SCALAR_2D = [
    lambda row, prev, x: 0,
    lambda row, prev, x: int(row[x-1]) if x > 0 else 0,
    lambda row, prev, x: int(prev[x]) if prev is not None else 0,
    lambda row, prev, x: ((int(row[x-1]) if x > 0 else 0) + (int(prev[x]) if prev is not None else 0)) >> 1,
    lambda row, prev, x: _paeth(row[x-1] if x > 0 else 0, prev[x] if prev is not None else 0, prev[x-1] if prev is not None and x > 0 else 0),
    lambda row, prev, x: _loco(row[x-1] if x > 0 else 0, prev[x] if prev is not None else 0, prev[x-1] if prev is not None and x > 0 else 0),
    lambda row, prev, x: int(prev[x-1]) if prev is not None and x > 0 else 0,
    _scalar_ne,
    _scalar_cross,
    _scalar_gradient,
]

def _unfilter_row_2d(filtered, prev, filt_id, W):
    """Reverse 2D filter."""
    pred_fn = SCALAR_2D[filt_id]
    row = bytearray(W)
    for x in range(W):
        pred = pred_fn(row, prev, x)
        row[x] = (filtered[x] + pred) & 0xFF
    return bytes(row)

def _unfilter_line_1d(filtered, filt_id, W):
    """Reverse 1D filter on a scan line."""
    vals = bytearray(W)
    for x in range(W):
        if filt_id == 0:
            vals[x] = filtered[x]
        elif filt_id == 1:
            pred = vals[x-1] if x > 0 else 0
            vals[x] = (filtered[x] + pred) & 0xFF
        elif filt_id == 2:
            pred = vals[x-1] if x == 1 else ((vals[x-1] + vals[x-2]) >> 1 if x >= 2 else 0)
            vals[x] = (filtered[x] + pred) & 0xFF
        elif filt_id == 3:
            if x == 0:
                pred = 0
            elif x == 1:
                pred = vals[0]
            else:
                pred = max(0, min(255, 2 * vals[x-1] - vals[x-2]))
            vals[x] = (filtered[x] + pred) & 0xFF
    return bytes(vals)

# ============================================================
# COLOR SPACE TRANSFORMS
# ============================================================

def _ct_none(img): return img
def _ct_none_inv(img): return img

def _ct_grd(img):
    r, g, b = img[:,:,0].astype(np.int16), img[:,:,1].astype(np.int16), img[:,:,2].astype(np.int16)
    return np.stack([g & 0xFF, (r - g) & 0xFF, (b - g) & 0xFF], axis=-1).astype(np.uint8)

def _ct_grd_inv(img):
    ch = img.astype(np.int16)
    g, rg, bg = ch[:,:,0], ch[:,:,1], ch[:,:,2]
    return np.stack([(rg + g) & 0xFF, g, (bg + g) & 0xFF], axis=-1).astype(np.uint8)

def _ct_rbd(img):
    r, g, b = img[:,:,0].astype(np.int16), img[:,:,1].astype(np.int16), img[:,:,2].astype(np.int16)
    return np.stack([r & 0xFF, (g - r) & 0xFF, (b - r) & 0xFF], axis=-1).astype(np.uint8)

def _ct_rbd_inv(img):
    ch = img.astype(np.int16)
    r, gr, br = ch[:,:,0], ch[:,:,1], ch[:,:,2]
    return np.stack([r, (gr + r) & 0xFF, (br + r) & 0xFF], axis=-1).astype(np.uint8)

CT_FWD = [_ct_none, _ct_grd, _ct_rbd]
CT_INV = [_ct_none_inv, _ct_grd_inv, _ct_rbd_inv]

# ============================================================
# COMPRESSION BACKENDS
# ============================================================

def _fast_compress(data):
    return zlib.compress(data, 9)

def _best_backend(data):
    results = [(0, zlib.compress(data, 9))]
    # Z_FILTERED often better for prediction residuals
    try:
        obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, zlib.Z_FILTERED)
        results.append((0, obj.compress(data) + obj.flush()))
    except: pass
    try: results.append((1, lzma.compress(data, preset=9)))
    except: pass
    return min(results, key=lambda x: len(x[1]))

def _decompress(data, backend_id):
    if backend_id == 0: return zlib.decompress(data)
    if backend_id == 1: return lzma.decompress(data)
    raise ValueError(f"Unknown backend {backend_id}")

# ============================================================
# SCORING
# ============================================================

def _score(filtered):
    f = filtered.astype(np.int16)
    return int(np.sum(np.minimum(f, 256 - f)))

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
        # Pre-screen color transforms
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

    # For each color transform, try each scan mode, pick overall best
    all_candidates = []  # (scan_id, ct_id, payload_bytes)

    for planes, ct_id in planes_list:
        # SCAN 0: Raster with per-row adaptive 2D filtering
        buf = bytearray()
        for plane in planes:
            prev = None
            for y in range(H):
                row = plane[y]
                best_f, best_score, best_out = 0, W * 256, None
                for f in range(N_FILTERS_2D):
                    filt = VFILTERS_2D[f](row, prev, W)
                    s = _score(filt)
                    if s < best_score:
                        best_score = s; best_f = f; best_out = filt
                buf.append(best_f)
                buf.extend(best_out.tobytes())
                prev = row
        all_candidates.append((0, ct_id, bytes(buf)))

        # SCAN 0x10+f: Whole-plane single 2D filter (no per-row IDs — saves H bytes)
        for f in range(min(N_FILTERS_2D, 7)):
            buf = bytearray()
            for plane in planes:
                prev = None
                for y in range(H):
                    filt = VFILTERS_2D[f](plane[y], prev, W)
                    buf.extend(filt.tobytes())
                    prev = plane[y]
            all_candidates.append((0x10 + f, ct_id, bytes(buf)))

        # SCAN 1: Transpose with per-row adaptive 2D filtering
        buf = bytearray()
        for plane in planes:
            pt = plane.T.copy()
            prev = None
            for y in range(W):
                row = pt[y]
                best_f, best_score, best_out = 0, H * 256, None
                for f in range(N_FILTERS_2D):
                    filt = VFILTERS_2D[f](row, prev, H)
                    s = _score(filt)
                    if s < best_score:
                        best_score = s; best_f = f; best_out = filt
                buf.append(best_f)
                buf.extend(best_out.tobytes())
                prev = row
        all_candidates.append((1, ct_id, bytes(buf)))

        # SCAN 2/3: Diagonal/Anti-diagonal with global filter (no per-line overhead)
        # For each direction, try raw + sub filter, pick best
        for scan_id, scan_lines_fn in [(2, _scan_diagonal), (3, _scan_antidiag)]:
            lines = scan_lines_fn(H, W)
            for filt_id in range(N_FILTERS_1D):
                buf = bytearray()
                for plane in planes:
                    for line_coords in lines:
                        line_vals = np.array([plane[r, c] for r, c in line_coords], dtype=np.uint8)
                        L = len(line_vals)
                        if L == 0: continue
                        filt = VFILTERS_1D[filt_id](line_vals, L)
                        buf.extend(filt.tobytes())
                # Encode filter ID into scan_id: scan_id*16 + filt_id
                all_candidates.append((scan_id * 16 + filt_id, ct_id, bytes(buf)))

        # SCAN 4: Ring with average-neighbor prediction
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

        # SCAN 5: Delta-row (fast, good baseline)
        buf = bytearray()
        for plane in planes:
            buf.extend(plane[0].tobytes())
            for y in range(1, H):
                delta = ((plane[y].astype(np.int16) - plane[y-1].astype(np.int16)) & 0xFF).astype(np.uint8)
                buf.extend(delta.tobytes())
        all_candidates.append((5, ct_id, bytes(buf)))

    # Phase 2: Fast screen with zlib
    screened = []
    for scan_id, ct_id, payload in all_candidates:
        zdata = _fast_compress(payload)
        screened.append((scan_id, ct_id, payload, len(zdata)))

    # Phase 3: Try lzma on top 3
    screened.sort(key=lambda x: x[3])
    best_total = 1 << 60
    best_result = None
    for scan_id, ct_id, payload, _ in screened[:3]:
        bid, compressed = _best_backend(payload)
        if len(compressed) < best_total:
            best_total = len(compressed)
            best_result = (scan_id, ct_id, bid, compressed)

    scan_id, ct_id, backend_id, compressed = best_result
    header = struct.pack('<4sBHHBBBB', MAGIC, 3, H, W, channels, scan_id, backend_id, ct_id)
    return header + compressed


# ============================================================
# DECODE
# ============================================================

def decode(data):
    magic = data[:4]
    if magic == b'TIC2':
        # Backward compat — delegate to tic2
        import tic2
        return tic2.decode(data)

    assert magic == MAGIC, f"Not TIC3 (magic={magic})"
    _, version, H, W, channels, scan_id, backend_id, ct_id = struct.unpack_from('<4sBHHBBBB', data, 0)
    hdr_size = struct.calcsize('<4sBHHBBBB')
    payload = _decompress(data[hdr_size:], backend_id)

    planes = _decode_planes(payload, H, W, channels, scan_id)

    if channels == 1:
        return planes[0]

    img = np.stack(planes, axis=2)
    return CT_INV[ct_id](img)


def _decode_planes(payload, H, W, channels, scan_id):
    # Whole-plane single filter: scan_id = 0x10 + filter_id
    if 0x10 <= scan_id < 0x10 + N_FILTERS_2D:
        f_id = scan_id - 0x10
        planes = []
        offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            prev = None
            for y in range(H):
                filtered = payload[offset:offset + W]; offset += W
                row = _unfilter_row_2d(filtered, prev, f_id, W)
                plane[y] = np.frombuffer(row, dtype=np.uint8)
                prev = plane[y]
            planes.append(plane)
        return planes

    if scan_id == 0:
        # Raster with per-row adaptive 2D
        planes = []
        offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            prev = None
            for y in range(H):
                f_id = payload[offset]; offset += 1
                filtered = payload[offset:offset + W]; offset += W
                row = _unfilter_row_2d(filtered, prev, f_id, W)
                plane[y] = np.frombuffer(row, dtype=np.uint8)
                prev = plane[y]
            planes.append(plane)
        return planes

    if scan_id == 1:
        # Transpose with per-row adaptive 2D
        planes = []
        offset = 0
        for ch in range(channels):
            plane_t = np.zeros((W, H), dtype=np.uint8)
            prev = None
            for y in range(W):
                f_id = payload[offset]; offset += 1
                filtered = payload[offset:offset + H]; offset += H
                row = _unfilter_row_2d(filtered, prev, f_id, H)
                plane_t[y] = np.frombuffer(row, dtype=np.uint8)
                prev = plane_t[y]
            planes.append(plane_t.T)
        return planes

    # Diagonal/Anti-diagonal with global filter: scan_id = direction*16 + filt_id
    if 32 <= scan_id < 32 + N_FILTERS_1D or 48 <= scan_id < 48 + N_FILTERS_1D:
        if scan_id >= 48:
            direction = 3  # anti-diagonal
            filt_id = scan_id - 48
            lines_fn = _scan_antidiag
        else:
            direction = 2  # diagonal
            filt_id = scan_id - 32
            lines_fn = _scan_diagonal
        lines = lines_fn(H, W)
        planes = []
        offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            for line_coords in lines:
                L = len(line_coords)
                if L == 0: continue
                filtered = payload[offset:offset + L]; offset += L
                vals = _unfilter_line_1d(filtered, filt_id, L)
                for i, (r, c) in enumerate(line_coords):
                    plane[r, c] = vals[i]
            planes.append(plane)
        return planes

    if scan_id == 4:
        # Ring with avg-neighbor prediction
        ring_lines = _scan_ring(H, W)
        planes = []
        offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            known = set()
            for ring_pixels in ring_lines:
                for r, c in ring_pixels:
                    pred = _predict_ring_pixel(plane, r, c, known, H, W)
                    plane[r, c] = (int(payload[offset]) + pred) & 0xFF
                    offset += 1
                for r, c in ring_pixels:
                    known.add((r, c))
            planes.append(plane)
        return planes

    if scan_id == 5:
        # Delta-row
        planes = []
        offset = 0
        for ch in range(channels):
            plane = np.zeros((H, W), dtype=np.uint8)
            plane[0] = np.frombuffer(payload[offset:offset + W], dtype=np.uint8)
            offset += W
            for y in range(1, H):
                delta = np.frombuffer(payload[offset:offset + W], dtype=np.uint8)
                offset += W
                plane[y] = ((plane[y-1].astype(np.int16) + delta.astype(np.int16)) & 0xFF).astype(np.uint8)
            planes.append(plane)
        return planes

    raise ValueError(f"Unknown scan_id {scan_id}")


# ============================================================
# PNG COMPARISON
# ============================================================

def png_size(img):
    from PIL import Image
    mode = 'L' if img.ndim == 2 else 'RGB'
    pil = Image.fromarray(img, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()


# ============================================================
# TEST SUITE
# ============================================================

def _test():
    np.random.seed(42)
    print("=" * 90)
    print(f"  TIC v{__version__} — Structure-Aligned Scanning")
    print("=" * 90)

    passed = 0; failed = 0

    def check(name, original, recovered, tic_size, scan_name=""):
        nonlocal passed, failed
        if np.array_equal(original, recovered):
            passed += 1
            raw = original.size
            ratio = raw / tic_size if tic_size > 0 else 0
            print(f"  OK {name:<35} {raw:>8}B -> {tic_size:>8}B  {ratio:>6.1f}:1  [{scan_name}]")
        else:
            failed += 1
            diff = np.sum(original != recovered)
            print(f"  FAIL {name:<35} {diff} pixels differ!")

    filt2d_names = ["none","sub","up","avg","paeth","loco","diag","NE","cross","grad"]
    def get_scan_name(data):
        _, _, _, _, _, scan_id, _, _ = struct.unpack_from('<4sBHHBBBB', data, 0)
        if scan_id < len(SCAN_MODES):
            return SCAN_MODES[scan_id][0]
        if 0x10 <= scan_id < 0x10 + N_FILTERS_2D:
            return f"whole/{filt2d_names[scan_id-0x10]}"
        if 32 <= scan_id < 36:
            return f"diag/{FILTER_NAMES_1D[scan_id-32]}"
        if 48 <= scan_id < 52:
            return f"adiag/{FILTER_NAMES_1D[scan_id-48]}"
        return f"scan{scan_id}"

    # ---- Grayscale tests ----
    print("\n--- Grayscale ---")
    for N in [64, 128]:
        for name, img in [
            (f"gradient_{N}", np.tile(np.linspace(0, 255, N, dtype=np.uint8), (N, 1))),
            (f"diagonal_grad_{N}", np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"diag_edge_{N}", np.array([[255 if r > c else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"antidiag_edge_{N}", np.array([[255 if r + c > N else 0 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"circle_{N}", np.array([[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"noise_{N}", np.random.randint(0, 256, (N, N), dtype=np.uint8)),
            (f"solid_{N}", np.full((N, N), 128, dtype=np.uint8)),
        ]:
            compressed = encode(img)
            recovered = decode(compressed)
            check(name, img, recovered, len(compressed), get_scan_name(compressed))

    # ---- Diagonal structure tests (the KEY test) ----
    print("\n--- Diagonal Structure (where this codec should excel) ---")
    for N in [128, 256]:
        # 45-degree stripes
        img = np.array([[(r + c) // 8 % 2 * 255 for c in range(N)] for r in range(N)], dtype=np.uint8)
        compressed = encode(img)
        recovered = decode(compressed)
        check(f"diag_stripes_{N}", img, recovered, len(compressed), get_scan_name(compressed))

        # 135-degree stripes (anti-diagonal)
        img = np.array([[(r - c) // 8 % 2 * 255 for c in range(N)] for r in range(N)], dtype=np.uint8)
        compressed = encode(img)
        recovered = decode(compressed)
        check(f"antidiag_stripes_{N}", img, recovered, len(compressed), get_scan_name(compressed))

        # Diagonal sine wave
        img = np.array([[int(128 + 127 * np.sin((r + c) / 10)) for c in range(N)] for r in range(N)], dtype=np.uint8)
        compressed = encode(img)
        recovered = decode(compressed)
        check(f"diag_sine_{N}", img, recovered, len(compressed), get_scan_name(compressed))

    # ---- RGB ----
    print("\n--- RGB ---")
    for N in [64, 128]:
        img = np.random.randint(0, 256, (N, N, 3), dtype=np.uint8)
        compressed = encode(img)
        recovered = decode(compressed)
        check(f"rgb_noise_{N}", img, recovered, len(compressed), get_scan_name(compressed))

    # ---- vs PNG ----
    print("\n--- TIC v3 vs PNG ---")
    comparisons = [
        ("gradient_128", np.tile(np.linspace(0, 255, 128, dtype=np.uint8), (128, 1))),
        ("diag_grad_128", np.array([[int(255*(r+c)/254) for c in range(128)] for r in range(128)], dtype=np.uint8)),
        ("diag_stripe_128", np.array([[(r+c)//8%2*255 for c in range(128)] for r in range(128)], dtype=np.uint8)),
        ("circle_128", np.array([[min(255, int(255*np.exp(-((r-64)**2+(c-64)**2)/2048))) for c in range(128)] for r in range(128)], dtype=np.uint8)),
        ("noise_128", np.random.randint(0, 256, (128, 128), dtype=np.uint8)),
    ]
    wins = 0
    for name, img in comparisons:
        tic_data = encode(img)
        tic_sz = len(tic_data)
        png_sz = png_size(img)
        ratio = png_sz / tic_sz
        win = "WIN" if ratio > 1.001 else ("TIE" if ratio > 0.999 else "LOSS")
        if ratio > 1.001: wins += 1
        scan = get_scan_name(tic_data)
        print(f"  {name:<25} TIC={tic_sz:>7} PNG={png_sz:>7} {ratio:.3f}x {win:>4} [{scan}]")

    # ---- Real photos ----
    print("\n--- Real Photos ---")
    import os
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg', '/home/e/t.jpg']:
        if os.path.exists(path):
            img = np.array(PILImage.open(path).convert('RGB'))
            t0 = time.time()
            tic_data = encode(img)
            t_enc = time.time() - t0
            recovered = decode(tic_data)
            ok = np.array_equal(img, recovered)
            if ok: passed += 1
            else: failed += 1
            tic_sz = len(tic_data)
            png_sz = png_size(img)
            ratio = png_sz / tic_sz
            win = "WIN" if ratio > 1.001 else "LOSS"
            if ratio > 1.001: wins += 1
            scan = get_scan_name(tic_data)
            print(f"  {os.path.basename(path):<25} TIC={tic_sz:>7} PNG={png_sz:>7} {ratio:.3f}x {win:>4} {'OK' if ok else 'FAIL'} {t_enc:.1f}s [{scan}]")

    # ---- Speed ----
    print("\n--- Speed ---")
    test_img = np.random.randint(0, 256, (128, 128), dtype=np.uint8)
    t0 = time.time()
    for _ in range(3):
        encode(test_img)
    t = (time.time() - t0) / 3
    print(f"  128x128 gray: {t*1000:.0f}ms ({1/t:.0f} img/s)")

    print(f"\n  RESULTS: {passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] in ("test", "bench"):
        _test()
    elif sys.argv[1] == "encode" and len(sys.argv) >= 3:
        from PIL import Image
        img = np.array(Image.open(sys.argv[2]).convert('RGB'))
        compressed = encode(img)
        outfile = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2].rsplit('.', 1)[0] + '.tic'
        with open(outfile, 'wb') as f:
            f.write(compressed)
        print(f"Encoded: {img.size}B -> {len(compressed)}B ({img.size/len(compressed):.1f}:1)")
    elif sys.argv[1] == "decode" and len(sys.argv) >= 3:
        with open(sys.argv[2], 'rb') as f:
            compressed = f.read()
        image = decode(compressed)
        outfile = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2].rsplit('.', 1)[0] + '.png'
        from PIL import Image
        if image.ndim == 2:
            Image.fromarray(image, 'L').save(outfile)
        else:
            Image.fromarray(image, 'RGB').save(outfile)
        print(f"Decoded: {len(compressed)}B -> {image.size}B")
    else:
        print(f"TIC v{__version__}")
        print("  python3 tic3.py test")
        print("  python3 tic3.py encode image.png out.tic")
        print("  python3 tic3.py decode out.tic image.png")
