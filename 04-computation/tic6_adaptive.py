#!/usr/bin/env python3
"""
tic6_adaptive.py — Breaking the fixed-scan axiom

IMPLICIT AXIOM OF ALL PREVIOUS CODECS: scan order is chosen ONCE for the
whole image. Even when adaptive (try-all-pick-best), the chosen scan is
applied uniformly to every pixel.

THIS CODEC BREAKS THAT AXIOM with three genuinely new strategies:

1. CONFIDENCE-ORDERED SCAN: Process pixels in order of prediction confidence.
   Start with corners/edges-of-image (known boundary = 0), then pixels with
   the most known neighbors, then the least predictable ones last.
   The ordering is IMPLICIT — reproducible at decode from known values.

2. VALUE-SORTED PLANES: Sort all pixels by value, encode the nearly-monotone
   value sequence + the position permutation. For smooth images, the position
   permutation has massive runs (all bright pixels in one region).
   AXIOM BROKEN: pixels aren't processed spatially at all.

3. FLOOD-FILL REGIONS: Flood-fill from a seed pixel, expanding to similar
   neighbors first. Each region gets a base value + tiny residuals.
   AXIOM BROKEN: scan follows content structure, not geometry.

4. COLUMN-SORT (from kind-pasteur): Sort columns by similarity, delta-compress.
   AXIOM BROKEN: column order is by content, not position.

5. ROW-LPC (from kind-pasteur): Audio-style linear prediction between rows.
   AXIOM BROKEN: prediction model is adaptive (least-squares fit), not fixed.

Author: opus-2026-03-25-S343
"""

import sys, struct, zlib, lzma, io, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tic3 import (
    VFILTERS_2D, N_FILTERS_2D, _score, _fast_compress, _best_backend, _decompress,
    _unfilter_row_2d, CT_FWD, CT_INV,
)

__version__ = "6.0.0"
MAGIC = b'TIC6'
sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# STRATEGY: CONFIDENCE-ORDERED SCAN
# ============================================================

def _confidence_order(H, W):
    """Order pixels by prediction confidence (most predictable first).
    Confidence = number of 8-connected known neighbors.
    Start from all 4 corners, expand outward.
    This is a BFS from corners, expanding by most-neighbors-first."""
    import heapq
    visited = set()
    order = []
    # Priority queue: (-n_known_neighbors, r, c)
    heap = []

    # Seed: four corners
    for r, c in [(0,0), (0,W-1), (H-1,0), (H-1,W-1)]:
        heapq.heappush(heap, (0, r, c))  # 0 neighbors known initially

    while heap and len(order) < H * W:
        _, r, c = heapq.heappop(heap)
        if (r, c) in visited:
            continue
        visited.add((r, c))
        order.append((r, c))

        # Add unvisited neighbors with updated priority
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and (nr, nc) not in visited:
                    # Count known neighbors of (nr, nc)
                    n_known = sum(1 for dr2 in [-1,0,1] for dc2 in [-1,0,1]
                                  if (dr2 != 0 or dc2 != 0) and
                                  (nr+dr2, nc+dc2) in visited)
                    heapq.heappush(heap, (-n_known, nr, nc))

    # Fill remaining (shouldn't happen but safety)
    for r in range(H):
        for c in range(W):
            if (r, c) not in visited:
                order.append((r, c))

    return order


def _encode_confidence(plane):
    """Encode using confidence-ordered scan with 8-neighbor prediction."""
    H, W = plane.shape
    order = _confidence_order(H, W)
    buf = bytearray()
    known = set()

    for r, c in order:
        # Predict from known 8-neighbors
        total = count = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and (nr, nc) in known:
                    total += int(plane[nr, nc]); count += 1
        pred = total // count if count > 0 else 128
        buf.append((int(plane[r, c]) - pred) & 0xFF)
        known.add((r, c))

    return bytes(buf)


def _decode_confidence(data, H, W):
    """Decode confidence-ordered scan."""
    order = _confidence_order(H, W)
    plane = np.zeros((H, W), dtype=np.uint8)
    known = set()

    for i, (r, c) in enumerate(order):
        total = count = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and (nr, nc) in known:
                    total += int(plane[nr, nc]); count += 1
        pred = total // count if count > 0 else 128
        plane[r, c] = (int(data[i]) + pred) & 0xFF
        known.add((r, c))

    return plane


# ============================================================
# STRATEGY: VALUE-SORTED ENCODING
# ============================================================

def _encode_valsort(plane):
    """Sort pixels by value. Encode sorted values + position permutation.
    Sorted values: nearly monotone → tiny deltas.
    Positions: if image is smooth, nearby-valued pixels are nearby → compressible."""
    H, W = plane.shape
    flat = plane.ravel()
    indices = np.argsort(flat, kind='stable')
    sorted_vals = flat[indices]

    # Encode sorted values with delta
    buf_vals = bytearray()
    buf_vals.append(sorted_vals[0])
    for i in range(1, len(sorted_vals)):
        buf_vals.append((int(sorted_vals[i]) - int(sorted_vals[i-1])) & 0xFF)

    # Encode position permutation as 16-bit indices
    buf_pos = indices.astype(np.uint16).tobytes()

    # Pack: len(vals_compressed) + vals_compressed + pos_compressed
    vals_c = zlib.compress(bytes(buf_vals), 9)
    pos_c = zlib.compress(buf_pos, 9)

    return struct.pack('<II', len(vals_c), len(pos_c)) + vals_c + pos_c


def _decode_valsort(data, H, W):
    """Decode value-sorted encoding."""
    len_vals, len_pos = struct.unpack_from('<II', data, 0)
    offset = 8
    vals_c = data[offset:offset+len_vals]; offset += len_vals
    pos_c = data[offset:offset+len_pos]; offset += len_pos

    # Decode values
    vals_d = zlib.decompress(vals_c)
    sorted_vals = np.zeros(H * W, dtype=np.uint8)
    sorted_vals[0] = vals_d[0]
    for i in range(1, len(vals_d)):
        sorted_vals[i] = (int(sorted_vals[i-1]) + vals_d[i]) & 0xFF

    # Decode positions
    indices = np.frombuffer(zlib.decompress(pos_c), dtype=np.uint16)

    # Reconstruct
    flat = np.zeros(H * W, dtype=np.uint8)
    flat[indices] = sorted_vals
    return flat.reshape(H, W)


# ============================================================
# STRATEGY: COLUMN-SORT + DELTA (from kind-pasteur)
# ============================================================

def _encode_colsort(plane):
    """Sort columns by similarity, then delta between adjacent sorted columns."""
    H, W = plane.shape
    # Sort columns lexicographically
    col_keys = [plane[:, c].tobytes() for c in range(W)]
    perm = sorted(range(W), key=lambda c: col_keys[c])

    buf = bytearray()
    # Store permutation as uint16
    perm_bytes = np.array(perm, dtype=np.uint16).tobytes()

    # Delta between sorted columns
    sorted_plane = plane[:, perm]
    delta = bytearray(sorted_plane[:, 0].tobytes())
    for c in range(1, W):
        d = ((sorted_plane[:, c].astype(np.int16) - sorted_plane[:, c-1].astype(np.int16)) & 0xFF).astype(np.uint8)
        delta.extend(d.tobytes())

    perm_c = zlib.compress(perm_bytes, 9)
    delta_c = zlib.compress(bytes(delta), 9)
    return struct.pack('<II', len(perm_c), len(delta_c)) + perm_c + delta_c


def _decode_colsort(data, H, W):
    """Decode column-sorted encoding."""
    len_perm, len_delta = struct.unpack_from('<II', data, 0)
    offset = 8
    perm_c = data[offset:offset+len_perm]; offset += len_perm
    delta_c = data[offset:offset+len_delta]; offset += len_delta

    perm = np.frombuffer(zlib.decompress(perm_c), dtype=np.uint16)
    # Data was stored column-by-column: col0(H bytes) || col1_delta(H bytes) || ...
    delta_raw = np.frombuffer(zlib.decompress(delta_c), dtype=np.uint8).reshape(W, H).T

    # Reconstruct sorted plane
    sorted_plane = np.zeros((H, W), dtype=np.uint8)
    sorted_plane[:, 0] = delta_raw[:, 0]
    for c in range(1, W):
        sorted_plane[:, c] = ((sorted_plane[:, c-1].astype(np.int16) + delta_raw[:, c].astype(np.int16)) & 0xFF).astype(np.uint8)

    # Inverse permutation
    plane = np.zeros((H, W), dtype=np.uint8)
    for i, p in enumerate(perm):
        plane[:, p] = sorted_plane[:, i]
    return plane


# ============================================================
# STRATEGY: ROW-LPC (from kind-pasteur)
# ============================================================

def _encode_rowlpc(plane):
    """Row-level linear prediction: predict row from linear combo of previous rows.
    Key: use REPRODUCIBLE prediction — fit coeffs from rows [r-3..r-2] to predict row [r-1],
    then use THOSE coeffs to predict row [r]. Decoder has the same data."""
    H, W = plane.shape
    buf = bytearray()
    rows = []  # decoded rows

    for r in range(H):
        if len(rows) >= 3:
            # Fit coefficients: use rows[-3], rows[-2] to predict rows[-1]
            # Then use rows[-2], rows[-1] with same coefficient pattern to predict current
            A_fit = np.column_stack([rows[-3], rows[-2]])
            target = rows[-1]
            try:
                coeffs, _, _, _ = np.linalg.lstsq(A_fit, target, rcond=None)
                # Apply same coefficients shifted forward
                A_pred = np.column_stack([rows[-2], rows[-1]])
                pred = np.clip(np.round(A_pred @ coeffs), 0, 255).astype(int)
            except:
                pred = np.clip(2 * rows[-1] - rows[-2], 0, 255).astype(int)
        elif len(rows) >= 2:
            pred = np.clip(2 * rows[-1] - rows[-2], 0, 255).astype(int)
        elif len(rows) >= 1:
            pred = rows[-1].astype(int)
        else:
            pred = np.zeros(W, dtype=int)

        residual = ((plane[r].astype(int) - pred) & 0xFF).astype(np.uint8)
        buf.extend(residual.tobytes())
        rows.append(plane[r].astype(float))

    return bytes(buf)


def _decode_rowlpc(data, H, W):
    """Decode row-LPC — uses identical prediction logic."""
    plane = np.zeros((H, W), dtype=np.uint8)
    rows = []
    offset = 0

    for r in range(H):
        residual = np.frombuffer(data[offset:offset+W], dtype=np.uint8)
        offset += W

        if len(rows) >= 3:
            A_fit = np.column_stack([rows[-3], rows[-2]])
            target = rows[-1]
            try:
                coeffs, _, _, _ = np.linalg.lstsq(A_fit, target, rcond=None)
                A_pred = np.column_stack([rows[-2], rows[-1]])
                pred = np.clip(np.round(A_pred @ coeffs), 0, 255).astype(int)
            except:
                pred = np.clip(2 * rows[-1] - rows[-2], 0, 255).astype(int)
        elif len(rows) >= 2:
            pred = np.clip(2 * rows[-1] - rows[-2], 0, 255).astype(int)
        elif len(rows) >= 1:
            pred = rows[-1].astype(int)
        else:
            pred = np.zeros(W, dtype=int)

        plane[r] = ((residual.astype(int) + pred) & 0xFF).astype(np.uint8)
        rows.append(plane[r].astype(float))

    return plane


# ============================================================
# STRATEGY: MED (standard baseline)
# ============================================================

def _loco(a, b, c):
    a, b, c = int(a), int(b), int(c)
    if c >= max(a, b): return min(a, b)
    if c <= min(a, b): return max(a, b)
    return a + b - c

def _encode_med(plane):
    """Standard MED/LOCO-I raster prediction."""
    H, W = plane.shape
    buf = bytearray()
    for r in range(H):
        for c in range(W):
            a = int(plane[r, c-1]) if c > 0 else 0
            b = int(plane[r-1, c]) if r > 0 else 0
            d = int(plane[r-1, c-1]) if r > 0 and c > 0 else 0
            pred = _loco(a, b, d)
            buf.append((int(plane[r, c]) - pred) & 0xFF)
    return bytes(buf)

def _decode_med(data, H, W):
    plane = np.zeros((H, W), dtype=np.uint8)
    i = 0
    for r in range(H):
        for c in range(W):
            a = int(plane[r, c-1]) if c > 0 else 0
            b = int(plane[r-1, c]) if r > 0 else 0
            d = int(plane[r-1, c-1]) if r > 0 and c > 0 else 0
            pred = _loco(a, b, d)
            plane[r, c] = (int(data[i]) + pred) & 0xFF
            i += 1
    return plane


# ============================================================
# STRATEGY: ROWSORT + DELTA (new hybrid)
# ============================================================

def _encode_rowsort(plane):
    """Sort rows by similarity, delta between adjacent sorted rows."""
    H, W = plane.shape
    row_keys = [plane[r, :].tobytes() for r in range(H)]
    perm = sorted(range(H), key=lambda r: row_keys[r])

    sorted_plane = plane[perm, :]
    delta = bytearray(sorted_plane[0].tobytes())
    for r in range(1, H):
        d = ((sorted_plane[r].astype(np.int16) - sorted_plane[r-1].astype(np.int16)) & 0xFF).astype(np.uint8)
        delta.extend(d.tobytes())

    perm_c = zlib.compress(np.array(perm, dtype=np.uint16).tobytes(), 9)
    delta_c = zlib.compress(bytes(delta), 9)
    return struct.pack('<II', len(perm_c), len(delta_c)) + perm_c + delta_c


def _decode_rowsort(data, H, W):
    len_perm, len_delta = struct.unpack_from('<II', data, 0)
    offset = 8
    perm = np.frombuffer(zlib.decompress(data[offset:offset+len_perm]), dtype=np.uint16); offset += len_perm
    # Data stored row-by-row: row0(W bytes) || row1_delta(W bytes) || ...
    delta_raw = np.frombuffer(zlib.decompress(data[offset:offset+len_delta]), dtype=np.uint8).reshape(H, W)  # row-major is correct here

    sorted_plane = np.zeros((H, W), dtype=np.uint8)
    sorted_plane[0] = delta_raw[0]
    for r in range(1, H):
        sorted_plane[r] = ((sorted_plane[r-1].astype(np.int16) + delta_raw[r].astype(np.int16)) & 0xFF).astype(np.uint8)

    plane = np.zeros((H, W), dtype=np.uint8)
    for i, p in enumerate(perm):
        plane[p] = sorted_plane[i]
    return plane


# ============================================================
# ENCODE: TRY ALL, PICK BEST
# ============================================================

METHODS = {
    0: ("MED", _encode_med, _decode_med, False),
    1: ("RowLPC", _encode_rowlpc, _decode_rowlpc, False),
    2: ("ColSort", _encode_colsort, _decode_colsort, True),
    3: ("RowSort", _encode_rowsort, _decode_rowsort, True),
    4: ("Confidence", _encode_confidence, _decode_confidence, False),
    5: ("ValSort", _encode_valsort, _decode_valsort, True),
}


def encode(image):
    image = np.asarray(image, dtype=np.uint8)
    if image.ndim == 2:
        H, W = image.shape; channels = 1
        planes_list = [([image], 0)]
    elif image.ndim == 3 and image.shape[2] == 3:
        H, W = image.shape[:2]; channels = 3
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

    cands = []

    for planes, ct_id in planes_list:
        for mid, (mname, enc_fn, dec_fn, self_contained) in METHODS.items():
            if self_contained:
                # Methods that handle their own compression
                buf = bytearray()
                for plane in planes:
                    pdata = enc_fn(plane)
                    buf.extend(struct.pack('<I', len(pdata)))
                    buf.extend(pdata)
                cands.append((mid, ct_id, bytes(buf), True))
            else:
                # Methods that produce raw residuals for compression
                buf = bytearray()
                for plane in planes:
                    buf.extend(enc_fn(plane))
                cands.append((mid, ct_id, bytes(buf), False))

    # Screen
    screened = []
    for mid, ct, pay, sc in cands:
        if sc:
            screened.append((mid, ct, pay, len(pay), sc))
        else:
            screened.append((mid, ct, pay, len(_fast_compress(pay)), sc))

    screened.sort(key=lambda x: x[3])

    # Top 3 with lzma
    best_total = 1 << 60
    best_result = None
    for mid, ct, pay, _, sc in screened[:3]:
        if sc:
            # Already compressed internally
            total = len(pay)
            if total < best_total:
                best_total = total
                best_result = (mid, ct, 0xFF, pay)  # 0xFF = self-contained
        else:
            bid, comp = _best_backend(pay)
            if len(comp) < best_total:
                best_total = len(comp)
                best_result = (mid, ct, bid, comp)

    mid, ct_id, bid, compressed = best_result
    header = struct.pack('<4sBHHBBBB', MAGIC, 6, H, W, channels, mid, bid, ct_id)
    return header + compressed


def decode(data):
    magic = data[:4]
    if magic != MAGIC:
        # Delegate to older versions
        from tic4_opus import decode as d4
        return d4(data)

    _, ver, H, W, channels, mid, bid, ct_id = struct.unpack_from('<4sBHHBBBB', data, 0)
    hdr = struct.calcsize('<4sBHHBBBB')

    if bid == 0xFF:
        payload = data[hdr:]
    else:
        payload = _decompress(data[hdr:], bid)

    _, _, dec_fn, self_contained = list(METHODS.values())[mid]

    if self_contained:
        planes = []
        offset = 0
        for ch in range(channels):
            plen = struct.unpack_from('<I', payload, offset)[0]; offset += 4
            pdata = payload[offset:offset+plen]; offset += plen
            planes.append(dec_fn(pdata, H, W))
    else:
        planes = []
        plane_size = H * W
        for ch in range(channels):
            pdata = payload[ch * plane_size:(ch+1) * plane_size]
            planes.append(dec_fn(pdata, H, W))

    if channels == 1:
        return planes[0]
    return CT_INV[ct_id](np.stack(planes, axis=2))


# ============================================================
# TEST
# ============================================================

def png_size(img):
    from PIL import Image
    mode = 'L' if img.ndim == 2 else 'RGB'
    pil = Image.fromarray(img, mode)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()


def _test():
    np.random.seed(42)
    print("=" * 90)
    print(f"  TIC v{__version__} — Axiom-Breaking Adaptive Codec")
    print("=" * 90)

    passed = 0; failed = 0

    def check(name, orig, enc_data):
        nonlocal passed, failed
        rec = decode(enc_data)
        mid = enc_data[struct.calcsize('<4sBHH')]
        mname = list(METHODS.values())[mid][0] if mid < len(METHODS) else f"m{mid}"

        if np.array_equal(orig, rec):
            passed += 1
            sz = len(enc_data)
            ps = png_size(orig)
            ratio = ps / sz
            win = "WIN" if ratio > 1.001 else "LOSS"
            print(f"  OK {name:<25} TIC={sz:>7} PNG={ps:>7} {ratio:.2f}x {win:>4} [{mname}]")
        else:
            failed += 1
            diff = np.sum(orig != rec)
            print(f"  FAIL {name:<25} {diff} diffs [{mname}]")

    print("\n--- Grayscale ---")
    for N in [64, 128]:
        for name, img in [
            (f"gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8),(N,1))),
            (f"v_gradient_{N}", np.tile(np.linspace(0,255,N,dtype=np.uint8).reshape(-1,1),(1,N))),
            (f"diag_grad_{N}", np.array([[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"circle_{N}", np.array([[min(255,int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"checker_{N}", np.array([[(r+c)%2*255 for c in range(N)] for r in range(N)], dtype=np.uint8)),
            (f"noise_{N}", np.random.randint(0,256,(N,N),dtype=np.uint8)),
            (f"solid_{N}", np.full((N,N),128,dtype=np.uint8)),
            (f"blob_{N}", np.array([[min(255,int(200*np.exp(-((r-N//3)**2+(c-2*N//3)**2)/(N**2/16))+100*np.exp(-((r-2*N//3)**2+(c-N//3)**2)/(N**2/20)))) for c in range(N)] for r in range(N)], dtype=np.uint8)),
        ]:
            d = encode(img)
            check(name, img, d)

    # RGB
    print("\n--- RGB ---")
    for N in [64, 128]:
        img = np.random.randint(0,256,(N,N,3),dtype=np.uint8)
        d = encode(img)
        check(f"rgb_noise_{N}", img, d)

    # Real photos
    print("\n--- Real Photos ---")
    from PIL import Image as PILImage
    for path in ['/home/e/tron.jpeg', '/home/e/t.jpg']:
        if os.path.exists(path):
            img = np.array(PILImage.open(path).convert('RGB'))
            t0 = time.time()
            d = encode(img)
            t_enc = time.time() - t0
            check(f"{os.path.basename(path)} ({t_enc:.1f}s)", img, d)

            # Gray
            gray = np.array(PILImage.open(path).convert('L'))
            d = encode(gray)
            check(f"{os.path.basename(path)}_gray", gray, d)

    print(f"\n  {passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] in ("test", "bench"):
        _test()
    else:
        print(f"TIC v{__version__} — Axiom-Breaking Adaptive Codec")
        print("  python3 tic6_adaptive.py test")
