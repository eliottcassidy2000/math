#!/usr/bin/env python3
"""
RIGOROUS OVERHEAD ANALYSIS: Is block classification fundamentally limited?

The probe codec failed because "classification overhead" ate the gains.
But is this FUNDAMENTAL or just a bad implementation?

To prove it's fundamental, I need to show:
  For ANY block classification with k classes on B blocks,
  the classification map costs at least B * log2(k) / 8 bytes,
  and the maximum prediction gain from knowing the class is bounded.

To prove it's FIXABLE, I need to show:
  There EXISTS a classification scheme where:
  classification_cost < prediction_gain_from_classification

Let me measure EXACTLY:
1. How many bits does the classification map cost? (after compression)
2. How many bits does knowing the classification SAVE in prediction?
3. Is the net (savings - cost) positive?

If positive: the probe codec CAN work with better implementation.
If negative: it's fundamental and we should abandon it.

THE KEY INSIGHT: The classification map is ITSELF an image (of block types).
It has spatial structure (adjacent blocks are usually the same type).
If it compresses well, the overhead shrinks.

kind-pasteur-2026-03-25-S14
"""
import sys, io, zlib, math, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

try: import brotli; HAS_BROTLI = True
except: HAS_BROTLI = False

def bc(data):
    r = [zlib.compress(data, 9)]
    if HAS_BROTLI:
        try: r.append(brotli.compress(data, quality=11))
        except: pass
    return min(r, key=len)

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
# STEP 1: Measure the EXACT overhead of block classification
# ============================================================

def measure_classification_overhead(plane, bs=8, n_classes=4):
    """Compute the compressed size of the classification map.
    This is the MINIMUM cost of block-level routing."""
    h, w = plane.shape
    bh, bw = (h + bs - 1) // bs, (w + bs - 1) // bs

    # Classify blocks
    classes = np.zeros((bh, bw), dtype=np.uint8)
    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by * bs, bx * bs
            r1, c1 = min(r0 + bs, h), min(c0 + bs, w)
            block = plane[r0:r1, c0:c1].astype(float)
            var_val = block.var()
            if var_val < 2: classes[by, bx] = 0      # flat
            elif var_val < 25: classes[by, bx] = 1    # smooth
            elif var_val < 100: classes[by, bx] = 2   # edge
            else: classes[by, bx] = 3                  # texture

    # Compressed size of classification map
    raw_map_bytes = bh * bw  # 1 byte per block
    compressed_map = bc(classes.tobytes())

    # Entropy of classification map
    unique, counts = np.unique(classes, return_counts=True)
    total = bh * bw
    entropy_bits = -sum((c/total) * math.log2(c/total) for c in counts if c > 0) * total

    return len(compressed_map), entropy_bits / 8, raw_map_bytes, classes

# ============================================================
# STEP 2: Measure EXACT prediction gain from knowing the class
# ============================================================

def measure_prediction_gain(plane, classes, bs=8):
    """For each block, compute:
    - Size with MED (one method for all)
    - Size with OPTIMAL method per class
    The difference is the GAIN from classification."""
    h, w = plane.shape
    bh, bw = classes.shape

    # Method 0: MED on everything (baseline)
    med_residuals = bytearray()
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            med_residuals.append((int(plane[r,c]) - med(a,b,d)) & 0xFF)
    med_size = len(bc(bytes(med_residuals)))

    # Method 1: Per-class encoding
    # Collect residuals by class, compress each class separately
    class_residuals = {0: bytearray(), 1: bytearray(), 2: bytearray(), 3: bytearray()}

    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by * bs, bx * bs
            r1, c1 = min(r0 + bs, h), min(c0 + bs, w)
            block = plane[r0:r1, c0:c1]
            cls = classes[by, bx]

            if cls == 0:  # flat: store delta from mean
                mean_val = int(round(block.mean()))
                for r in range(r0, r1):
                    for c in range(c0, c1):
                        class_residuals[0].append((int(plane[r,c]) - mean_val) & 0xFF)
            elif cls == 1:  # smooth: delta from left
                for r in range(r0, r1):
                    for c in range(c0, c1):
                        pred = int(plane[r, c-1]) if c > c0 else int(plane[r-1, c]) if r > r0 else 0
                        class_residuals[1].append((int(plane[r,c]) - pred) & 0xFF)
            elif cls == 2:  # edge: MED
                for r in range(r0, r1):
                    for c in range(c0, c1):
                        a = int(plane[r,c-1]) if c > 0 else 0
                        b = int(plane[r-1,c]) if r > 0 else 0
                        d2 = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
                        class_residuals[2].append((int(plane[r,c]) - med(a,b,d2)) & 0xFF)
            else:  # texture: raw
                for r in range(r0, r1):
                    for c in range(c0, c1):
                        class_residuals[3].append(int(plane[r,c]))

    per_class_size = sum(len(bc(bytes(v))) for v in class_residuals.values() if v)

    return med_size, per_class_size

# ============================================================
# STEP 3: The fundamental question
# ============================================================

SZ = 128

def make_tests():
    T = {}; np.random.seed(42)
    x, y = np.meshgrid(np.arange(SZ, dtype=float), np.arange(SZ, dtype=float))
    r = np.sqrt((x-SZ/2)**2 + (y-SZ/2)**2)
    T["solid"] = np.full((SZ,SZ), 128, dtype=np.uint8)
    T["grad_h"] = np.tile(np.linspace(0, 255, SZ, dtype=np.uint8), (SZ, 1))
    T["circles"] = (np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"] = (np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"] = np.random.randint(0, 256, (SZ,SZ), dtype=np.uint8)
    sm = np.random.randint(0, 256, (SZ//16, SZ//16), dtype=np.uint8)
    T["natural"] = np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ), Image.BILINEAR)).astype(float)
                           + np.random.normal(0, 10, (SZ,SZ)), 0, 255).astype(np.uint8)
    # THE CRITICAL TEST: mixed content (where classification should help)
    mixed = np.zeros((SZ, SZ), dtype=np.uint8)
    mixed[:SZ//2, :SZ//2] = 128  # flat
    mixed[:SZ//2, SZ//2:] = np.tile(np.linspace(0, 255, SZ//2, dtype=np.uint8), (SZ//2, 1))
    mixed[SZ//2:, :SZ//2] = np.random.randint(0, 256, (SZ//2, SZ//2), dtype=np.uint8)
    r2 = np.sqrt((x[SZ//2:, SZ//2:]-SZ*0.75)**2 + (y[SZ//2:, SZ//2:]-SZ*0.75)**2)
    mixed[SZ//2:, SZ//2:] = (np.sin(r2/5)*127+128).astype(np.uint8)
    T["mixed"] = mixed
    return T

print("=" * 90)
print("  RIGOROUS OVERHEAD ANALYSIS: Is block classification fundamentally limited?")
print("  kind-pasteur-2026-03-25-S14")
print("=" * 90)

tests = make_tests()

print(f"\n  STEP 1: Classification map overhead")
print(f"  {'Image':<12} {'raw_map':>8} {'comp_map':>9} {'entropy':>8} {'#blocks':>8} {'class_dist':>30}")
print("  " + "-" * 80)

for name, img in sorted(tests.items()):
    comp_map, entropy_bytes, raw_bytes, classes = measure_classification_overhead(img)
    unique, counts = np.unique(classes, return_counts=True)
    dist = {0:'F',1:'S',2:'E',3:'T'}
    dist_str = ' '.join(f"{dist[u]}:{c}" for u, c in zip(unique, counts))
    print(f"  {name:<12} {raw_bytes:>8} {comp_map:>9} {entropy_bytes:>8.1f} {classes.size:>8} {dist_str:>30}")

print(f"\n  STEP 2: Prediction gain from classification")
print(f"  {'Image':<12} {'MED_all':>8} {'per_class':>10} {'gain':>6} {'map_cost':>9} {'NET':>6} {'verdict':>10}")
print("  " + "-" * 70)

for name, img in sorted(tests.items()):
    comp_map, _, _, classes = measure_classification_overhead(img)
    med_sz, class_sz = measure_prediction_gain(img, classes)
    gain = med_sz - class_sz
    net = gain - comp_map
    verdict = "WORTH IT" if net > 0 else "NOT WORTH"
    print(f"  {name:<12} {med_sz:>8} {class_sz:>10} {gain:>+6} {comp_map:>9} {net:>+6} {verdict:>10}")

print(f"""
  STEP 3: THE FUNDAMENTAL QUESTION

  Is block classification fundamentally limited?

  THEOREM ATTEMPT: For an image where all blocks have the same statistics
  (homogeneous image), classification provides ZERO prediction gain.
  But the map still costs > 0 bytes. So net < 0 always.

  COUNTER-THEOREM: For a mixed image (flat + gradient + noise + circles),
  classification provides prediction gain proportional to the entropy
  difference between the best single predictor and the per-class predictors.
  If this exceeds the map cost, net > 0.

  THE ANSWER depends on the image:
  - Homogeneous images: classification is PROVABLY useless
  - Mixed images: classification CAN help IF the map compresses well
    AND the prediction gain exceeds the map cost

  THE CRITICAL RATIO:
    classification_gain / map_cost = prediction_entropy_reduction / map_entropy

  For this to be > 1, the prediction gain per block must exceed
  the cost of identifying that block's class (~2 bits for 4 classes).

  A block of 8×8 = 64 pixels. 2 bits of class info per block =
  2/64 = 0.03 bits per pixel overhead. If knowing the class saves
  more than 0.03 bpp in prediction, it's worth it.
""")

# ============================================================
# STEP 4: What does this teach us about FIXING the probe codec?
# ============================================================

print("  STEP 4: How to FIX the probe codec")
print()

# The key: if per-block LPC is bad because blocks are too small,
# what about per-REGION LPC? Contiguous regions of the same class
# can be encoded as one long signal.

# Test: concatenate all blocks of the same class into one signal,
# apply LPC to the concatenated signal.

for name in ["mixed", "natural", "blob"]:
    img = tests[name]
    h, w = img.shape
    _, _, _, classes = measure_classification_overhead(img, bs=8)
    bh, bw = classes.shape

    print(f"  {name}: region-level LPC")

    # Concatenate blocks by class
    class_signals = {0: [], 1: [], 2: [], 3: []}
    for by in range(bh):
        for bx in range(bw):
            r0, c0 = by * 8, bx * 8
            r1, c1 = min(r0+8, h), min(c0+8, w)
            cls = classes[by, bx]
            class_signals[cls].extend(img[r0:r1, c0:c1].ravel())

    # LPC on each class's concatenated signal
    from scipy.signal import lfilter

    total_class = 0
    for cls in [0, 1, 2, 3]:
        sig = class_signals[cls]
        if not sig:
            continue
        sig = np.array(sig, dtype=float)
        n = len(sig)

        if n < 10:
            total_class += len(bc(np.array(sig, dtype=np.uint8).tobytes()))
            continue

        # LPC order 4
        order = min(4, n-1)
        r2 = np.correlate(sig - sig.mean(), sig - sig.mean(), mode='full')
        r2 = r2[n-1:n+order]
        if r2[0] < 1e-10:
            res = sig.astype(np.uint8)
        else:
            a = np.zeros(order); e = r2[0]
            for i in range(order):
                if e < 1e-10: break
                lam = r2[i+1]
                for j in range(i): lam -= a[j] * r2[i-j]
                k = lam / e
                a_new = np.zeros(order); a_new[i] = k
                for j in range(i): a_new[j] = a[j] - k * a[i-1-j]
                a = a_new; e *= (1 - k*k)

            res = np.empty(n, dtype=np.uint8)
            for i in range(n):
                pred = sum(a[j] * sig[i-1-j] for j in range(min(i, order)))
                res[i] = (int(sig[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF

        csz = len(bc(res.tobytes()))
        cls_names = {0: 'FLAT', 1: 'SMOOTH', 2: 'EDGE', 3: 'TEXTURE'}
        print(f"    {cls_names[cls]:>7}: {n:>6} pixels → {csz:>6} bytes ({csz*8/n:.2f} bpp)")
        total_class += csz

    map_cost, _, _, _ = measure_classification_overhead(img, bs=8)
    total_region = total_class + map_cost

    med_sz = 0
    res2 = bytearray()
    for r in range(h):
        for c in range(w):
            a2 = int(img[r,c-1]) if c>0 else 0; b2 = int(img[r-1,c]) if r>0 else 0
            d2 = int(img[r-1,c-1]) if r>0 and c>0 else 0
            res2.append((int(img[r,c])-med(a2,b2,d2))&0xFF)
    med_sz = len(bc(bytes(res2)))

    print(f"    TOTAL: class_data={total_class} + map={map_cost} = {total_region} vs MED={med_sz}")
    print(f"    {'REGION LPC WINS' if total_region < med_sz else 'MED STILL WINS'} (diff={total_region - med_sz:+d})")
    print()
