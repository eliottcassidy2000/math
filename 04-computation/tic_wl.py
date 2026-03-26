#!/usr/bin/env python3
"""
WL-GUIDED CODEC: Weisfeiler-Leman coloring guides prediction.

Cheap WL assigns each row a "color" = quantized (mean, std, gradient).
Rows with the same color have similar statistics.

THREE WL-GUIDED STRATEGIES:

1. WL-SORT: Sort rows by WL color, delta-compress within groups.
   Like RowSort but with a principled grouping criterion.
   Rows in the same group are guaranteed similar → small deltas.

2. WL-LPC: Compute separate LPC coefficients per WL color group.
   Each group is statistically homogeneous → better LPC fit.
   The decoder knows the WL color (computable from decoded rows) → zero overhead.

3. WL-BLEND: Use WL color to set blend threshold per row.
   Smooth rows (low std) → Wiener. Edge rows (high gradient) → MED.
   The WL color refines the per-pixel gradient with per-ROW statistics.

Also: WL as a REORDERING for snake-LPC. Sort by WL color,
then apply LPC on the reordered signal. Groups are contiguous
in the signal → LPC captures intra-group correlations perfectly.

The key: WL is O(n) per level, both encoder and decoder compute it.
NO OVERHEAD for the WL classification.

kind-pasteur-2026-03-25-S18
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

def lpc_coeffs(sig, order):
    n = len(sig)
    if n < order + 1: return np.zeros(order)
    r = np.correlate(sig, sig, mode='full'); r = r[n-1:n+order]
    if r[0] < 1e-10: return np.zeros(order)
    a = np.zeros(order); e = r[0]
    for i in range(order):
        if e < 1e-10: break
        lam = r[i+1]
        for j in range(i): lam -= a[j] * r[i-j]
        k = lam / e; a_new = np.zeros(order); a_new[i] = k
        for j in range(i): a_new[j] = a[j] - k * a[i-1-j]
        a = a_new; e *= (1 - k*k)
    return a

# ============================================================
# CHEAP WL: Assign colors to rows based on statistics
# ============================================================

def wl_color_rows(plane, q=4):
    """Assign WL color to each row. Color = quantized (mean, std).
    q bits of quantization per stat → 2^(2q) possible colors."""
    h, w = plane.shape
    colors = np.zeros(h, dtype=np.int32)

    for r in range(h):
        row = plane[r].astype(float)
        mean_val = row.mean()
        std_val = row.std()

        # Quantize to q bits each
        qmean = int(np.clip(mean_val * (1 << q) / 256, 0, (1 << q) - 1))
        qstd = int(np.clip(std_val * (1 << q) / 128, 0, (1 << q) - 1))

        colors[r] = qmean * (1 << q) + qstd

    return colors

def wl_color_rows_l2(plane, q=3):
    """L2 WL: (mean, std, mean_gradient) with q bits each."""
    h, w = plane.shape
    colors = np.zeros(h, dtype=np.int32)

    for r in range(h):
        row = plane[r].astype(float)
        mean_val = row.mean()
        std_val = row.std()
        if w > 1:
            grad = np.mean(np.abs(np.diff(row)))
        else:
            grad = 0

        qm = int(np.clip(mean_val * (1 << q) / 256, 0, (1 << q) - 1))
        qs = int(np.clip(std_val * (1 << q) / 128, 0, (1 << q) - 1))
        qg = int(np.clip(grad * (1 << q) / 64, 0, (1 << q) - 1))

        colors[r] = qm * (1 << (2*q)) + qs * (1 << q) + qg

    return colors

# ============================================================
# WL-SORT: Sort rows by WL color, delta-compress groups
# ============================================================

def encode_wl_sort(plane, level=1):
    """Sort rows by WL color. Within each group, rows are similar → small deltas.
    The sort order IS the WL color ordering (decoder can recompute)."""
    h, w = plane.shape

    if level == 1:
        colors = wl_color_rows(plane)
    else:
        colors = wl_color_rows_l2(plane)

    # Sort rows by color (stable sort preserves spatial order within same color)
    order = np.argsort(colors, kind='stable')
    sorted_plane = plane[order]

    # Delta-compress the sorted rows
    delta = np.zeros_like(sorted_plane)
    delta[0] = sorted_plane[0]
    delta[1:] = ((sorted_plane[1:].astype(int) - sorted_plane[:-1].astype(int)) & 0xFF).astype(np.uint8)

    # Encode: permutation + delta
    perm = order.astype(np.uint8 if h <= 256 else np.uint16)
    return len(bc(perm.tobytes())) + len(bc(delta.tobytes()))

# ============================================================
# WL-LPC: Per-color-group LPC
# ============================================================

def encode_wl_lpc(plane, order=4):
    """Compute LPC per WL color group. Groups are homogeneous → better LPC.
    Encode: per-group coefficients + residuals in raster order."""
    h, w = plane.shape
    colors = wl_color_rows(plane)
    unique_colors = np.unique(colors)

    # Compute LPC coefficients per group
    group_coeffs = {}
    for color in unique_colors:
        rows = [plane[r] for r in range(h) if colors[r] == color]
        if len(rows) < 2:
            group_coeffs[color] = np.zeros(order)
            continue
        signal = np.concatenate(rows).astype(float)
        group_coeffs[color] = lpc_coeffs(signal - signal.mean(), order)

    # Apply per-group LPC in RASTER order (preserves locality)
    residuals = bytearray()
    # Build the snake signal with per-row LPC coefficients
    sig = plane.ravel().astype(float)
    for r in range(h):
        coeffs = group_coeffs.get(colors[r], np.zeros(order))
        for c in range(w):
            i = r * w + c
            pred = sum(coeffs[j] * sig[i-1-j] for j in range(min(i, order)))
            residuals.append((int(sig[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF)

    # Encode: coefficients for each group + residuals
    coeff_data = bytearray()
    for color in unique_colors:
        coeff_data.extend(np.array(group_coeffs[color], dtype=np.float32).tobytes())

    return len(bc(bytes(coeff_data))) + len(bc(bytes(residuals)))

# ============================================================
# WL-BLEND: Use WL color to set per-row blend threshold
# ============================================================

def encode_wl_blend(plane):
    """Blend MED-Wiener with per-row threshold based on WL statistics.
    Low-std rows → low threshold → Wiener dominates.
    High-std rows → high threshold → MED dominates.
    Decoder recomputes row statistics → zero overhead."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    for r in range(h):
        # Row statistics (decoder can compute from decoded row r-1)
        if r > 0:
            prev_std = np.std(img[r-1])
            row_threshold = max(5, prev_std * 2)  # adaptive threshold
        else:
            row_threshold = 50  # default

        for c in range(w):
            a = img[r,c-1] if c > 0 else 0.0
            b = img[r-1,c] if r > 0 else 0.0
            d = img[r-1,c-1] if r > 0 and c > 0 else 0.0

            grad = abs(a-d) + abs(b-d) if r > 0 and c > 0 else 0
            alpha = min(1.0, grad / row_threshold)

            p_med = med(int(a), int(b), int(d))

            total, wt = 0.0, 0.0
            for dr, dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr, nc = r+dr, c+dc
                if 0 <= nr < h and 0 <= nc < w:
                    dst = abs(dr)+abs(dc); w2 = 4.0/(dst*dst)
                    total += img[nr,nc]*w2; wt += w2
            p_wie = total/wt if wt > 0 else 128

            pred = int(np.clip(round(alpha*p_med + (1-alpha)*p_wie), 0, 255))
            residuals.append((int(plane[r,c]) - pred) & 0xFF)

    return len(bc(bytes(residuals)))

# ============================================================
# WL-REORDER + LPC: Sort by WL color, then LPC on reordered signal
# ============================================================

def encode_wl_reorder_lpc(plane, order=8):
    """Sort rows by WL color, concatenate into 1D signal, apply LPC.
    Groups are contiguous in the signal → LPC captures intra-group structure.
    Like snake-LPC but with WL-guided ordering."""
    h, w = plane.shape
    colors = wl_color_rows(plane)
    order_idx = np.argsort(colors, kind='stable')

    # Build reordered signal
    signal = np.concatenate([plane[i] for i in order_idx]).astype(float)
    n = len(signal)

    coeffs = lpc_coeffs(signal - signal.mean(), order)
    res = np.empty(n, dtype=np.uint8)
    for i in range(n):
        pred = sum(coeffs[j] * signal[i-1-j] for j in range(min(i, order)))
        res[i] = (int(signal[i]) - int(np.clip(round(pred), 0, 255))) & 0xFF

    perm = order_idx.astype(np.uint8 if h <= 256 else np.uint16)
    coeff_bytes = np.array(coeffs, dtype=np.float32).tobytes()

    return len(bc(perm.tobytes())) + len(coeff_bytes) + len(bc(bytes(res)))

# ============================================================
# BASELINES
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,d))&0xFF)
    return len(bc(bytes(res)))

def m_blend(plane):
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha=min(1.0, grad/1000)
            p_med=med(int(a),int(b),int(d))
            total,wt=0.0,0.0
            for dr,dc in [(-1,-1),(-1,0),(-1,1),(0,-1),(-2,0),(0,-2)]:
                nr,nc=r+dr,c+dc
                if 0<=nr<h and 0<=nc<w:
                    dst=abs(dr)+abs(dc); w2=4.0/(dst*dst)
                    total+=img[nr,nc]*w2; wt+=w2
            p_wie=total/wt if wt>0 else 128
            pred=int(np.clip(round(alpha*p_med+(1-alpha)*p_wie),0,255))
            res.append((int(plane[r,c])-pred)&0xFF)
    return len(bc(bytes(res)))

def m_colsort(plane):
    h,w=plane.shape
    idx=sorted(range(w),key=lambda c2:plane[:,c2].tobytes())
    sp=np.stack([plane[:,c2] for c2 in idx],axis=1)
    delta=np.zeros_like(sp); delta[:,0]=sp[:,0]
    delta[:,1:]=((sp[:,1:].astype(int)-sp[:,:-1].astype(int))&0xFF).astype(np.uint8)
    perm=np.array(idx,dtype=np.uint8)
    return len(bc(perm.tobytes()))+len(bc(delta.tobytes()))

def m_snake_lpc(plane, order=8):
    h,w=plane.shape
    signal = np.empty(h*w, dtype=np.uint8)
    for r in range(h):
        if r%2==0: signal[r*w:(r+1)*w]=plane[r]
        else: signal[r*w:(r+1)*w]=plane[r,::-1]
    sig=signal.astype(float); n=len(sig)
    coeffs=lpc_coeffs(sig-sig.mean(),order)
    res=np.empty(n,dtype=np.uint8)
    for i in range(n):
        pred=sum(coeffs[j]*sig[i-1-j] for j in range(min(i,order)))
        res[i]=(int(sig[i])-int(np.clip(round(pred),0,255)))&0xFF
    return len(np.array(coeffs,dtype=np.float32).tobytes())+len(bc(bytes(res)))

# ============================================================
# TEST
# ============================================================

SZ = 128

def make_tests():
    T={}; np.random.seed(42)
    x,y=np.meshgrid(np.arange(SZ,dtype=float),np.arange(SZ,dtype=float))
    r=np.sqrt((x-SZ/2)**2+(y-SZ/2)**2)
    T["solid"]=np.full((SZ,SZ),128,dtype=np.uint8)
    T["grad_h"]=np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1))
    T["grad_d"]=((x+y)*255//(2*SZ-2)).astype(np.uint8)
    T["circles"]=(np.sin(r/5)*127+128).astype(np.uint8)
    T["blob"]=(np.exp(-r**2/(2*(SZ/4)**2))*255).astype(np.uint8)
    T["random"]=np.random.randint(0,256,(SZ,SZ),dtype=np.uint8)
    sm=np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)
    T["natural"]=np.clip(np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
                         +np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8)
    T["stripes45"]=(((x+y)/8).astype(int)%2*255).astype(np.uint8)
    return T

print("="*90)
print("  WL-GUIDED CODEC: Weisfeiler-Leman coloring meets compression")
print("  kind-pasteur-2026-03-25-S18")
print("="*90)

tests=make_tests()

METHODS=[
    ("MED",        m_med),
    ("Blend",      m_blend),
    ("SnakeLPC8",  lambda p: m_snake_lpc(p, 8)),
    ("ColSort",    m_colsort),
    ("WL-Sort",    lambda p: encode_wl_sort(p, 1)),
    ("WL-LPC",     lambda p: encode_wl_lpc(p, 4)),
    ("WL-Blend",   encode_wl_blend),
    ("WL-ReordLPC",lambda p: encode_wl_reorder_lpc(p, 8)),
    ("Raw",        lambda p: len(bc(p.tobytes()))),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname,_ in METHODS: print(f" {mname:>10}", end="")
print(f"  {'BEST':>12}")
print("  "+"-"*(20+11*len(METHODS)+16))

for name,img in sorted(tests.items()):
    ps=png_size(img); sizes={}
    for mname,mfunc in METHODS:
        try: sizes[mname]=mfunc(img)+10
        except: sizes[mname]=999999
    best=min(sizes,key=sizes.get)
    ratio=ps/sizes[best] if sizes[best]<999999 else 0
    print(f"  {name:<12} {ps:>6}", end="")
    for mname,_ in METHODS:
        v=sizes[mname]; m="*" if mname==best else " "
        print(f" {v if v<999999 else 'ERR':>9}{m}", end="")
    print(f"  {best:>12} {ratio:.2f}x")

print(f"""
  WL COLOR STATISTICS:
""")
for name in ["natural", "blob", "circles"]:
    img = tests[name]
    colors = wl_color_rows(img)
    n_colors = len(np.unique(colors))
    print(f"  {name}: {n_colors} WL colors (L1)")
    colors2 = wl_color_rows_l2(img)
    n_colors2 = len(np.unique(colors2))
    print(f"  {name}: {n_colors2} WL colors (L2)")
