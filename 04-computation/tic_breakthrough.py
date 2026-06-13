#!/usr/bin/env python3
"""
BREAKTHROUGH ATTEMPT: Merge tournament theory insights with signal processing.

THREE IDEAS COMBINED:

1. ADAPTIVE WIENER FILTER: Maintain running autocorrelation estimate.
   Update per-pixel. Always use the LOCALLY OPTIMAL linear predictor.
   Combines Wiener's optimality with LPC's adaptivity.

2. SCORE-CONDITIONED RESIDUAL ENCODING: After prediction, residuals
   are nearly Laplacian. The SUM of absolute residuals per block
   constrains the pattern. Encode sum first (cheap), then pattern
   within that sum class (combinadic-like via zlib on pre-sorted data).

3. CONTEXT-ADAPTIVE PREDICTION ORDER: Use gradient magnitude to decide
   HOW MANY neighbors to use. Smooth regions: use many (Wiener with
   large radius = lower variance prediction). Edge regions: use few
   (MED with 3 neighbors = nonlinear edge detection).

The tournament theory connection:
- Score conditioning = Hamming weight prior on residuals (from tournament theory)
- Adaptive prediction = waggly layer concept (different scales for different d)
- Context selection = silent vs expressive mutation (smooth vs edge regions)

kind-pasteur-2026-03-25-S10
"""
import sys, io, struct, zlib, time, math
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

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            c2=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,c2))&0xFF)
    return len(bc(bytes(res)))

# ============================================================
# IDEA 1: ADAPTIVE WIENER (sliding window autocorrelation)
# ============================================================

def encode_adaptive_wiener(plane, win=32):
    """Adaptive Wiener filter with sliding window.
    Maintains running autocorrelation estimate, updates per pixel.
    Uses the CURRENT local statistics to predict each pixel."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    # Causal 2D neighborhood offsets (within scan order)
    offsets = []
    for dr in range(-2, 1):
        for dc in range(-2, 3):
            if dr == 0 and dc >= 0: continue
            if abs(dr) + abs(dc) <= 3:
                offsets.append((dr, dc))
    nf = len(offsets)

    # Running autocorrelation matrix (R) and cross-correlation (p)
    # Updated with exponential decay
    alpha = 0.99  # decay factor (higher = more memory)
    R = np.eye(nf) * 0.01  # small regularization
    p = np.zeros(nf)
    weights = np.zeros(nf)

    for r in range(h):
        for c in range(w):
            # Gather neighbor values
            feats = np.zeros(nf)
            valid = True
            for k, (dr, dc) in enumerate(offsets):
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w:
                    feats[k] = img[nr, nc]
                else:
                    feats[k] = 128  # default

            # Predict using current weights
            pred = np.dot(weights, feats)
            pred = np.clip(round(pred), 0, 255)

            actual = int(plane[r, c])
            residuals.append((actual - int(pred)) & 0xFF)

            # Update running autocorrelation
            R = alpha * R + np.outer(feats, feats)
            p = alpha * p + actual * feats

            # Recompute weights periodically (every w pixels = once per row)
            if c == w - 1 or (c % 32 == 31):
                try:
                    weights = np.linalg.solve(R + 1e-6 * np.eye(nf), p)
                except:
                    pass  # keep previous weights

    return len(bc(bytes(residuals)))

# ============================================================
# IDEA 2: CONTEXT-SWITCHING PREDICTOR (smooth → Wiener, edge → MED)
# ============================================================

def encode_context_switch(plane):
    """Switch between Wiener (smooth) and MED (edge) per pixel.
    Use gradient magnitude to classify context."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    # Precompute gradient magnitude
    grad = np.zeros((h, w))
    for r in range(1, h):
        for c in range(1, w):
            gx = abs(img[r, c-1] - img[r-1, c-1])
            gy = abs(img[r-1, c] - img[r-1, c-1])
            grad[r, c] = gx + gy

    # Precompute Wiener weights from global autocorrelation
    from scipy.signal import correlate2d
    acf = correlate2d(img - img.mean(), img - img.mean(), mode='full')
    acf /= (acf.max() + 1e-10)
    cr2, cc2 = acf.shape[0]//2, acf.shape[1]//2

    offs = []
    for dr in range(-2, 1):
        for dc in range(-2, 3):
            if dr == 0 and dc >= 0: continue
            if dr*dr + dc*dc <= 4:
                offs.append((dr, dc))
    nf = len(offs)

    R_w = np.zeros((nf, nf))
    p_w = np.zeros(nf)
    for i, (d1, e1) in enumerate(offs):
        p_w[i] = acf[cr2 - d1, cc2 - e1]
        for j, (d2, e2) in enumerate(offs):
            R_w[i, j] = acf[cr2 + d1 - d2, cc2 + e1 - e2]
    try:
        w_weights = np.linalg.solve(R_w + 1e-6 * np.eye(nf), p_w)
    except:
        w_weights = np.zeros(nf)
        w_weights[0] = 1.0

    # Compute threshold for switching (median gradient)
    grad_flat = grad[1:, 1:].ravel()
    threshold = np.median(grad_flat[grad_flat > 0]) if np.any(grad_flat > 0) else 10

    for r in range(h):
        for c in range(w):
            if grad[r, c] > threshold:
                # EDGE: use MED
                a = int(plane[r,c-1]) if c>0 else 0
                b = int(plane[r-1,c]) if r>0 else 0
                c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
                pred = med(a, b, c2)
            else:
                # SMOOTH: use Wiener
                pred_val = 0.0
                for k, (dr, dc) in enumerate(offs):
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < h and 0 <= nc < w:
                        pred_val += w_weights[k] * img[nr, nc]
                pred = int(np.clip(round(pred_val), 0, 255))

            residuals.append((int(plane[r, c]) - pred) & 0xFF)

    return len(bc(bytes(residuals)))

# ============================================================
# IDEA 3: ZIGZAG + SCORE-SORTED RESIDUAL ENCODING
# ============================================================

def encode_score_sorted(plane):
    """After MED prediction: zigzag-encode residuals (0,-1,1,-2,2,...),
    then sort by magnitude within each row. The sorted sequence compresses
    better because large residuals cluster at the end.
    Store the permutation for each row as overhead."""
    h, w = plane.shape
    sorted_data = bytearray()
    perm_data = bytearray()

    for r in range(h):
        row_res = []
        for c in range(w):
            a = int(plane[r,c-1]) if c>0 else 0
            b = int(plane[r-1,c]) if r>0 else 0
            c2 = int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res = (int(plane[r,c]) - med(a, b, c2)) & 0xFF
            # Zigzag: unsigned to magnitude order
            signed = res if res <= 128 else res - 256
            zz = 2*abs(signed) - (1 if signed > 0 else 0) if signed != 0 else 0
            row_res.append(zz)

        # Sort by magnitude, keep permutation
        indexed = sorted(enumerate(row_res), key=lambda x: x[1])
        perm = [idx for idx, _ in indexed]
        sorted_vals = [val for _, val in indexed]

        # Encode permutation as delta-coded indices
        for p in perm:
            perm_data.append(p & 0xFF)
        for v in sorted_vals:
            sorted_data.append(v & 0xFF)

    total = len(bc(bytes(sorted_data))) + len(bc(bytes(perm_data)))
    return total

# ============================================================
# IDEA 4: LPC + WIENER HYBRID (use LPC globally, Wiener locally)
# ============================================================

def lpc_coeffs(sig, order):
    n=len(sig)
    if n<order+1: return np.zeros(order)
    r2=np.correlate(sig,sig,mode='full'); r2=r2[n-1:n+order]
    if r2[0]<1e-10: return np.zeros(order)
    a=np.zeros(order); e=r2[0]
    for i in range(order):
        if e<1e-10: break
        lam=r2[i+1]
        for j in range(i): lam-=a[j]*r2[i-j]
        k=lam/e; a_new=np.zeros(order)
        a_new[i]=k
        for j in range(i): a_new[j]=a[j]-k*a[i-1-j]
        a=a_new; e*=(1-k*k)
    return a

def encode_lpc_wiener_hybrid(plane):
    """Use LPC for global trend, Wiener for local deviation from LPC."""
    h, w = plane.shape
    signal = plane.ravel().astype(float)

    # Step 1: LPC order 4 for global prediction
    coeffs = lpc_coeffs(signal - signal.mean(), 4)
    lpc_pred = np.zeros(len(signal))
    for i in range(len(signal)):
        for j in range(min(i, 4)):
            lpc_pred[i] += coeffs[j] * signal[i-1-j]
    lpc_pred = np.clip(np.round(lpc_pred), 0, 255)

    # Step 2: Compute LPC residual
    lpc_residual = signal - lpc_pred  # float residual

    # Step 3: Wiener prediction on the residual (local correction)
    # Reshape back to 2D
    res_2d = lpc_residual.reshape(h, w)

    # Simple 2D Wiener correction: predict residual from neighbors' residuals
    correction = np.zeros((h, w))
    for r in range(1, h):
        for c in range(1, w):
            # Average of 3 known neighbor residuals
            neighbors = [res_2d[r, c-1], res_2d[r-1, c], res_2d[r-1, c-1]]
            correction[r, c] = sum(neighbors) / 3

    final_pred = lpc_pred.reshape(h, w) + correction
    final_pred = np.clip(np.round(final_pred), 0, 255).astype(int)

    final_res = ((plane.astype(int) - final_pred) & 0xFF).astype(np.uint8)

    coeff_bytes = np.array(coeffs, dtype=np.float32).tobytes()
    total = len(coeff_bytes) + len(bc(final_res.tobytes()))
    return total

# ============================================================
# IDEA 5: MULTI-SCALE RESIDUAL (encode at 3 scales simultaneously)
# ============================================================

def encode_multiscale(plane):
    """Encode image at 3 scales:
    Scale 0: 4x-downsampled image (captures global structure)
    Scale 1: 2x-downsampled residual (captures medium structure)
    Scale 2: Full-resolution residual (captures fine detail)
    Each scale predicted from the previous."""
    h, w = plane.shape
    img = plane.astype(float)

    # Scale 0: 4x downsample (average 4x4 blocks)
    s0h, s0w = h//4, w//4
    s0 = np.zeros((s0h, s0w))
    for r in range(s0h):
        for c in range(s0w):
            s0[r, c] = img[4*r:4*r+4, 4*c:4*c+4].mean()
    s0_uint8 = np.clip(np.round(s0), 0, 255).astype(np.uint8)

    # Scale 1: 2x downsample (average 2x2 blocks), predict from upsampled s0
    s1h, s1w = h//2, w//2
    s0_up = np.repeat(np.repeat(s0, 2, axis=0), 2, axis=1)[:s1h, :s1w]
    s1 = np.zeros((s1h, s1w))
    for r in range(s1h):
        for c in range(s1w):
            s1[r, c] = img[2*r:2*r+2, 2*c:2*c+2].mean()
    s1_res = ((np.clip(np.round(s1), 0, 255).astype(int) -
               np.clip(np.round(s0_up), 0, 255).astype(int)) & 0xFF).astype(np.uint8)

    # Scale 2: full resolution, predict from upsampled s1
    s1_up = np.repeat(np.repeat(np.clip(np.round(s1), 0, 255), 2, axis=0), 2, axis=1)[:h, :w]
    s2_res = ((plane.astype(int) - s1_up.astype(int)) & 0xFF).astype(np.uint8)

    total = len(bc(s0_uint8.tobytes())) + len(bc(s1_res.tobytes())) + len(bc(s2_res.tobytes()))
    return total

# ============================================================
# BASELINE METHODS
# ============================================================

def m_lpc4(plane):
    h,w=plane.shape; sig=plane.ravel().astype(float)
    coeffs=lpc_coeffs(sig-sig.mean(),4)
    res=np.zeros(len(sig),dtype=np.uint8)
    for i in range(len(sig)):
        pred=0.0
        for j in range(min(i,4)): pred+=coeffs[j]*sig[i-1-j]
        res[i]=(int(sig[i])-int(np.clip(round(pred),0,255)))&0xFF
    cb=np.array(coeffs,dtype=np.float32).tobytes()
    return len(cb)+len(bc(bytes(res)))

def m_wiener(plane, radius=3):
    h,w=plane.shape; img=plane.astype(float)
    from scipy.signal import correlate2d
    acf=correlate2d(img-img.mean(),img-img.mean(),mode='full')
    acf/=(acf.max()+1e-10)
    cr2,cc2=acf.shape[0]//2,acf.shape[1]//2
    offs=[(dr,dc) for dr in range(-radius,1) for dc in range(-radius,radius+1)
          if not(dr==0 and dc>=0) and dr*dr+dc*dc<=radius*radius]
    nf=len(offs)
    R=np.zeros((nf,nf)); pv=np.zeros(nf)
    for i,(d1,e1) in enumerate(offs):
        pv[i]=acf[cr2-d1,cc2-e1]
        for j,(d2,e2) in enumerate(offs):
            R[i,j]=acf[cr2+d1-d2,cc2+e1-e2]
    try: wts=np.linalg.solve(R+1e-8*np.eye(nf),pv)
    except: wts=np.zeros(nf); wts[0]=1.0
    res=bytearray()
    for r in range(h):
        for c in range(w):
            pred=sum(wts[k]*img[r+dr,c+dc] for k,(dr,dc) in enumerate(offs)
                     if 0<=r+dr<h and 0<=c+dc<w)
            res.append((int(plane[r,c])-int(np.clip(round(pred),0,255)))&0xFF)
    return len(np.array(wts,dtype=np.float32).tobytes())+len(bc(bytes(res)))

def m_colsort(plane):
    h,w=plane.shape
    idx=sorted(range(w),key=lambda c:plane[:,c].tobytes())
    sp=np.stack([plane[:,c] for c in idx],axis=1)
    delta=np.zeros_like(sp); delta[:,0]=sp[:,0]
    delta[:,1:]=((sp[:,1:].astype(int)-sp[:,:-1].astype(int))&0xFF).astype(np.uint8)
    perm=np.array(idx,dtype=np.uint8)
    return len(bc(perm.tobytes()))+len(bc(delta.tobytes()))

# ============================================================
# COMPREHENSIVE TEST
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
print("  BREAKTHROUGH: Adaptive Wiener + Context Switch + Score Conditioning")
print("  Merging tournament theory with signal processing")
print("  kind-pasteur-2026-03-25-S10")
print("="*90)

tests=make_tests()

METHODS=[
    ("MED",          lambda p: m_med(p)),
    ("LPC4",         lambda p: m_lpc4(p)),
    ("Wiener",       lambda p: m_wiener(p)),
    ("ColSort",      lambda p: m_colsort(p)),
    ("AdaptWiener",  lambda p: encode_adaptive_wiener(p)),
    ("CtxSwitch",    lambda p: encode_context_switch(p)),
    ("ScoreSort",    lambda p: encode_score_sorted(p)),
    ("LPC+Wiener",   lambda p: encode_lpc_wiener_hybrid(p)),
    ("MultiScale",   lambda p: encode_multiscale(p)),
    ("Raw",          lambda p: len(bc(p.tobytes()))),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname, _ in METHODS:
    print(f" {mname:>11}", end="")
print(f"  {'BEST':>12}")
print("  "+"-"*(20+12*len(METHODS)+16))

for name,img in sorted(tests.items()):
    ps=png_size(img); sizes={}
    for mname,mfunc in METHODS:
        try: sizes[mname]=mfunc(img)+10
        except Exception as e: sizes[mname]=999999
    best=min(sizes,key=sizes.get)
    ratio=ps/sizes[best] if sizes[best]>0 and sizes[best]<999999 else 0
    print(f"  {name:<12} {ps:>6}", end="")
    for mname,_ in METHODS:
        marker="*" if mname==best else " "
        v=sizes[mname]
        print(f" {v if v<999999 else 'ERR':>10}{marker}", end="")
    print(f"  {best:>12} {ratio:.2f}x")

print(f"""
  THE KEY QUESTIONS:
  1. Does Adaptive Wiener beat fixed Wiener? (adapts to local statistics)
  2. Does Context Switch (MED on edges, Wiener on smooth) beat either alone?
  3. Does LPC+Wiener hybrid capture both global and local structure?
  4. Does score-sorted residual encoding help zlib?
  5. Does multi-scale decomposition capture structure missed by prediction?
""")
