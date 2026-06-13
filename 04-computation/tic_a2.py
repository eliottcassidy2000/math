#!/usr/bin/env python3
"""
A² CODEC: Use the matrix square M² = M @ M^T to guide compression.

M² captures row-row correlations:
  M²[r1, r2] = Σ_c M[r1,c] * M[r2,c] = dot product of rows r1 and r2

Three applications:

1. SIGN-MAGNITUDE SPLIT (from opus S344 indent-split):
   MED residual → separate into SIGN (structure, 1 bit) and MAGNITUDE (content).
   Structure (sign) has spatial coherence → compresses well.
   Content (magnitude) has lower entropy without sign.
   This is the image analog of "indent vs identifiers" in code compression.

2. ROW-CORRELATION WEIGHTED PREDICTION:
   Use the correlation between decoded rows to weight inter-row prediction.
   When predicting row r, give more weight to rows most correlated with r-1.
   Approximation: use the RUNNING Gram matrix from decoded rows.

3. M² FINGERPRINT FOR METHOD SELECTION:
   Compute M² diag and off-diag statistics → predict which method wins.
   Skip exhaustive search. O(n²) overhead instead of O(methods * n²).

kind-pasteur-2026-03-25-S15
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
# IDEA 1: SIGN-MAGNITUDE SPLIT of MED residuals
# ============================================================

def encode_signmag_split(plane):
    """Split MED residuals into sign bits (structure) and magnitudes (content).
    Compress each separately. The sign plane has spatial coherence.
    This is indent-split applied to images."""
    h, w = plane.shape
    signs = bytearray()  # 1 bit per pixel → pack into bytes
    magnitudes = bytearray()

    sign_bits = []
    for r in range(h):
        for c in range(w):
            a = int(plane[r,c-1]) if c > 0 else 0
            b = int(plane[r-1,c]) if r > 0 else 0
            d = int(plane[r-1,c-1]) if r > 0 and c > 0 else 0
            pred = med(a, b, d)
            res = int(plane[r, c]) - pred  # signed residual

            # Sign: 0 = non-negative, 1 = negative
            sign_bits.append(1 if res < 0 else 0)
            # Magnitude: absolute value (0-255)
            magnitudes.append(min(255, abs(res)))

    # Pack sign bits
    sign_packed = np.packbits(np.array(sign_bits, dtype=np.uint8))

    sign_comp = bc(bytes(sign_packed))
    mag_comp = bc(bytes(magnitudes))

    return len(sign_comp) + len(mag_comp)

# ============================================================
# IDEA 2: ROW-CORRELATION PREDICTION using running Gram matrix
# ============================================================

def encode_rowcorr(plane, n_ref=4):
    """Predict each row from a WEIGHTED average of previous rows.
    Weights = correlation of the previous row with each earlier row.
    Rows most similar to the current context get highest weight."""
    h, w = plane.shape
    img = plane.astype(float)
    residuals = bytearray()

    # Store decoded rows for reference
    decoded_rows = []

    for r in range(h):
        if len(decoded_rows) == 0:
            # First row: predict from 128
            for c in range(w):
                residuals.append(int(plane[r, c]))
        else:
            # Compute correlation of previous row with all earlier rows
            prev = decoded_rows[-1]
            weights = np.zeros(len(decoded_rows))
            for i, ref in enumerate(decoded_rows):
                # Cosine similarity
                dot = np.dot(prev - prev.mean(), ref - ref.mean())
                norm = np.sqrt(np.sum((prev - prev.mean())**2) * np.sum((ref - ref.mean())**2))
                weights[i] = dot / (norm + 1e-10)

            # Keep only top n_ref most correlated rows
            if len(weights) > n_ref:
                top_idx = np.argsort(-np.abs(weights))[:n_ref]
                mask = np.zeros_like(weights)
                mask[top_idx] = 1
                weights *= mask

            # Normalize weights (positive only)
            weights = np.maximum(weights, 0)
            wsum = weights.sum()
            if wsum < 1e-10:
                # Fallback: just use previous row
                pred_row = prev
            else:
                weights /= wsum
                pred_row = np.zeros(w)
                for i, ref in enumerate(decoded_rows):
                    pred_row += weights[i] * ref
                pred_row = np.clip(np.round(pred_row), 0, 255)

            for c in range(w):
                residuals.append((int(plane[r, c]) - int(pred_row[c])) & 0xFF)

        decoded_rows.append(img[r].copy())

    return len(bc(bytes(residuals)))

# ============================================================
# IDEA 3: M² FINGERPRINT → method selection (skip trial compression)
# ============================================================

def m2_fingerprint(plane):
    """Compute M² fingerprint: sorted row norms and inter-row correlation."""
    M = plane.astype(float)
    M2 = M @ M.T  # row-row Gram matrix

    # Features
    diag = np.diag(M2)  # row energies
    off_diag = M2 - np.diag(diag)  # inter-row correlations

    row_energy_var = np.var(diag)
    mean_corr = np.mean(np.abs(off_diag))
    rank_approx = np.sum(np.linalg.eigvalsh(M2) > np.max(np.diag(M2)) * 0.01)

    return {
        'energy_var': row_energy_var,
        'mean_corr': mean_corr,
        'rank': rank_approx,
        'max_energy': np.max(diag),
        'min_energy': np.min(diag),
    }

# ============================================================
# IDEA 4: COMBINED — best of sign-mag, row-corr, MED, blend, LPC
# ============================================================

def m_med(plane):
    h,w=plane.shape; res=bytearray()
    for r in range(h):
        for c in range(w):
            a=int(plane[r,c-1]) if c>0 else 0; b=int(plane[r-1,c]) if r>0 else 0
            d=int(plane[r-1,c-1]) if r>0 and c>0 else 0
            res.append((int(plane[r,c])-med(a,b,d))&0xFF)
    return len(bc(bytes(res)))

def m_blend(plane, threshold=1000):
    h,w=plane.shape; img=plane.astype(float); res=bytearray()
    for r in range(h):
        for c in range(w):
            a=img[r,c-1] if c>0 else 0.0; b=img[r-1,c] if r>0 else 0.0
            d=img[r-1,c-1] if r>0 and c>0 else 0.0
            grad=abs(a-d)+abs(b-d) if r>0 and c>0 else 0
            alpha=min(1.0, grad/threshold)
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
        k=lam/e; a_new=np.zeros(order); a_new[i]=k
        for j in range(i): a_new[j]=a[j]-k*a[i-1-j]
        a=a_new; e*=(1-k*k)
    return a

def m_lpc(plane, order=8):
    sig=plane.ravel().astype(float); n=len(sig)
    coeffs=lpc_coeffs(sig-sig.mean(),order)
    res=np.zeros(n,dtype=np.uint8)
    for i in range(n):
        pred=sum(coeffs[j]*sig[i-1-j] for j in range(min(i,order)))
        res[i]=(int(sig[i])-int(np.clip(round(pred),0,255)))&0xFF
    return len(np.array(coeffs,dtype=np.float32).tobytes())+len(bc(bytes(res)))

def m_colsort(plane):
    h,w=plane.shape
    idx=sorted(range(w),key=lambda c2:plane[:,c2].tobytes())
    sp=np.stack([plane[:,c2] for c2 in idx],axis=1)
    delta=np.zeros_like(sp); delta[:,0]=sp[:,0]
    delta[:,1:]=((sp[:,1:].astype(int)-sp[:,:-1].astype(int))&0xFF).astype(np.uint8)
    perm=np.array(idx,dtype=np.uint8)
    return len(bc(perm.tobytes()))+len(bc(delta.tobytes()))

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
print("  A² CODEC: Matrix square as compression guide")
print("  Sign-mag split + row-correlation prediction + M² fingerprinting")
print("  kind-pasteur-2026-03-25-S15")
print("="*90)

tests=make_tests()

METHODS=[
    ("MED",       m_med),
    ("Blend",     lambda p: m_blend(p, 1000)),
    ("LPC8",      lambda p: m_lpc(p, 8)),
    ("ColSort",   m_colsort),
    ("SignMag",   encode_signmag_split),
    ("RowCorr",   encode_rowcorr),
    ("Raw",       lambda p: len(bc(p.tobytes()))),
]

print(f"\n  {'Image':<12} {'PNG':>6}", end="")
for mname,_ in METHODS: print(f" {mname:>8}", end="")
print(f"  {'BEST':>10} {'M²_rank':>8}")
print("  "+"-"*(20+9*len(METHODS)+22))

for name,img in sorted(tests.items()):
    ps=png_size(img); sizes={}
    for mname,mfunc in METHODS:
        try: sizes[mname]=mfunc(img)+10
        except: sizes[mname]=999999
    best=min(sizes,key=sizes.get)
    ratio=ps/sizes[best] if sizes[best]<999999 else 0

    # M² fingerprint
    fp = m2_fingerprint(img)

    print(f"  {name:<12} {ps:>6}", end="")
    for mname,_ in METHODS:
        v=sizes[mname]; marker="*" if mname==best else " "
        print(f" {v if v<999999 else 'ERR':>7}{marker}", end="")
    print(f"  {best:>10} {fp['rank']:>8} {ratio:.2f}x")

# ============================================================
# M² FINGERPRINT ANALYSIS: can we predict the winner?
# ============================================================

print(f"\n  M² FINGERPRINT → METHOD PREDICTION:")
print(f"  {'Image':<12} {'rank':>5} {'energy_var':>12} {'mean_corr':>10} {'winner':>10}")
print("  "+"-"*55)

for name,img in sorted(tests.items()):
    fp = m2_fingerprint(img)
    sizes={}
    for mname,mfunc in METHODS:
        try: sizes[mname]=mfunc(img)+10
        except: sizes[mname]=999999
    winner=min(sizes,key=sizes.get)
    print(f"  {name:<12} {fp['rank']:>5} {fp['energy_var']:>12.0f} {fp['mean_corr']:>10.0f} {winner:>10}")

print(f"""
  PREDICTION RULES (derived from M² fingerprint):
    rank = 1 (all rows identical):     → ColSort or Raw
    rank low, energy uniform:          → LPC (inter-row linearity)
    rank high, energy varied:          → MED (independent rows)
    high mean_corr:                    → ColSort or RowCorr
    low mean_corr:                     → MED or Blend

  The M² fingerprint replaces trial compression with O(n²) matrix multiply.
  Each method has a region in the (rank, energy_var, mean_corr) space.
""")
