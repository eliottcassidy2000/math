#!/usr/bin/env python3
"""
tc_compose_s315.py — Compositional Codec: transforms as a monoid
opus-2026-03-24-S315

THE HASKELL INSIGHT:
  Instead of: if gradient → delta; if text → stride; if photo → paeth
  Think:      compress = minimize (backend . fold transforms)

  A Transform is a pair (encode, decode) forming a MONOID under composition:
    identity: (id, id)
    compose:  (f, f⁻¹) ∘ (g, g⁻¹) = (f∘g, g⁻¹∘f⁻¹)

  The CODEC is: search over the monoid for the element that minimizes output size.

  type Transform = (ByteString -> ByteString, ByteString -> ByteString)
  type Backend   = ByteString -> ByteString  -- zlib, bz2, lzma
  type Codec     = [Transform] -> Backend -> ByteString -> ByteString

  compress data = minimize [backend (fold ts data) | ts <- powerlist transforms,
                                                      backend <- backends]

  The key insight: transforms COMPOSE. delta∘stride is different from stride∘delta.
  The optimal composition depends on the data. We search the composition space.

FOR VIDEO:
  A video frame is img(t). The delta frame is img(t) - img(t-1).
  Temporal prediction IS a transform: delta_t(frame) = frame ⊕ prev_frame.
  It composes with spatial transforms: backend(spatial(delta_t(frame))).

  The full video codec: backend ∘ spatial ∘ temporal ∘ frame

FOR IMAGES:
  backend ∘ filter_chain ∘ channel_transform ∘ pixel_data

  where filter_chain = foldl (∘) id [f1, f2, f3, ...]
  and we search for the best chain.
"""

import sys, os, struct, zlib, bz2, lzma, time, io
import numpy as np
from PIL import Image
from itertools import product as cartesian
from functools import reduce

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# THE TRANSFORM MONOID
# ============================================================================

class T:
    """A reversible transform: (encode, decode, name)."""
    __slots__ = ['enc', 'dec', 'name']

    def __init__(self, enc, dec, name):
        self.enc = enc
        self.dec = dec
        self.name = name

    def __call__(self, data):
        return self.enc(data)

    def inverse(self, data):
        return self.dec(data)

    def __mul__(self, other):
        """Compose: (self * other)(x) = self(other(x))."""
        return T(
            lambda d, f=self.enc, g=other.enc: f(g(d)),
            lambda d, f=self.dec, g=other.dec: g(f(d)),
            f"{self.name}∘{other.name}"
        )

    def __repr__(self):
        return self.name

# Identity
ID = T(lambda d: d, lambda d: d, "id")

# ============================================================================
# PRIMITIVE TRANSFORMS
# ============================================================================

def _delta_e(d):
    if len(d) == 0: return d
    o = bytearray(len(d)); o[0] = d[0]
    for i in range(1, len(d)): o[i] = (d[i] - d[i-1]) & 0xFF
    return bytes(o)

def _delta_d(d):
    if len(d) == 0: return d
    o = bytearray(len(d)); o[0] = d[0]
    for i in range(1, len(d)): o[i] = (o[i-1] + d[i]) & 0xFF
    return bytes(o)

def _xor_e(d):
    if len(d) == 0: return d
    o = bytearray(len(d)); o[0] = d[0]
    for i in range(1, len(d)): o[i] = d[i] ^ d[i-1]
    return bytes(o)

def _xor_d(d):
    if len(d) == 0: return d
    o = bytearray(len(d)); o[0] = d[0]
    for i in range(1, len(d)): o[i] = o[i-1] ^ d[i]
    return bytes(o)

def _gray_e(d): return bytes(b ^ (b >> 1) for b in d)
def _gray_d(d):
    o = bytearray(len(d))
    for i, g in enumerate(d):
        v = g; m = g >> 1
        while m: v ^= m; m >>= 1
        o[i] = v
    return bytes(o)

def _make_stride(s):
    def enc(d):
        n = len(d)
        if s <= 1 or s >= n: return d
        o = bytearray(n); idx = 0
        for off in range(s):
            for pos in range(off, n, s): o[idx] = d[pos]; idx += 1
        return bytes(o)
    def dec(d):
        n = len(d)
        if s <= 1 or s >= n: return d
        o = bytearray(n); idx = 0
        for off in range(s):
            for pos in range(off, n, s): o[pos] = d[idx]; idx += 1
        return bytes(o)
    return enc, dec

# Primitive transforms
DELTA  = T(_delta_e, _delta_d, "Δ")
XOR    = T(_xor_e, _xor_d, "⊕")
GRAY   = T(_gray_e, _gray_d, "γ")
DELTA2 = DELTA * DELTA  # second-order delta

STRIDES = {}
for s in [2, 3, 4, 8, 16]:
    e, d = _make_stride(s)
    STRIDES[s] = T(e, d, f"S{s}")

# ============================================================================
# IMAGE-SPECIFIC TRANSFORMS
# ============================================================================

def make_paeth_transform(H, W, bpp):
    """Create Paeth filter as a Transform on flat byte arrays."""
    def paeth_pred(a, b, c):
        p = int(a)+int(b)-int(c)
        pa,pb,pc = abs(p-int(a)),abs(p-int(b)),abs(p-int(c))
        if pa<=pb and pa<=pc: return int(a)
        elif pb<=pc: return int(b)
        return int(c)

    row_len = W * bpp

    def enc(data):
        flat = np.frombuffer(data[:H*row_len], dtype=np.uint8).reshape(H, row_len)
        above = np.zeros(row_len, dtype=np.uint8)
        out = bytearray()
        for i in range(H):
            row = flat[i]
            for j in range(row_len):
                left = int(row[j-bpp]) if j >= bpp else 0
                up = int(above[j])
                ul = int(above[j-bpp]) if j >= bpp else 0
                out.append((int(row[j]) - paeth_pred(left, up, ul)) & 0xFF)
            above = row
        return bytes(out)

    def dec(data):
        above = np.zeros(row_len, dtype=np.uint8)
        flat = np.zeros((H, row_len), dtype=np.uint8)
        for i in range(H):
            for j in range(row_len):
                left = int(flat[i, j-bpp]) if j >= bpp else 0
                up = int(above[j])
                ul = int(above[j-bpp]) if j >= bpp else 0
                flat[i,j] = (int(data[i*row_len+j]) + paeth_pred(left, up, ul)) & 0xFF
            above = flat[i].copy()
        return flat.tobytes()

    return T(enc, dec, "P")

def make_sub_transform(H, W, bpp):
    row_len = W * bpp
    def enc(data):
        flat = np.frombuffer(data[:H*row_len], dtype=np.uint8).reshape(H, row_len)
        out = bytearray()
        for i in range(H):
            for j in range(row_len):
                left = int(flat[i, j-bpp]) if j >= bpp else 0
                out.append((int(flat[i,j]) - left) & 0xFF)
        return bytes(out)
    def dec(data):
        flat = np.zeros((H, row_len), dtype=np.uint8)
        for i in range(H):
            for j in range(row_len):
                left = int(flat[i, j-bpp]) if j >= bpp else 0
                flat[i,j] = (int(data[i*row_len+j]) + left) & 0xFF
        return flat.tobytes()
    return T(enc, dec, "Sub")

def make_up_transform(H, W, bpp):
    row_len = W * bpp
    def enc(data):
        flat = np.frombuffer(data[:H*row_len], dtype=np.uint8).reshape(H, row_len)
        above = np.zeros(row_len, dtype=np.uint8)
        out = bytearray()
        for i in range(H):
            for j in range(row_len):
                out.append((int(flat[i,j]) - int(above[j])) & 0xFF)
            above = flat[i]
        return bytes(out)
    def dec(data):
        above = np.zeros(row_len, dtype=np.uint8)
        flat = np.zeros((H, row_len), dtype=np.uint8)
        for i in range(H):
            for j in range(row_len):
                flat[i,j] = (int(data[i*row_len+j]) + int(above[j])) & 0xFF
            above = flat[i].copy()
        return flat.tobytes()
    return T(enc, dec, "Up")

def make_channel_sep(H, W, C):
    """Separate RGB → RRR...GGG...BBB..."""
    if C == 1: return ID
    def enc(data):
        arr = np.frombuffer(data[:H*W*C], dtype=np.uint8).reshape(H, W, C)
        return b''.join(arr[:,:,c].tobytes() for c in range(C))
    def dec(data):
        per_ch = H * W
        channels = [np.frombuffer(data[c*per_ch:(c+1)*per_ch], dtype=np.uint8).reshape(H, W)
                    for c in range(C)]
        return np.stack(channels, axis=-1).tobytes()
    return T(enc, dec, "ChSep")

# ============================================================================
# BACKENDS
# ============================================================================

BACKENDS = [
    ('zlib', lambda d: zlib.compress(d, 9)),
    ('bz2',  lambda d: bz2.compress(d, 9)),
    ('lzma', lambda d: lzma.compress(d)),
]

# ============================================================================
# THE COMPOSITIONAL SEARCH
# ============================================================================

def compose_chain(transforms):
    """Compose a list of transforms left-to-right: t1 then t2 then t3..."""
    if not transforms: return ID
    return reduce(lambda a, b: b * a, transforms)  # b∘a means apply a first, then b

def search_best(data, transforms_pool, max_depth=3):
    """Search for the best composition of up to max_depth transforms + backend.

    Returns (compressed_bytes, chain_name, backend_name).
    """
    best_size = len(data) + 100
    best_result = None
    best_chain = "raw"
    best_backend = "zlib"

    # Depth 0: just backends
    for bname, bfn in BACKENDS:
        try:
            c = bfn(data)
            if len(c) < best_size:
                best_size = len(c)
                best_result = c
                best_chain = "id"
                best_backend = bname
        except: pass

    # Depth 1: single transforms
    for t in transforms_pool:
        try:
            transformed = t(data)
            for bname, bfn in BACKENDS:
                c = bfn(transformed)
                if len(c) < best_size:
                    best_size = len(c)
                    best_result = c
                    best_chain = str(t)
                    best_backend = bname
        except: pass

    # Depth 2: pairs of transforms
    if max_depth >= 2:
        for t1 in transforms_pool:
            for t2 in transforms_pool:
                if t1 is t2: continue
                try:
                    chain = t1 * t2  # apply t2 first, then t1
                    transformed = chain(data)
                    # Only try the cheapest backend (zlib) for depth-2 speed
                    c = zlib.compress(transformed, 9)
                    if len(c) < best_size:
                        best_size = len(c)
                        best_result = c
                        best_chain = str(chain)
                        best_backend = "zlib"
                    # Also try lzma if zlib was close to best
                    if len(c) < best_size * 1.2:
                        c2 = lzma.compress(transformed)
                        if len(c2) < best_size:
                            best_size = len(c2)
                            best_result = c2
                            best_chain = str(chain)
                            best_backend = "lzma"
                except: pass

    return best_result, best_chain, best_backend, best_size

# ============================================================================
# IMAGE COMPRESSOR
# ============================================================================

def compress_image(img_array):
    """Compositional image compression."""
    if img_array.ndim == 2:
        H, W = img_array.shape; C = 1; bpp = 1
    else:
        H, W, C = img_array.shape; bpp = C

    data = img_array.tobytes()

    # Build the transform pool for this image
    pool = [DELTA, XOR, GRAY, DELTA2]

    # Image-specific transforms
    pool.append(make_paeth_transform(H, W, bpp))
    pool.append(make_sub_transform(H, W, bpp))
    pool.append(make_up_transform(H, W, bpp))

    # Stride transforms
    for s in [2, 3, 4, bpp, W*bpp]:
        if 1 < s < len(data):
            pool.append(STRIDES.get(s, T(*_make_stride(s), f"S{s}")))

    # Channel separation (for color)
    if C > 1:
        pool.append(make_channel_sep(H, W, C))

    return search_best(data, pool, max_depth=2)

# ============================================================================
# VIDEO COMPRESSOR
# ============================================================================

def compress_video_frame(frame, prev_frame=None):
    """Compress a video frame, optionally using temporal prediction."""
    if prev_frame is None:
        # Keyframe: just compress the image
        return compress_image(frame)

    # Delta frame: XOR with previous
    delta = frame ^ prev_frame

    # Compress the delta
    result_delta = compress_image(delta)

    # Also try compressing the frame directly (in case delta is worse)
    result_direct = compress_image(frame)

    if result_delta[3] < result_direct[3]:
        return result_delta[0], "Δt∘" + result_delta[1], result_delta[2], result_delta[3]
    return result_direct

# ============================================================================
# BENCHMARK
# ============================================================================

def gen_img(name, W=256, H=256):
    rng = np.random.RandomState(42)
    if name == 'gradient_h':
        return np.array([[int(j*255/max(W-1,1)) for j in range(W)] for i in range(H)], dtype=np.uint8)
    elif name == 'gradient_d':
        return np.array([[int((i+j)*255/max(H+W-2,1)) for j in range(W)] for i in range(H)], dtype=np.uint8)
    elif name == 'circle':
        cx,cy,r = W//2,H//2,min(W,H)//3
        img = np.zeros((H,W), dtype=np.uint8)
        for i in range(H):
            for j in range(W):
                d = ((i-cy)**2+(j-cx)**2)**0.5
                img[i,j] = max(0,min(255,int(255*(1-d/r)))) if d < r else 0
        return img
    elif name == 'checker':
        return np.array([[(255 if (i//16+j//16)%2==0 else 0) for j in range(W)] for i in range(H)], dtype=np.uint8)
    elif name == 'smooth':
        x = np.linspace(0,4*np.pi,W); y = np.linspace(0,4*np.pi,H)
        xx,yy = np.meshgrid(x,y)
        return np.clip(128+100*np.sin(xx)*np.cos(yy),0,255).astype(np.uint8)
    elif name == 'photo':
        img = np.zeros((H,W), dtype=np.float64)
        for _ in range(25):
            cx,cy = rng.randint(W),rng.randint(H)
            r = rng.randint(8,max(9,min(W,H)//3))
            yy,xx = np.ogrid[:H,:W]
            img[(xx-cx)**2+(yy-cy)**2<r**2] = rng.randint(256)
        img += rng.normal(0,10,(H,W))
        return np.clip(img,0,255).astype(np.uint8)
    elif name == 'text':
        img = np.ones((H,W),dtype=np.uint8)*245
        for _ in range(H//8):
            y = rng.randint(H)
            for x in range(0,W,rng.randint(3,10)):
                if rng.random()<0.7:
                    h,w = rng.randint(6,14),rng.randint(2,7)
                    img[max(0,y):min(H,y+h),max(0,x):min(W,x+w)] = rng.randint(0,50)
        return img
    elif name == 'noise':
        return rng.randint(0,256,(H,W),dtype=np.uint8)
    elif name == 'mandelbrot':
        img = np.zeros((H,W),dtype=np.uint8)
        for i in range(H):
            for j in range(W):
                c = complex(-2.5+3.5*j/W,-1.25+2.5*i/H)
                z = 0; n = 0
                while abs(z)<2 and n<255: z=z*z+c; n+=1
                img[i,j] = n
        return img
    elif name == 'satellite':
        base = np.zeros((H,W),dtype=np.float64)
        for freq in [1,2,4,8,16]:
            x = np.linspace(0,freq*np.pi,W); y = np.linspace(0,freq*np.pi,H)
            xx,yy = np.meshgrid(x,y)
            base += (20/freq)*np.sin(xx+rng.random()*np.pi)*np.cos(yy+rng.random()*np.pi)
        for _ in range(8):
            y = rng.randint(H)
            base[max(0,y-1):min(H,y+2),:] = 200
        base += rng.normal(0,5,(H,W))
        return np.clip(base+128,0,255).astype(np.uint8)
    elif name == 'medical':
        img = np.zeros((H,W),dtype=np.float64)
        yy,xx = np.ogrid[:H,:W]
        img[((xx-W//2)**2/(W//3)**2+(yy-H//2)**2/(H//3)**2)<1] = 100
        for _ in range(5):
            cy,cx = H//2+rng.randint(-H//4,H//4),W//2+rng.randint(-W//4,W//4)
            img[max(0,cy-30):min(H,cy+30),max(0,cx-3):min(W,cx+3)] = 220
        img += rng.normal(0,5,(H,W))
        return np.clip(img,0,255).astype(np.uint8)
    return np.zeros((H,W),dtype=np.uint8)


print("=" * 100)
print("  TC COMPOSE — Compositional Codec: transforms as a monoid")
print("  opus-2026-03-24-S315")
print("=" * 100)

# IMAGE BENCHMARK
types = ['gradient_h', 'gradient_d', 'circle', 'checker', 'smooth',
         'photo', 'text', 'noise', 'mandelbrot', 'satellite', 'medical']

print(f"\n  GRAYSCALE 256×256:")
print(f"  {'Image':>12} {'Raw':>7} {'TC':>7} {'PNG':>7} {'TC/PNG':>6} {'Chain':>25} {'Backend':>7}")
print("  " + "-" * 85)

wins = 0; total = 0; tc_sum = 0; png_sum = 0

for name in types:
    img = gen_img(name, 256, 256)
    raw = img.size
    tc_result, chain, backend, tc_size = compress_image(img)

    pil = Image.fromarray(img, mode='L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    win = tc_size <= png_size
    if win: wins += 1
    total += 1; tc_sum += tc_size; png_sum += png_size

    print(f"  {name:>12} {raw:>6} {tc_size:>6} {png_size:>6} {tc_size/png_size:>5.2f}x "
          f"{chain:>25} {backend:>7} {'✓' if win else '✗'}")

print(f"\n  Wins: {wins}/{total}  Total: TC={tc_sum:,} PNG={png_sum:,} ({tc_sum/png_sum:.3f})")

# COLOR
print(f"\n  COLOR 256×256 RGB:")
color_wins = 0; color_total = 0
for name in ['gradient_d', 'circle', 'smooth', 'photo', 'checker', 'mandelbrot', 'satellite']:
    gray = gen_img(name, 256, 256)
    g2 = gen_img('smooth' if name != 'smooth' else 'circle', 256, 256)
    g3 = gen_img('circle' if name != 'circle' else 'gradient_h', 256, 256)
    img = np.stack([gray, g2, g3], axis=-1)

    tc_result, chain, backend, tc_size = compress_image(img)
    pil = Image.fromarray(img, mode='RGB')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    win = tc_size <= png_size
    if win: color_wins += 1
    color_total += 1

    print(f"  {name:>12} {img.size:>7} {tc_size:>7} {png_size:>7} {tc_size/png_size:>5.2f}x "
          f"{chain:>25} {'✓' if win else '✗'}")

print(f"  Color wins: {color_wins}/{color_total}")

# VIDEO
print(f"\n  VIDEO (256×256, 10 frames, motion+noise):")
print(f"  {'Frame':>8} {'Raw':>7} {'TC':>7} {'PNG':>7} {'TC/PNG':>6}")

base = gen_img('photo', 256, 256)
rng = np.random.RandomState(42)
prev = None
tc_video_total = 0; png_video_total = 0

for f in range(10):
    if f == 0:
        frame = base
    else:
        sx, sy = rng.randint(-5, 6), rng.randint(-5, 6)
        frame = np.roll(np.roll(base, sx, axis=1), sy, axis=0)
        frame = np.clip(frame.astype(np.int16) + rng.randint(-8, 9, frame.shape), 0, 255).astype(np.uint8)

    # TC
    if prev is None:
        _, chain, backend, tc_size = compress_image(frame)
        label = "key"
    else:
        delta = (frame.astype(np.int16) - prev.astype(np.int16)).astype(np.uint8)
        _, chain_d, backend_d, tc_delta = compress_image(delta)
        _, chain_k, backend_k, tc_key = compress_image(frame)
        if tc_delta < tc_key:
            tc_size = tc_delta; label = f"Δ{f}"
        else:
            tc_size = tc_key; label = f"K{f}"

    # PNG of the delta (or keyframe)
    to_png = delta if f > 0 and tc_size == tc_delta else frame
    pil = Image.fromarray(to_png if f > 0 else frame, mode='L')
    buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
    png_size = len(buf.getvalue())

    tc_video_total += tc_size; png_video_total += png_size
    print(f"  {label:>8} {frame.size:>6} {tc_size:>6} {png_size:>6} {tc_size/png_size:>5.2f}x "
          f"{'✓' if tc_size <= png_size else ''}")
    prev = frame

print(f"  Video total: TC={tc_video_total:,} PNG={png_video_total:,} ({tc_video_total/png_video_total:.3f})")

print(f"\n{'='*100}")
print(f"  THE COMPOSITIONAL ARCHITECTURE")
print(f"{'='*100}")
print("""
  THE MATHEMATICAL STRUCTURE:

  Transform = (encode :: ByteString → ByteString,
               decode :: ByteString → ByteString)

  The set of transforms forms a MONOID under composition:
    identity:    id = (\\x → x, \\x → x)
    composition: (f, f⁻¹) ∘ (g, g⁻¹) = (f∘g, g⁻¹∘f⁻¹)

  Compression = argmin    |backend(chain(data))|
                chain ∈ T*
                backend ∈ B

  where T* = free monoid generated by the primitive transforms,
  and B = {zlib, bz2, lzma}.

  This is NOT a case-switch. It's a SEARCH over a mathematical space.
  The optimal element of the monoid depends on the data.

  WHY COMPOSITION MATTERS:
    Paeth alone: good for smooth images
    Stride alone: good for structured records
    Paeth ∘ Stride: good for color images (Paeth removes spatial
                    correlation, stride groups similar residuals)
    Delta ∘ Gray: good for slowly-changing values (Gray reduces
                  bit-flip count, delta removes trends)

  The composition space is exponential in depth. We search depth 1-2,
  which gives ~100 candidates per image. A C implementation could
  search depth 3-4 (10,000+ candidates) in milliseconds.

  FOR VIDEO: Temporal delta is just another transform in the monoid.
    video_codec = argmin |backend(spatial(temporal(frame)))|
    where temporal ∈ {id, Δt, Δt², ...}
    and spatial ∈ T*
""")

print("DONE.")
