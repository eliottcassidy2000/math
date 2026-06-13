#!/usr/bin/env python3
"""
Benchmark TIC-GR4 vs all previous codecs and PNG.
opus-2026-04-02

Generates synthetic + loads real test images.
Compares: PNG (Pillow optimize), TIC-GR2, TIC-GR3, TIC-GR4.
"""
import sys, io, os, time, subprocess, tempfile, struct
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

DIR = os.path.dirname(os.path.abspath(__file__))
CODECS = {}
for name, exe_name in [("GR2", "tic_gr2_mac"), ("GR3", "tic_gr3_mac"), ("GR4", "tic_gr4")]:
    path = os.path.join(DIR, exe_name)
    if os.path.exists(path):
        CODECS[name] = path

def png_size(img_array):
    """Get PNG file size via Pillow with max compression."""
    pil = Image.fromarray(img_array)
    buf = io.BytesIO()
    pil.save(buf, format='PNG', optimize=True, compress_level=9)
    return buf.tell()

def tic_bench(exe, rgb, w, h, iters=20):
    """Run a TIC codec benchmark, return (compressed_size, roundtrip_ok)."""
    raw_path = os.path.join(tempfile.gettempdir(), "tic_bench_in.raw")
    rgb.tofile(raw_path)
    try:
        r = subprocess.run([exe, "bench", raw_path, str(w), str(h), str(iters)],
                           capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return 0, False
    sz = 0; ok = False
    for line in r.stdout.split('\n'):
        if 'Compressed:' in line:
            sz = int(line.split()[1])
        if 'PASS' in line:
            ok = True
    if r.returncode != 0 and not ok:
        if r.stderr:
            print(f"      STDERR: {r.stderr.strip()[:200]}")
    try: os.remove(raw_path)
    except: pass
    return sz, ok

def make_test_images(w=256, h=256):
    """Generate a suite of test images."""
    imgs = {}
    rng = np.random.default_rng(42)

    # Solid
    imgs["solid"] = np.full((h, w, 3), 128, dtype=np.uint8)

    # Horizontal gradient
    grad = np.linspace(0, 255, w, dtype=np.uint8)
    imgs["grad_h"] = np.stack([np.tile(grad, (h, 1))] * 3, axis=-1)

    # Smooth gradient (2D)
    x = np.linspace(0, 255, w)
    y = np.linspace(0, 255, h)
    xx, yy = np.meshgrid(x, y)
    smooth = ((xx + yy) / 2).clip(0, 255).astype(np.uint8)
    imgs["grad_2d"] = np.stack([smooth, smooth, smooth], axis=-1)

    # Checker 8x8
    ch = np.zeros((h, w), dtype=np.uint8)
    for r in range(h):
        for c in range(w):
            ch[r, c] = 255 if ((r // 8) + (c // 8)) % 2 else 0
    imgs["checker8"] = np.stack([ch, ch, ch], axis=-1)

    # Random noise
    imgs["noise"] = rng.integers(0, 256, (h, w, 3), dtype=np.uint8)

    # Synthetic photo: smooth + noise + edges
    base = np.zeros((h, w, 3), dtype=np.float64)
    for _ in range(5):
        cx, cy = rng.integers(0, w), rng.integers(0, h)
        r_max = rng.integers(20, 100)
        color = rng.integers(50, 200, 3)
        for r in range(h):
            for c in range(w):
                d = np.sqrt((r - cy)**2 + (c - cx)**2)
                if d < r_max:
                    f = 1.0 - d / r_max
                    base[r, c] += color * f
    base = base.clip(0, 255)
    noise = rng.normal(0, 8, (h, w, 3))
    imgs["synth_photo"] = (base + noise).clip(0, 255).astype(np.uint8)

    # Natural-ish: Perlin approximation via low-pass filtered noise
    raw = rng.normal(0, 1, (h, w, 3))
    from scipy.ndimage import gaussian_filter
    for c in range(3):
        raw[:,:,c] = gaussian_filter(raw[:,:,c], sigma=8)
    raw = (raw - raw.min()) / (raw.max() - raw.min()) * 255
    imgs["perlin_rgb"] = raw.astype(np.uint8)

    # Dithered B/W
    dith = (rng.random((h, w)) < 0.3).astype(np.uint8) * 255
    imgs["dithered"] = np.stack([dith]*3, axis=-1)

    return imgs

def load_real_images():
    """Try to load Kodak images and the real photo."""
    imgs = {}
    base = os.path.join(DIR, "..", "..")

    # Kodak
    import glob
    for f in sorted(glob.glob(os.path.join(base, "test_images", "kodim*.png")))[:4]:
        name = os.path.basename(f)[:7]
        pil = Image.open(f).convert('RGB')
        w0, h0 = pil.size
        crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
        imgs[name] = np.array(crop)

    # Real photo
    real = os.path.join(base, "inbox", "vlcsnap-2026-03-10-17h12m34s408.png")
    if os.path.exists(real):
        pil = Image.open(real).convert('RGB')
        w0, h0 = pil.size
        crop = pil.crop((w0//2-128, h0//2-128, w0//2+128, h0//2+128))
        imgs["vlcsnap"] = np.array(crop)

    return imgs


print("=" * 100)
print("  TIC-GR4 BENCHMARK: JPEG-LS-grade context model vs predecessors")
print("  opus-2026-04-02")
print("=" * 100)

print(f"\n  Available codecs: {', '.join(CODECS.keys())}")

# Generate test images
print("\n  Generating synthetic test images (256×256)...")
try:
    imgs = make_test_images()
except ImportError:
    # No scipy — skip perlin
    print("  (scipy not available, generating without gaussian filter)")
    imgs = {}
    rng = np.random.default_rng(42)
    h, w = 256, 256
    imgs["solid"] = np.full((h, w, 3), 128, dtype=np.uint8)
    grad = np.linspace(0, 255, w, dtype=np.uint8)
    imgs["grad_h"] = np.stack([np.tile(grad, (h, 1))] * 3, axis=-1)
    imgs["noise"] = rng.integers(0, 256, (h, w, 3), dtype=np.uint8)
    base = np.zeros((h, w, 3), dtype=np.float64)
    for _ in range(5):
        cx, cy = rng.integers(0, w), rng.integers(0, h)
        r_max = rng.integers(20, 100)
        color = rng.integers(50, 200, 3)
        for rr in range(h):
            for cc in range(w):
                d = np.sqrt((rr - cy)**2 + (cc - cx)**2)
                if d < r_max:
                    f = 1.0 - d / r_max
                    base[rr, cc] += color * f
    noise = rng.normal(0, 8, (h, w, 3))
    imgs["synth_photo"] = (base + noise).clip(0, 255).astype(np.uint8)

# Load real images
print("  Loading real test images...")
real_imgs = load_real_images()
imgs.update(real_imgs)

# Run benchmarks
print(f"\n  {'Image':<16} {'Raw':>7} {'PNG':>7}", end="")
for name in CODECS: print(f" {name:>7}", end="")
print(f" {'best':>7} {'b/PNG':>7} {'RT':>4}")
print("  " + "-" * (40 + 8 * len(CODECS)))

totals = {'PNG': 0}
for name in CODECS: totals[name] = 0
wins = {n: 0 for n in CODECS}
all_ok = True
n_tests = 0

for iname in sorted(imgs.keys()):
    img = imgs[iname]
    h, w = img.shape[:2]
    raw_sz = h * w * 3

    ps = png_size(img)
    totals['PNG'] += ps

    rgb_flat = img.astype(np.uint8).flatten()
    sizes = {}; oks = {}
    for cname, exe in CODECS.items():
        sz, ok = tic_bench(exe, rgb_flat, w, h)
        sizes[cname] = sz
        oks[cname] = ok
        totals[cname] += sz
        if not ok: all_ok = False

    best_name = min(sizes, key=lambda k: sizes[k]) if sizes else ""
    best_sz = sizes[best_name] if best_name else ps
    ratio = ps / best_sz if best_sz > 0 else 0
    rt = "OK" if all(oks.values()) else "FAIL"
    if best_name: wins[best_name] += 1
    n_tests += 1

    print(f"  {iname:<16} {raw_sz:>7} {ps:>7}", end="")
    for cname in CODECS:
        flag = " " if sizes.get(cname, 0) == 0 else ("*" if cname == best_name else " ")
        print(f" {sizes.get(cname, 0):>6}{flag}", end="")
    print(f" {best_sz:>7} {ratio:>6.3f}x {rt:>4}")

# Aggregates
print(f"\n  {'AGGREGATE':<16} {'':>7} {totals['PNG']:>7}", end="")
for name in CODECS: print(f" {totals[name]:>7}", end="")
best_total = min(totals[name] for name in CODECS if totals[name] > 0)
agg_ratio = totals['PNG'] / best_total if best_total > 0 else 0
print(f" {best_total:>7} {agg_ratio:>6.3f}x {'OK' if all_ok else 'FAIL':>4}")

print(f"\n  Per-codec vs PNG:")
for name in sorted(CODECS.keys()):
    if totals[name] > 0:
        r = totals['PNG'] / totals[name]
        print(f"    {name}: {r:.3f}x ({(r-1)*100:+.1f}%)")

print(f"\n  Win counts (smallest on each image):")
for name in sorted(CODECS.keys()):
    print(f"    {name}: {wins[name]}/{n_tests}")

print(f"\n  Roundtrip: {'ALL PASS' if all_ok else '** FAILURES **'}")
