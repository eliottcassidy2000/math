#!/usr/bin/env python3
"""
TIC WEISSMAN SCORE: Measure compression ratio AND speed against PNG.

Weissman score = r * (T_ref / T) * α
  where r = compression_ratio / reference_ratio
        T_ref = reference time (PNG encode time)
        T = our encode time
        α ≈ 1 (simplification)

A Weissman score > 1 means we beat PNG in the ratio×speed product.

Uses C kernel (tic_fast.dll) for prediction speed.
Python handles zlib compression + file I/O.

kind-pasteur-2026-03-25-S22
"""
import sys, io, struct, zlib, time, os, glob, ctypes
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# Load C kernel
dll_path = os.path.join(os.path.dirname(__file__), "tic_fast.dll")
if not os.path.exists(dll_path):
    dll_path = os.path.join(os.path.dirname(__file__), "tic_fast.so")
lib = ctypes.CDLL(dll_path)

# Set up function signatures
lib.encode_med.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
lib.decode_med.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
lib.rgb_to_grd.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
lib.grd_to_rgb.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
lib.encode_rgb_full.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
lib.decode_rgb_full.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]

def tic_encode_c(rgb_array):
    """Encode RGB image using C kernel + Python zlib."""
    h, w = rgb_array.shape[:2]
    n = h * w
    rgb_bytes = rgb_array.tobytes()

    # C kernel: RGB → G-RG-BG → MED residuals
    out_buf = ctypes.create_string_buffer(n * 3)
    lib.encode_rgb_full(rgb_bytes, out_buf, h, w)

    # Python zlib: compress each plane
    residuals = bytes(out_buf)
    g_res = residuals[:n]
    rg_res = residuals[n:2*n]
    bg_res = residuals[2*n:3*n]

    g_comp = zlib.compress(g_res, 9)
    rg_comp = zlib.compress(rg_res, 9)
    bg_comp = zlib.compress(bg_res, 9)

    # Pack: header + 3 compressed planes
    header = struct.pack('<4sHHB', b'TIC1', w, h, 1)
    body = b''
    for comp in [g_comp, rg_comp, bg_comp]:
        body += struct.pack('<I', len(comp)) + comp

    return header + body

def tic_decode_c(data, w, h):
    """Decode TIC data using C kernel."""
    n = w * h
    pos = 9  # skip header

    planes_raw = []
    for _ in range(3):
        clen = struct.unpack_from('<I', data, pos)[0]; pos += 4
        raw = zlib.decompress(data[pos:pos+clen])
        planes_raw.append(raw)
        pos += clen

    # Concatenate residuals
    all_res = b''.join(planes_raw)
    res_buf = ctypes.create_string_buffer(all_res)
    rgb_buf = ctypes.create_string_buffer(n * 3)

    lib.decode_rgb_full(res_buf, rgb_buf, h, w)

    return np.frombuffer(bytes(rgb_buf), dtype=np.uint8).reshape(h, w, 3)

def png_encode_time(pil_img, repeats=3):
    """Measure PNG encode time."""
    buf = io.BytesIO()
    # Warmup
    pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    sz = buf.tell()
    # Measure
    t0 = time.perf_counter()
    for _ in range(repeats):
        buf = io.BytesIO()
        pil_img.save(buf, format='PNG', optimize=True, compress_level=9)
    t = (time.perf_counter() - t0) / repeats
    return sz, t

def tic_encode_time(rgb_array, repeats=3):
    """Measure TIC encode time."""
    # Warmup
    data = tic_encode_c(rgb_array)
    # Measure
    t0 = time.perf_counter()
    for _ in range(repeats):
        data = tic_encode_c(rgb_array)
    t = (time.perf_counter() - t0) / repeats
    return len(data), t

# ============================================================
# WEISSMAN SCORE
# ============================================================

def weissman_score(our_ratio, ref_ratio, our_time, ref_time):
    """Simplified Weissman score.
    W = (r / r_ref) * log(T_ref) / log(T)
    where r = compression ratio, T = time.
    W > 1 means we're better overall."""
    if our_time <= 0 or ref_time <= 0:
        return 0
    r_ratio = our_ratio / ref_ratio
    # Use the original Weissman formulation
    import math
    if our_time < 1e-10 or ref_time < 1e-10:
        return r_ratio
    W = r_ratio * math.log(ref_time) / math.log(our_time) if our_time != 1 and ref_time != 1 else r_ratio
    return W

# ============================================================
# BENCHMARK
# ============================================================

print("=" * 90)
print("  TIC WEISSMAN SCORE BENCHMARK: C kernel + Python zlib")
print("  kind-pasteur-2026-03-25-S22")
print("=" * 90)

test_files = sorted(glob.glob("test_images/kodim*.png"))
test_files.append("inbox/vlcsnap-2026-03-10-17h12m34s408.png")

print(f"\n  {'Image':<20} {'Size':>8} {'PNG_sz':>8} {'TIC_sz':>8} {'PNG_r':>6} {'TIC_r':>6} "
      f"{'PNG_ms':>7} {'TIC_ms':>7} {'ratio':>6} {'speed':>6} {'RT':>4}")
print("  " + "-" * 100)

total_png_sz = 0; total_tic_sz = 0; total_raw = 0
total_png_t = 0; total_tic_t = 0
all_pass = True

for fpath in test_files:
    name = os.path.basename(fpath)[:18]
    pil = Image.open(fpath).convert('RGB')

    # Crop to 256×256 center
    w0, h0 = pil.size
    cx, cy = w0//2, h0//2
    crop = pil.crop((cx-128, cy-128, cx+128, cy+128))
    rgb = np.array(crop)
    h, w = rgb.shape[:2]
    raw_sz = h * w * 3

    # PNG benchmark
    png_sz, png_t = png_encode_time(crop, repeats=5)

    # TIC benchmark
    tic_sz, tic_t = tic_encode_time(rgb, repeats=5)

    # Roundtrip verify
    decoded = tic_decode_c(tic_encode_c(rgb), w, h)
    rt_ok = np.array_equal(rgb, decoded)
    if not rt_ok: all_pass = False

    png_ratio = raw_sz / png_sz
    tic_ratio = raw_sz / tic_sz
    size_ratio = png_sz / tic_sz  # >1 means we're smaller
    speed_ratio = png_t / tic_t   # >1 means we're faster

    total_png_sz += png_sz; total_tic_sz += tic_sz; total_raw += raw_sz
    total_png_t += png_t; total_tic_t += tic_t

    print(f"  {name:<20} {raw_sz:>8} {png_sz:>8} {tic_sz:>8} {png_ratio:>5.2f}x {tic_ratio:>5.2f}x "
          f"{png_t*1000:>6.1f}ms {tic_t*1000:>6.1f}ms {size_ratio:>5.3f}x {speed_ratio:>5.2f}x {'OK' if rt_ok else 'FAIL':>4}")

# Aggregates
agg_png_r = total_raw / total_png_sz
agg_tic_r = total_raw / total_tic_sz
agg_size = total_png_sz / total_tic_sz
agg_speed = total_png_t / total_tic_t

print(f"\n  {'AGGREGATE':<20} {total_raw:>8} {total_png_sz:>8} {total_tic_sz:>8} "
      f"{agg_png_r:>5.2f}x {agg_tic_r:>5.2f}x "
      f"{total_png_t*1000:>6.1f}ms {total_tic_t*1000:>6.1f}ms "
      f"{agg_size:>5.3f}x {agg_speed:>5.2f}x {'ALL' if all_pass else 'FAIL':>4}")

print(f"""
  SUMMARY:
    Compression: TIC is {agg_size:.1%} of PNG size (= {agg_size:.3f}x smaller)
    Speed:       TIC is {agg_speed:.2f}x {'faster' if agg_speed > 1 else 'slower'} than PNG
    Roundtrip:   {'ALL PASS' if all_pass else 'FAILURES DETECTED'}

  Size ratio > 1: TIC compresses better
  Speed ratio > 1: TIC is faster (C kernel + Python zlib vs PIL's libpng)

  WEISSMAN-LIKE METRIC: compression_gain × speed_gain = {agg_size * agg_speed:.3f}
    > 1 means we beat PNG in the combined ratio×speed product
    = 1 means we trade compression for speed (or vice versa)
    < 1 means PNG wins overall

  NOTE: Our encode is split between C (MED prediction, ~10% of time)
  and Python zlib (compression, ~90% of time). A full C implementation
  with direct zlib calls would be significantly faster.
""")
