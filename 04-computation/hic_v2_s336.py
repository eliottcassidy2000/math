#!/usr/bin/env python3
"""
hic_v2_s336.py — Hilbert Image Codec v2: Block-Adaptive

Divides image into 16×16 blocks. For each block:
  Try 8 combinations: {Hilbert, raster} × {none, delta, paeth, avg}
  Pick the smallest compressed block.
  Store: 1 byte per block (which combination) + compressed data.

Author: opus-2026-03-24-S336
"""

import sys
import numpy as np
import zlib
import struct
import time
sys.stdout.reconfigure(line_buffering=True)

BLOCK = 16

def hilbert_d2xy(n, d):
    x = y = 0; s = 1
    while s < n:
        rx = 1 if (d & 2) else 0
        ry = 1 if ((d & 1) ^ rx) else 0
        if ry == 0:
            if rx == 1: x, y = s-1-x, s-1-y
            x, y = y, x
        x += s * rx; y += s * ry
        d >>= 2; s <<= 1
    return x, y

def paeth(a, b, c):
    p = a + b - c
    pa, pb, pc = abs(p-a), abs(p-b), abs(p-c)
    if pa <= pb and pa <= pc: return a
    elif pb <= pc: return b
    return c

def scan_block(block, scan_type):
    B = block.shape[0]
    if scan_type == 0:  # raster
        return block.flatten().tolist()
    else:  # Hilbert
        N = 1
        while N < B: N <<= 1
        padded = np.zeros((N, N), dtype=np.uint8)
        padded[:B, :B] = block[:B, :B]
        pixels = []
        for d in range(N*N):
            x, y = hilbert_d2xy(N, d)
            pixels.append(padded[y][x])
        return pixels[:B*B]

def filter_pixels(pixels, filt):
    n = len(pixels)
    out = bytearray(n)
    if filt == 0:  # none
        for i in range(n): out[i] = pixels[i]
    elif filt == 1:  # delta
        out[0] = pixels[0]
        for i in range(1, n):
            out[i] = (pixels[i] - pixels[i-1]) % 256
    elif filt == 2:  # paeth (stride 4)
        stride = 4
        out[0] = pixels[0]
        for i in range(1, n):
            a = pixels[i-1]
            b = pixels[i-stride] if i >= stride else 0
            c = pixels[i-stride-1] if i >= stride+1 else 0
            out[i] = (pixels[i] - paeth(a, b, c)) % 256
    elif filt == 3:  # avg
        out[0] = pixels[0]
        for i in range(1, n):
            a = pixels[i-1]
            b = pixels[i-4] if i >= 4 else 0
            out[i] = (pixels[i] - (a+b)//2) % 256
    return bytes(out)

def compress_block(block):
    """Try all scan × filter combos, return smallest."""
    best = None; best_mode = 0
    for scan in [0, 1]:  # raster, hilbert
        pixels = scan_block(block, scan)
        for filt in [0, 1, 2, 3]:
            filtered = filter_pixels(pixels, filt)
            mode = scan * 4 + filt
            compressed = zlib.compress(filtered, 9)
            if best is None or len(compressed) < len(best):
                best = compressed
                best_mode = mode
    return best_mode, best

def compress_image(image):
    """Compress 8-bit grayscale image using block-adaptive Hilbert codec."""
    H, W = image.shape
    bh = (H + BLOCK - 1) // BLOCK
    bw = (W + BLOCK - 1) // BLOCK

    # Pad image
    padded = np.zeros((bh * BLOCK, bw * BLOCK), dtype=np.uint8)
    padded[:H, :W] = image

    modes = bytearray()
    block_data = []

    for br in range(bh):
        for bc in range(bw):
            block = padded[br*BLOCK:(br+1)*BLOCK, bc*BLOCK:(bc+1)*BLOCK]
            mode, compressed = compress_block(block)
            modes.append(mode)
            block_data.append(compressed)

    # Pack: header + modes + block lengths + block data
    header = struct.pack('<HHH', H, W, BLOCK)
    n_blocks = bh * bw
    lengths = struct.pack('<' + 'H' * n_blocks, *[len(d) for d in block_data])
    all_blocks = b''.join(block_data)

    # Final zlib on the whole thing for any remaining redundancy
    payload = bytes(modes) + lengths + all_blocks
    return header + zlib.compress(payload, 9)

def decompress_image(data):
    """Decompress block-adaptive Hilbert codec."""
    H, W, B = struct.unpack_from('<HHH', data, 0)
    payload = zlib.decompress(data[6:])

    bh = (H + B - 1) // B
    bw = (W + B - 1) // B
    n_blocks = bh * bw

    modes = list(payload[:n_blocks])
    offset = n_blocks
    lengths = struct.unpack_from('<' + 'H' * n_blocks, payload, offset)
    offset += 2 * n_blocks

    image = np.zeros((bh * B, bw * B), dtype=np.uint8)

    for idx in range(n_blocks):
        br = idx // bw
        bc = idx % bw
        mode = modes[idx]
        scan = mode // 4
        filt = mode % 4
        block_compressed = payload[offset:offset + lengths[idx]]
        offset += lengths[idx]

        filtered = list(zlib.decompress(block_compressed))
        n = B * B

        # Reverse filter
        pixels = list(filtered)
        if filt == 1:
            for i in range(1, n):
                pixels[i] = (pixels[i-1] + pixels[i]) % 256
        elif filt == 2:
            stride = 4
            for i in range(1, n):
                a = pixels[i-1]
                b = pixels[i-stride] if i >= stride else 0
                c = pixels[i-stride-1] if i >= stride+1 else 0
                pixels[i] = (paeth(a, b, c) + pixels[i]) % 256
        elif filt == 3:
            for i in range(1, n):
                a = pixels[i-1]
                b = pixels[i-4] if i >= 4 else 0
                pixels[i] = ((a+b)//2 + pixels[i]) % 256

        # Reverse scan
        block = np.zeros((B, B), dtype=np.uint8)
        if scan == 0:  # raster
            block = np.array(pixels, dtype=np.uint8).reshape(B, B)
        else:  # hilbert
            N = 1
            while N < B: N <<= 1
            for d in range(min(N*N, n)):
                x, y = hilbert_d2xy(N, d)
                if y < B and x < B:
                    block[y][x] = pixels[d]

        image[br*B:(br+1)*B, bc*B:(bc+1)*B] = block

    return image[:H, :W]

# ============================================================
# BENCHMARK
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 70)
    print("  HIC v2 — Block-Adaptive Hilbert Image Codec")
    print("=" * 70)
    print()

    tests = {}
    for N in [64, 128, 256]:
        tests[f'gradient_{N}'] = np.array(
            [[int(255*(r+c)/(2*N-2)) for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f'circle_{N}'] = np.array(
            [[min(255, int(255*np.exp(-((r-N//2)**2+(c-N//2)**2)/(N**2/8))))
              for c in range(N)] for r in range(N)], dtype=np.uint8)
        tests[f'smooth_{N}'] = np.array(
            [[int(128+60*np.sin(r/5)*np.cos(c/4)) for c in range(N)] for r in range(N)], dtype=np.uint8).clip(0,255)

    tests['noise_128'] = np.random.randint(0, 256, (128, 128), dtype=np.uint8)
    tests['photo_128'] = np.array(
        [[int(128+60*np.sin(r/5)*np.cos(c/4)+np.random.normal(0,15))
          for c in range(128)] for r in range(128)], dtype=np.uint8).clip(0,255)

    print("  %-20s %-8s %-8s %-8s %-8s %-6s %-6s" %
          ("Image", "Raw", "HICv2", "HICv1", "zlib", "v2/z", "v2/v1"))
    print("  " + "-"*65)

    from hilbert_image_codec_s336 import hilbert_compress

    wins_v1 = 0; wins_zlib = 0; total = 0
    for name, img in sorted(tests.items()):
        raw = img.size
        total += 1

        v2 = compress_image(img)
        v2_size = len(v2)
        recovered = decompress_image(v2)
        ok = "✓" if np.array_equal(img, recovered) else "✗"

        v1_data, _ = hilbert_compress(img)
        v1_size = len(v1_data)

        z_size = len(zlib.compress(img.tobytes(), 9))

        vs_z = z_size / v2_size if v2_size > 0 else 0
        vs_v1 = v1_size / v2_size if v2_size > 0 else 0

        if v2_size < z_size: wins_zlib += 1
        if v2_size < v1_size: wins_v1 += 1

        print("  %-20s %-8d %-7d%s %-8d %-8d %-6.2f %-6.2f" %
              (name, raw, v2_size, ok, v1_size, z_size, vs_z, vs_v1))

    print(f"\n  Wins vs zlib: {wins_zlib}/{total}")
    print(f"  Wins vs HICv1: {wins_v1}/{total}")
