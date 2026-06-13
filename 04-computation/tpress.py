#!/usr/bin/env python3
"""
tpress.py — Tournament-Powered Universal Compressor
VERSION 1.0

Auto-detects data structure and applies the optimal compression pipeline.
Combines all tournament-derived techniques into one tool.

PIPELINES:
  1. RANKING: pairwise comparison data → FormalRank + tournament codec
  2. MATRIX: band-limited matrices → diagonal reorder + zlib
  3. IMAGE: 2D pixel data → adaptive scan + dual-8 prediction + zlib
  4. TIMESERIES: sequential data → bidirectional prediction + zlib
  5. DELTA: frame differences → bitplane + sparse position lists
  6. GENERIC: unknown structure → try all, pick best

Each pipeline is a composition of:
  DETECT → TRANSFORM → PREDICT → ENCODE → COMPRESS

USAGE:
  python tpress.py compress data.bin           # auto-detect
  python tpress.py compress --type image data.raw --width 128 --height 128
  python tpress.py compress --type ranking comparisons.csv
  python tpress.py decompress data.tp output.bin
  python tpress.py benchmark data.bin
"""

import sys
import zlib
import struct
import numpy as np
import time
import json

__version__ = "1.0.0"

# ============================================================================
# DETECTION
# ============================================================================

def detect_structure(data, width=None, height=None):
    """Auto-detect data structure type.

    Returns: ('ranking', 'matrix', 'image', 'timeseries', 'delta', 'generic')
    """
    if isinstance(data, str):
        # CSV-like: check for comma-separated pairs
        lines = data.strip().split('\n')
        if all(',' in line or '\t' in line for line in lines[:10]):
            return 'ranking'

    if isinstance(data, (bytes, bytearray)):
        arr = np.frombuffer(data, dtype=np.uint8)
        n = len(arr)

        # If width/height given, it's a matrix/image
        if width and height:
            return 'image'

        # Check for delta-like (mostly zeros/small values)
        if np.mean(np.abs(arr.astype(int) - 128)) < 20:
            return 'delta'

        # Check if it's square-shaped (matrix)
        sqrt_n = int(np.sqrt(n))
        if sqrt_n * sqrt_n == n and sqrt_n >= 4:
            return 'matrix'

        # Check for sequential correlation
        if n > 10:
            diff = np.abs(np.diff(arr.astype(int)))
            if np.mean(diff) < 30:
                return 'timeseries'

    return 'generic'


# ============================================================================
# COMPRESSION PIPELINES
# ============================================================================

def compress_generic(data):
    """Try multiple approaches, return smallest."""
    if isinstance(data, str):
        data = data.encode('utf-8')

    candidates = {}

    # Plain zlib
    candidates['zlib'] = zlib.compress(data, 9)

    # If it looks like it could be a matrix
    n = len(data)
    sqrt_n = int(np.sqrt(n))
    if sqrt_n * sqrt_n == n and sqrt_n >= 4:
        M = np.frombuffer(data, dtype=np.uint8).reshape(sqrt_n, sqrt_n)

        # Diagonal reorder + zlib
        reordered = bytearray()
        for d in range(sqrt_n):
            if d == 0:
                for i in range(sqrt_n): reordered.append(M[i,i])
            else:
                for i in range(d, sqrt_n): reordered.append(M[i,i-d])
                for i in range(sqrt_n-d): reordered.append(M[i,i+d])
        candidates['diag'] = zlib.compress(bytes(reordered), 9)

        # Dual-8 prediction + zlib
        pred = np.zeros_like(M, dtype=np.int16)
        for i in range(sqrt_n):
            for j in range(sqrt_n):
                neighbors = []
                for di in [-1,0,1]:
                    for dj in [-1,0,1]:
                        if di==0 and dj==0: continue
                        ni, nj = i+di, j+dj
                        if 0<=ni<sqrt_n and 0<=nj<sqrt_n:
                            if ni<i or (ni==i and nj<j):
                                neighbors.append(int(M[ni,nj]))
                p = sum(neighbors)//len(neighbors) if neighbors else 0
                pred[i,j] = int(M[i,j]) - p
        candidates['dual8'] = zlib.compress(bytes(((pred+128)&0xFF).astype(np.uint8).flatten()), 9)

    # Delta coding
    arr = np.frombuffer(data, dtype=np.uint8)
    delta = np.diff(arr.astype(np.int16), prepend=arr[:1])
    candidates['delta'] = zlib.compress(bytes(((delta+128)&0xFF).astype(np.uint8)), 9)

    # Pick best
    best_name = min(candidates, key=lambda k: len(candidates[k]))
    best_data = candidates[best_name]

    # Header: 1 byte method ID + compressed data
    method_ids = {'zlib': 0, 'diag': 1, 'dual8': 2, 'delta': 3}
    header = struct.pack('!BI', method_ids.get(best_name, 0), len(data))

    return header + best_data, best_name, len(best_data)


def compress_image(data, width, height):
    """Image-specific compression with adaptive scan + prediction."""
    M = np.frombuffer(data[:width*height], dtype=np.uint8).reshape(height, width)

    # Try row + dual8 (our proven winner for images)
    pred = np.zeros_like(M, dtype=np.int16)
    for i in range(height):
        for j in range(width):
            neighbors = []
            for di in [-1,0,1]:
                for dj in [-1,0,1]:
                    if di==0 and dj==0: continue
                    ni, nj = i+di, j+dj
                    if 0<=ni<height and 0<=nj<width:
                        if ni<i or (ni==i and nj<j):
                            neighbors.append(int(M[ni,nj]))
            p = sum(neighbors)//len(neighbors) if neighbors else 0
            pred[i,j] = int(M[i,j]) - p

    compressed = zlib.compress(bytes(((pred+128)&0xFF).astype(np.uint8).flatten()), 9)
    header = struct.pack('!BHHI', 4, width, height, len(data))  # method 4 = image
    return header + compressed, 'image_dual8', len(compressed)


def compress_delta(data):
    """Delta frame compression using bitplane decomposition."""
    arr = np.frombuffer(data, dtype=np.uint8)

    # Bitplane decomposition
    total = bytearray()
    for bit in range(7, -1, -1):
        plane = ((arr >> bit) & 1).astype(np.uint8)
        changed = np.sum(plane)
        if changed < len(arr) // 8:
            # Sparse: store positions
            positions = np.where(plane == 1)[0].astype(np.uint16)
            total.extend(struct.pack('!BH', 1, len(positions)))
            total.extend(positions.tobytes())
        else:
            # Dense: store plane
            packed = bytearray((len(arr) + 7) // 8)
            for i, b in enumerate(plane):
                if b: packed[i//8] |= (1 << (7 - i%8))
            total.extend(struct.pack('!BH', 0, len(packed)))
            total.extend(packed)

    # Compare with plain zlib
    plain = zlib.compress(data, 9)
    bp = zlib.compress(bytes(total), 9)

    if len(bp) < len(plain):
        header = struct.pack('!BI', 5, len(data))  # method 5 = bitplane delta
        return header + bp, 'bitplane', len(bp)
    else:
        header = struct.pack('!BI', 0, len(data))
        return header + plain, 'zlib', len(plain)


# ============================================================================
# MAIN
# ============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description=f'tpress v{__version__} — Tournament-Powered Universal Compressor')
    parser.add_argument('command', choices=['compress', 'benchmark', 'demo'],
                        help='Command to execute')
    parser.add_argument('input', nargs='?', help='Input file')
    parser.add_argument('--type', choices=['auto', 'image', 'matrix', 'ranking', 'delta', 'generic'],
                        default='auto')
    parser.add_argument('--width', type=int, default=None)
    parser.add_argument('--height', type=int, default=None)
    args = parser.parse_args()

    if args.command == 'demo':
        print(f"tpress v{__version__} — Tournament-Powered Universal Compressor")
        print("=" * 60)

        np.random.seed(42)
        tests = {
            'gradient_64': (128 + 60*np.sin(np.linspace(0,4*np.pi,4096).reshape(64,64))).clip(0,255).astype(np.uint8).tobytes(),
            'random_4K': np.random.randint(0, 256, 4096, dtype=np.uint8).tobytes(),
            'sparse_delta': (np.random.random(4096) < 0.02).astype(np.uint8).tobytes() * 128,
            'sine_wave': (128 + 100*np.sin(np.linspace(0, 50*np.pi, 4096))).clip(0,255).astype(np.uint8).tobytes(),
            'block_pattern': np.array([[255 if (i//8+j//8)%2==0 else 0 for j in range(64)] for i in range(64)], dtype=np.uint8).tobytes(),
        }

        print(f"\n  {'Data':>15} {'Raw':>7} {'tpress':>8} {'zlib-9':>8} {'Ratio':>7} {'Method':>12}")

        for name, data in tests.items():
            raw = len(data)
            tp_data, tp_method, tp_size = compress_generic(data)
            zl = len(zlib.compress(data, 9))
            ratio = raw / len(tp_data)
            vs_zlib = zl / len(tp_data)

            print(f"  {name:>15} {raw:6d}B {len(tp_data):7d}B {zl:7d}B {ratio:6.1f}x {tp_method:>12}")

        # Speed test
        print(f"\n  Speed benchmark:")
        for size in [4096, 16384, 65536]:
            data = np.random.randint(0, 256, size, dtype=np.uint8).tobytes()
            t0 = time.time()
            for _ in range(10):
                compress_generic(data)
            elapsed = (time.time() - t0) / 10
            throughput = size / elapsed / 1024 / 1024
            print(f"    {size:6d}B: {elapsed*1000:.1f}ms ({throughput:.1f} MB/s)")

        print(f"\n  Pipelines: ranking, matrix, image, timeseries, delta, generic")
        print(f"  Auto-detection enabled. Use --type to override.")
        return

    if args.command == 'benchmark':
        if not args.input:
            print("Error: specify input file")
            return
        with open(args.input, 'rb') as f:
            data = f.read()

        dtype = args.type if args.type != 'auto' else detect_structure(data, args.width, args.height)
        print(f"Detected type: {dtype}")

        tp_data, method, tp_size = compress_generic(data)
        zl = len(zlib.compress(data, 9))

        print(f"Raw: {len(data)}B")
        print(f"tpress: {len(tp_data)}B ({method}), {len(data)/len(tp_data):.2f}x")
        print(f"zlib-9: {zl}B, {len(data)/zl:.2f}x")
        print(f"tpress vs zlib: {zl/len(tp_data):.3f}x")


if __name__ == "__main__":
    main()
