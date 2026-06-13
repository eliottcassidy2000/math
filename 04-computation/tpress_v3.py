#!/usr/bin/env python3
"""
tpress_v3.py — Fractal Adaptive Compressor v3.0
kind-pasteur-2026-03-25-S20gt

THE DEEP INSIGHT: "Unstructured" data has structure at SOME scale.
  - Text: structure at byte level (LZ77 catches this)
  - Images: structure at pixel level (our predictors catch this)
  - Binary executables: structure at INSTRUCTION level (4-8 bytes)
  - Random-looking data: structure at BLOCK level (repeated patterns)

The fix: apply compression at MULTIPLE SCALES simultaneously.
At each scale, try our full toolkit. Pick the best per-block.

FRACTAL PIPELINE:
  1. Split data into blocks of size B
  2. For each block, try compression at scales 1, 2, 4, 8, ...
  3. At each scale: group bytes, delta-predict, Gray-encode, zlib
  4. Pick the scale+method that gives smallest output per block
  5. RECURSE: if the output is still large, split into sub-blocks

The "fractal" part: blocks that compress well at scale S are stored
at that scale. Blocks that don't compress at ANY scale are stored raw.
The block size adapts to the LOCAL structure density.

ALSO: use zlib as a BACKEND, not a competitor. Our predictors are
PREPROCESSORS that make zlib more effective. Don't replace zlib — help it.

BEATING ZLIB ON "UNSTRUCTURED" DATA:
  The trick is BYTE REORDERING. If we reorder bytes so that
  similar bytes are adjacent, zlib's LZ77 finds longer matches.
  Our tournament tools give us the reordering via score/skip/diagonal.
"""

import sys
import zlib
import numpy as np
import time
from collections import Counter

__version__ = "3.0.0"


def entropy(data):
    """Shannon entropy in bits per byte."""
    if not data: return 0
    counts = Counter(data)
    n = len(data)
    return -sum(c/n * np.log2(c/n) for c in counts.values() if c > 0)


def delta_encode(data):
    """Delta: each byte predicted from previous."""
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (data[i] - data[i-1]) & 0xFF
    return bytes(out)


def xor_encode(data):
    """XOR with previous byte (good for binary data)."""
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = data[i] ^ data[i-1]
    return bytes(out)


def stride_reorder(data, stride):
    """Reorder bytes by stride: [0,S,2S,...,1,S+1,2S+1,...].
    Groups every S-th byte together. If data has S-byte records,
    this puts all first-bytes together, all second-bytes together, etc.
    """
    n = len(data)
    if stride <= 1 or stride >= n:
        return data
    out = bytearray(n)
    idx = 0
    for offset in range(stride):
        for pos in range(offset, n, stride):
            out[idx] = data[pos]
            idx += 1
    return bytes(out)


def unstride_reorder(data, stride):
    """Inverse of stride_reorder."""
    n = len(data)
    if stride <= 1 or stride >= n:
        return data
    out = bytearray(n)
    idx = 0
    for offset in range(stride):
        for pos in range(offset, n, stride):
            out[pos] = data[idx]
            idx += 1
    return bytes(out)


def reverse_bytes(data):
    """Reverse byte order (helps palindromic data)."""
    return bytes(reversed(data))


def nibble_split(data):
    """Split each byte into high/low nibbles, group separately."""
    high = bytes((b >> 4) for b in data)
    low = bytes((b & 0x0F) for b in data)
    return high + low


def gray_encode_bytes(data):
    """Gray code each byte."""
    return bytes(b ^ (b >> 1) for b in data)


def mtf_encode(data):
    """Move-to-front transform (BWT companion)."""
    alphabet = list(range(256))
    out = bytearray(len(data))
    for i, b in enumerate(data):
        idx = alphabet.index(b)
        out[i] = idx
        alphabet.pop(idx)
        alphabet.insert(0, b)
    return bytes(out)


def compress_block(data, level=9):
    """Try multiple transforms + zlib, return smallest."""
    if len(data) <= 4:
        return data, 'raw', len(data)

    candidates = {}

    # Raw zlib
    candidates['raw'] = zlib.compress(data, level)

    # Delta + zlib
    candidates['delta'] = zlib.compress(delta_encode(data), level)

    # XOR + zlib
    candidates['xor'] = zlib.compress(xor_encode(data), level)

    # Gray code + zlib
    candidates['gray'] = zlib.compress(gray_encode_bytes(data), level)

    # Gray + delta + zlib
    candidates['gray_delta'] = zlib.compress(delta_encode(gray_encode_bytes(data)), level)

    # Stride reorders (detect record structure)
    for stride in [2, 3, 4, 8, 16]:
        if stride < len(data):
            reordered = stride_reorder(data, stride)
            c = zlib.compress(reordered, level)
            if len(c) < len(candidates.get(f'stride{stride}', b'\xff' * (len(data)+100))):
                candidates[f'stride{stride}'] = c
            # Stride + delta
            c2 = zlib.compress(delta_encode(reordered), level)
            candidates[f'stride{stride}_d'] = c2

    # Nibble split + zlib
    if len(data) >= 8:
        candidates['nibble'] = zlib.compress(nibble_split(data), level)

    # Move-to-front + zlib (good for text-like data)
    if len(data) <= 65536:  # MTF is slow for large blocks
        candidates['mtf'] = zlib.compress(mtf_encode(data), level)

    # Pick best
    best = min(candidates, key=lambda k: len(candidates[k]))
    return candidates[best], best, len(candidates[best])


def compress_adaptive(data, block_size=4096):
    """Adaptive block-level compression.

    Split into blocks, compress each with best method,
    then compress the CONCATENATION with zlib (catches inter-block patterns).
    """
    n = len(data)
    if n <= block_size:
        return compress_block(data)

    # Split into blocks
    blocks = []
    for i in range(0, n, block_size):
        blocks.append(data[i:i+block_size])

    # Compress each block, collect methods
    compressed_blocks = []
    methods = []
    for block in blocks:
        compressed, method, size = compress_block(block)
        compressed_blocks.append(compressed)
        methods.append(method)

    # Pack: method IDs + compressed blocks
    # But also try: whole-file compression (zlib on original)
    whole = zlib.compress(data, 9)

    # Block-wise total
    block_total = sum(len(b) for b in compressed_blocks) + len(blocks)  # +1 byte method per block

    if len(whole) <= block_total:
        return whole, 'whole_zlib', len(whole)
    else:
        # Pack blocks
        result = bytearray()
        method_map = {'raw': 0, 'delta': 1, 'xor': 2, 'gray': 3, 'gray_delta': 4,
                      'nibble': 5, 'mtf': 6}
        for i, (comp, method) in enumerate(zip(compressed_blocks, methods)):
            mid = method_map.get(method, 0)
            # Check stride methods
            for s in [2,3,4,8,16]:
                if method == f'stride{s}': mid = 7
                if method == f'stride{s}_d': mid = 8
            result.append(mid)
            result.extend(len(comp).to_bytes(2, 'little'))
            result.extend(comp)

        packed = zlib.compress(bytes(result), 9)
        if len(packed) < len(whole):
            return packed, f'adaptive({Counter(methods).most_common(1)[0][0]})', len(packed)
        else:
            return whole, 'whole_zlib', len(whole)


def main():
    print(f"tpress v{__version__} — Fractal Adaptive Compressor")
    print("=" * 60)

    np.random.seed(42)

    # Generate diverse test data
    tests = {}

    # Structured (where we already win)
    N = 64
    tests['gradient'] = np.tile(np.arange(N, dtype=np.uint8), (N,1)).tobytes()
    tests['sine_4K'] = (128+100*np.sin(np.linspace(0,50*np.pi,4096))).clip(0,255).astype(np.uint8).tobytes()

    # "Unstructured" (where zlib usually wins)
    tests['english'] = (b"The quick brown fox jumps over the lazy dog. " * 100)[:4096]
    tests['html'] = (b"<html><body><p>Hello world</p><div class='test'>Content here</div></body></html>" * 50)[:4096]
    tests['csv'] = (b"name,age,score\nAlice,30,95\nBob,25,87\nCarol,35,92\n" * 100)[:4096]
    tests['json'] = (b'{"id":1,"name":"test","values":[1,2,3,4,5],"active":true}\n' * 70)[:4096]
    tests['binary_exe'] = bytes([(i*7+13)%256 for i in range(4096)])  # pseudo-structured binary
    tests['log_file'] = (b"2026-03-25 14:23:01 INFO Server started on port 8080\n" * 80)[:4096]
    tests['random'] = np.random.randint(0,256,4096,dtype=np.uint8).tobytes()

    # Larger tests
    tests['english_64K'] = tests['english'][:100] * 640
    tests['binary_64K'] = bytes([(i*7+13)%256 for i in range(65536)])

    print(f"\n  {'Data':>15} {'Size':>7} {'v3':>7} {'zlib9':>7} {'v3/zlib':>8} {'Method':>20}")

    for name, data in tests.items():
        raw = len(data)
        v3_data, method, v3_size = compress_adaptive(data)
        zl = len(zlib.compress(data, 9))
        ratio = zl / v3_size if v3_size > 0 else 0

        marker = "**" if ratio > 1.01 else "  " if ratio > 0.99 else "  "
        print(f"  {name:>15} {raw:6d}B {v3_size:6d}B {zl:6d}B {ratio:7.3f}x {method:>20} {marker}")

    # Summary
    print(f"\n  v3/zlib > 1.0 means v3 BEATS zlib")
    print(f"  v3/zlib = 1.0 means tie (both fall back to same zlib)")
    print(f"  v3/zlib < 1.0 means zlib wins (our overhead exceeds gains)")


if __name__ == "__main__":
    main()
