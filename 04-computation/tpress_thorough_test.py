#!/usr/bin/env python3
"""
tpress_thorough_test.py — Thorough real-world testing of tpress v3
kind-pasteur-2026-03-25-S20gu

Test with REAL data from this repository + edge cases.
Find where we win, where we lose, and improve.
"""

import sys
import os
import zlib
import bz2
import lzma
import numpy as np
import time
from collections import Counter

sys.path.insert(0, '04-computation')
from tpress_v3 import compress_adaptive, compress_block, delta_encode, stride_reorder, gray_encode_bytes, xor_encode, mtf_encode, entropy

sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("  THOROUGH REAL-WORLD TESTING OF tpress v3")
print("=" * 70)


def test_data(name, data):
    """Test one data sample against multiple compressors."""
    raw = len(data)
    if raw == 0:
        return

    # Our compressor
    v3_data, v3_method, v3_size = compress_adaptive(data)

    # Industry compressors
    zl1 = len(zlib.compress(data, 1))
    zl9 = len(zlib.compress(data, 9))
    try:
        bz = len(bz2.compress(data, 9))
    except:
        bz = raw
    try:
        lz = len(lzma.compress(data))
    except:
        lz = raw

    # Best industry
    best_ind = min(zl1, zl9, bz, lz)
    best_name = ['zlib1','zlib9','bz2','lzma'][[zl1,zl9,bz,lz].index(best_ind)]

    # Our ratio vs best industry
    ratio = best_ind / len(v3_data) if len(v3_data) > 0 else 0

    ent = entropy(data)

    marker = "WIN" if ratio > 1.005 else "TIE" if ratio > 0.995 else "LOSE"

    print(f"  {name:>25} {raw:>8}B {len(v3_data):>8}B {best_ind:>8}B({best_name:>5}) {ratio:>7.3f}x {ent:>5.2f}b {v3_method:>15} {marker}")

    return ratio


print(f"\n  {'Data':>25} {'Raw':>8} {'v3':>8} {'BestInd':>13} {'v3/ind':>7} {'Ent':>5} {'Method':>15} {'Result'}")
print("  " + "-" * 110)

# ================================================================
# 1. REAL FILES FROM THIS REPO
# ================================================================
print("\n  --- REAL FILES ---")

real_files = [
    ('CLAUDE.md', 'CLAUDE.md'),
    ('SESSION-LOG', '00-navigation/SESSION-LOG.md'),
    ('formalrank.py', '04-computation/formalrank.py'),
    ('tpress_v3.py', '04-computation/tpress_v3.py'),
    ('explorer.html', '03-artifacts/visualizations/tournament-tiling-explorer.html'),
]

for name, path in real_files:
    full = os.path.join('C:/Users/Eliott/Documents/GitHub/math', path)
    if os.path.exists(full):
        with open(full, 'rb') as f:
            data = f.read()
        # Test first 16K and full file
        if len(data) > 16384:
            test_data(f"{name}[:16K]", data[:16384])
        test_data(name, data)

# ================================================================
# 2. SYNTHETIC STRUCTURED DATA
# ================================================================
print("\n  --- STRUCTURED ---")
np.random.seed(42)

test_data("zeros_4K", bytes(4096))
test_data("ones_4K", bytes([0xFF]*4096))
test_data("counter_4K", bytes([i%256 for i in range(4096)]))
test_data("sawtooth_4K", bytes([i%32 for i in range(4096)]))
test_data("sine_4K", (128+100*np.sin(np.linspace(0,50*np.pi,4096))).clip(0,255).astype(np.uint8).tobytes())
test_data("gradient_2d", np.array([[(i+j)%256 for j in range(64)] for i in range(64)], dtype=np.uint8).tobytes())
test_data("blocks_2d", np.array([[255 if (i//8+j//8)%2==0 else 0 for j in range(64)] for i in range(64)], dtype=np.uint8).tobytes())

# ================================================================
# 3. SYNTHETIC UNSTRUCTURED DATA
# ================================================================
print("\n  --- UNSTRUCTURED ---")

test_data("random_4K", np.random.randint(0,256,4096,dtype=np.uint8).tobytes())
test_data("random_64K", np.random.randint(0,256,65536,dtype=np.uint8).tobytes())

# Text variants
test_data("english_rep", (b"the quick brown fox " * 205)[:4096])
test_data("english_var", b"Tournament theory connects combinatorics algebra and geometry through the study of oriented complete graphs called tournaments. Each tournament on n vertices encodes a total ordering preference or dominance relation among n items.")
test_data("python_code", b"def compress(data):\n    best = zlib.compress(data, 9)\n    for stride in [2,4,8]:\n        c = zlib.compress(reorder(data, stride), 9)\n        if len(c) < len(best): best = c\n    return best\n" * 10)
test_data("xml_data", b'<?xml version="1.0"?><root><item id="1" name="test"><value>42</value></item><item id="2" name="demo"><value>99</value></item></root>' * 10)
test_data("csv_numbers", b"1.23,4.56,7.89\n2.34,5.67,8.90\n3.45,6.78,9.01\n" * 100)
test_data("log_entries", b"[2026-03-25T14:23:01Z] INFO: Request processed in 42ms\n" * 60)

# ================================================================
# 4. EDGE CASES
# ================================================================
print("\n  --- EDGE CASES ---")

test_data("single_byte", b'\x42')
test_data("two_bytes", b'\x42\x43')
test_data("all_same_64K", bytes([0x42]*65536))
test_data("alternating", bytes([0x00,0xFF]*2048))
test_data("near_random", bytes([(i*131+17)%256 for i in range(4096)]))
test_data("mostly_zero", bytes([0]*3900 + [i%256 for i in range(196)]))
test_data("high_entropy", bytes(np.random.randint(200,256,4096,dtype=np.uint8)))
test_data("low_entropy", bytes(np.random.choice([0,1,2,3], 4096).astype(np.uint8)))
test_data("bimodal", bytes(np.random.choice([0,255], 4096, p=[0.9,0.1]).astype(np.uint8)))

# Binary-like patterns
test_data("struct_4byte", bytes([0x00,0x00,i%256,(i*7)%256] for i in range(1024) for _ in [0])[:4096])
test_data("ieee_float", np.random.randn(1024).astype(np.float32).tobytes())
test_data("int16_signal", (10000*np.sin(np.linspace(0,20*np.pi,2048))).astype(np.int16).tobytes())

# ================================================================
# 5. MIXED / COMPOUND DATA
# ================================================================
print("\n  --- MIXED ---")

test_data("header+payload", b'\x89PNG\r\n\x1a\n' + bytes(range(256))*16)
test_data("text+binary", b"Hello World! " * 100 + np.random.randint(0,256,2700,dtype=np.uint8).tobytes())
test_data("sparse+dense", bytes([0]*2048) + np.random.randint(0,256,2048,dtype=np.uint8).tobytes())

print("\nDONE.")
print("=" * 70)
