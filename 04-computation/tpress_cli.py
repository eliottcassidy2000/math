#!/usr/bin/env python3
"""
tpress_cli.py -- Tournament-Powered Compressor: Production CLI
kind-pasteur-2026-03-24-S20cq

USAGE:
  python tpress_cli.py compress input.dat              # -> input.dat.tp
  python tpress_cli.py compress input.dat -o out.tp    # custom output
  python tpress_cli.py decompress input.dat.tp         # -> input.dat
  python tpress_cli.py decompress input.dat.tp -o out  # custom output
  python tpress_cli.py analyze input.dat               # show best pipeline
  python tpress_cli.py bench                            # run full benchmark

FILE FORMAT (.tp):
  Magic:    4 bytes  "TP50"
  Version:  1 byte   0x05
  Flags:    1 byte   (bit 0: has original size)
  OrigSize: 8 bytes  (uint64 LE, original uncompressed size)
  TxID:     1 byte   (transform chain ID)
  BkID:     1 byte   (backend ID)
  Payload:  rest     (compressed data)

Total header: 16 bytes. Decompression needs TxID + BkID to reverse.
"""

import sys
import zlib
import bz2
import lzma
import struct
import time
import os
import argparse
from collections import Counter

__version__ = "5.0.0"
MAGIC = b"TP50"

# ============================================================================
# TRANSFORMS (forward only -- inverse needed for decompress)
# ============================================================================

def t_identity(data): return data
def t_identity_inv(data): return data

def t_delta(data):
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (data[i] - data[i-1]) & 0xFF
    return bytes(out)

def t_delta_inv(data):
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (out[i-1] + data[i]) & 0xFF
    return bytes(out)

def t_delta2(data): return t_delta(t_delta(data))
def t_delta2_inv(data): return t_delta_inv(t_delta_inv(data))

def t_xor(data):
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = data[i] ^ data[i-1]
    return bytes(out)

def t_xor_inv(data):
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = data[i] ^ out[i-1]
    return bytes(out)

def t_gray(data):
    return bytes(b ^ (b >> 1) for b in data)

def t_gray_inv(data):
    def inv_gray(b):
        mask = b >> 1
        while mask:
            b ^= mask
            mask >>= 1
        return b
    return bytes(inv_gray(b) for b in data)

def t_ungray(data): return t_gray_inv(data)
def t_ungray_inv(data): return t_gray(data)

def _stride_fwd(data, s):
    n = len(data)
    if s <= 1 or s >= n: return data
    out = bytearray(n)
    idx = 0
    for offset in range(s):
        for pos in range(offset, n, s):
            out[idx] = data[pos]
            idx += 1
    return bytes(out)

def _stride_inv(data, s):
    n = len(data)
    if s <= 1 or s >= n: return data
    out = bytearray(n)
    idx = 0
    for offset in range(s):
        for pos in range(offset, n, s):
            out[pos] = data[idx]
            idx += 1
    return bytes(out)

def s2(d): return _stride_fwd(d, 2)
def s3(d): return _stride_fwd(d, 3)
def s4(d): return _stride_fwd(d, 4)
def s8(d): return _stride_fwd(d, 8)
def s16(d): return _stride_fwd(d, 16)
def s2_inv(d): return _stride_inv(d, 2)
def s3_inv(d): return _stride_inv(d, 3)
def s4_inv(d): return _stride_inv(d, 4)
def s8_inv(d): return _stride_inv(d, 8)
def s16_inv(d): return _stride_inv(d, 16)

def t_nibble_split(data):
    high = bytes((b >> 4) for b in data)
    low = bytes((b & 0x0F) for b in data)
    return high + low

def t_nibble_split_inv(data):
    n = len(data) // 2
    high = data[:n]
    low = data[n:]
    return bytes(((high[i] << 4) | low[i]) for i in range(n))

def t_rle(data):
    if not data: return data
    counts = Counter(data)
    escape = min(range(256), key=lambda b: counts.get(b, 0))
    out = bytearray([escape])
    i = 0
    n = len(data)
    while i < n:
        run_start = i
        while i < n - 1 and data[i] == data[i+1] and i - run_start < 254:
            i += 1
        run_len = i - run_start + 1
        if run_len >= 4:
            out.append(escape)
            out.append(run_len)
            out.append(data[run_start])
            i += 1
        else:
            for j in range(run_start, run_start + run_len):
                if data[j] == escape:
                    out.append(escape)
                    out.append(0)
                else:
                    out.append(data[j])
            i = run_start + run_len
    result = bytes(out)
    return result if len(result) < len(data) else data

def t_rle_inv(data):
    if not data or len(data) < 2: return data
    escape = data[0]
    out = bytearray()
    i = 1
    n = len(data)
    while i < n:
        if data[i] == escape:
            i += 1
            if i >= n: break
            count = data[i]
            if count == 0:
                out.append(escape)
                i += 1
            else:
                i += 1
                if i >= n: break
                val = data[i]
                out.extend([val] * count)
                i += 1
        else:
            out.append(data[i])
            i += 1
    return bytes(out)

def t_bwt(data):
    n = len(data)
    if n > 16384 or n == 0: return data
    doubled = data + data
    indices = sorted(range(n), key=lambda i: doubled[i:i+n])
    bwt = bytes(doubled[i + n - 1] for i in indices)
    orig_idx = indices.index(0)
    return struct.pack('!H', orig_idx) + bwt

def t_bwt_inv(data):
    if len(data) <= 2: return data
    orig_idx = struct.unpack('!H', data[:2])[0]
    bwt = data[2:]
    n = len(bwt)
    # Standard inverse BWT
    table = sorted(range(n), key=lambda i: bwt[i])
    # Build first column order
    first = sorted(range(n), key=lambda i: bwt[i])
    # LF mapping
    counts = [0] * 256
    lf = [0] * n
    cumul = [0] * 256
    for b in bwt:
        counts[b] += 1
    total = 0
    for i in range(256):
        cumul[i] = total
        total += counts[i]
    counts2 = [0] * 256
    for i in range(n):
        b = bwt[i]
        lf[i] = cumul[b] + counts2[b]
        counts2[b] += 1
    # Reconstruct
    out = bytearray(n)
    idx = orig_idx
    for i in range(n - 1, -1, -1):
        out[i] = bwt[idx]
        idx = lf[idx]
    return bytes(out)

def t_mtf(data):
    if not data: return data
    alphabet = list(range(256))
    out = bytearray(len(data))
    for i, b in enumerate(data):
        idx = alphabet.index(b)
        out[i] = idx
        alphabet.pop(idx)
        alphabet.insert(0, b)
    return bytes(out)

def t_mtf_inv(data):
    if not data: return data
    alphabet = list(range(256))
    out = bytearray(len(data))
    for i, idx in enumerate(data):
        b = alphabet[idx]
        out[i] = b
        alphabet.pop(idx)
        alphabet.insert(0, b)
    return bytes(out)

def t_bwt_mtf(data):
    n = len(data)
    if n > 16384 or n == 0: return data
    bwt_data = t_bwt(data)
    if len(bwt_data) <= 2: return bwt_data
    header = bwt_data[:2]
    return header + t_mtf(bwt_data[2:])

def t_bwt_mtf_inv(data):
    if len(data) <= 2: return data
    header = data[:2]
    payload = t_mtf_inv(data[2:])
    return t_bwt_inv(header + payload)

def t_rev(data): return bytes(reversed(data))
def t_rev_inv(data): return bytes(reversed(data))

def t_sub128(data): return bytes((b - 128) & 0xFF for b in data)
def t_sub128_inv(data): return bytes((b + 128) & 0xFF for b in data)


# ============================================================================
# CHAIN BUILDER with inverse
# ============================================================================

def make_chain(fns_fwd, fns_inv):
    """Return (forward_fn, inverse_fn) for a chain."""
    def fwd(data):
        for fn in fns_fwd:
            data = fn(data)
        return data
    def inv(data):
        for fn in reversed(fns_inv):
            data = fn(data)
        return data
    return fwd, inv


# Transform registry: (name, forward_fn, inverse_fn)
TRANSFORM_REGISTRY = []

def _reg(name, fwd, inv):
    TRANSFORM_REGISTRY.append((name, fwd, inv))

# Singles
_reg('raw',      t_identity,  t_identity_inv)
_reg('delta',    t_delta,     t_delta_inv)
_reg('delta2',   t_delta2,    t_delta2_inv)
_reg('xor',      t_xor,       t_xor_inv)
_reg('gray',     t_gray,      t_gray_inv)
_reg('ungray',   t_ungray,    t_ungray_inv)
_reg('s2',       s2,          s2_inv)
_reg('s3',       s3,          s3_inv)
_reg('s4',       s4,          s4_inv)
_reg('s8',       s8,          s8_inv)
_reg('s16',      s16,         s16_inv)
_reg('bwt',      t_bwt,       t_bwt_inv)
_reg('bwtmtf',   t_bwt_mtf,   t_bwt_mtf_inv)
_reg('mtf',      t_mtf,       t_mtf_inv)
_reg('nibsplit', t_nibble_split, t_nibble_split_inv)
_reg('rle',      t_rle,       t_rle_inv)
_reg('rev',      t_rev,       t_rev_inv)
_reg('sub128',   t_sub128,    t_sub128_inv)

# Chains
def _reg_chain(name, fwd_list, inv_list):
    fwd, inv = make_chain(fwd_list, inv_list)
    TRANSFORM_REGISTRY.append((name, fwd, inv))

_reg_chain('s2+d',      [s2, t_delta],           [s2_inv, t_delta_inv])
_reg_chain('s3+d',      [s3, t_delta],           [s3_inv, t_delta_inv])
_reg_chain('s4+d',      [s4, t_delta],           [s4_inv, t_delta_inv])
_reg_chain('s8+d',      [s8, t_delta],           [s8_inv, t_delta_inv])
_reg_chain('s16+d',     [s16, t_delta],          [s16_inv, t_delta_inv])
_reg_chain('s2+xor',    [s2, t_xor],            [s2_inv, t_xor_inv])
_reg_chain('s3+xor',    [s3, t_xor],            [s3_inv, t_xor_inv])
_reg_chain('s4+xor',    [s4, t_xor],            [s4_inv, t_xor_inv])
_reg_chain('d+gray',    [t_delta, t_gray],       [t_delta_inv, t_gray_inv])
_reg_chain('gray+d',    [t_gray, t_delta],       [t_gray_inv, t_delta_inv])
_reg_chain('gray+s2',   [t_gray, s2],           [t_gray_inv, s2_inv])
_reg_chain('gray+s4',   [t_gray, s4],           [t_gray_inv, s4_inv])
_reg_chain('rle+d',     [t_rle, t_delta],        [t_rle_inv, t_delta_inv])
_reg_chain('nib+d',     [t_nibble_split, t_delta], [t_nibble_split_inv, t_delta_inv])
_reg_chain('nib+s2',    [t_nibble_split, s2],    [t_nibble_split_inv, s2_inv])
_reg_chain('rev+d',     [t_rev, t_delta],        [t_rev_inv, t_delta_inv])
_reg_chain('sub+d',     [t_sub128, t_delta],     [t_sub128_inv, t_delta_inv])
_reg_chain('d+s2',      [t_delta, s2],           [t_delta_inv, s2_inv])
_reg_chain('d+s4',      [t_delta, s4],           [t_delta_inv, s4_inv])
_reg_chain('s2+d+gray', [s2, t_delta, t_gray],   [s2_inv, t_delta_inv, t_gray_inv])
_reg_chain('s4+d+gray', [s4, t_delta, t_gray],   [s4_inv, t_delta_inv, t_gray_inv])
_reg_chain('bwt+d',     [t_bwt, t_delta],        [t_bwt_inv, t_delta_inv])
_reg_chain('mtf+d',     [t_mtf, t_delta],        [t_mtf_inv, t_delta_inv])
_reg_chain('xor+s2',    [t_xor, s2],            [t_xor_inv, s2_inv])
_reg_chain('xor+s4',    [t_xor, s4],            [t_xor_inv, s4_inv])

# Build ID maps
TX_NAME_TO_ID = {entry[0]: i for i, entry in enumerate(TRANSFORM_REGISTRY)}
TX_ID_TO_ENTRY = {i: entry for i, entry in enumerate(TRANSFORM_REGISTRY)}


# ============================================================================
# BACKENDS with inverse
# ============================================================================

BACKEND_REGISTRY = [
    ('zlib1', lambda d: zlib.compress(d, 1),          lambda d: zlib.decompress(d)),
    ('zlib6', lambda d: zlib.compress(d, 6),          lambda d: zlib.decompress(d)),
    ('zlib9', lambda d: zlib.compress(d, 9),          lambda d: zlib.decompress(d)),
    ('bz2',   lambda d: bz2.compress(d, 9),           lambda d: bz2.decompress(d)),
    ('lzma',  lambda d: lzma.compress(d, preset=6),   lambda d: lzma.decompress(d)),
    ('lzma9', lambda d: lzma.compress(d, preset=9|lzma.PRESET_EXTREME), lambda d: lzma.decompress(d)),
]

BK_NAME_TO_ID = {entry[0]: i for i, entry in enumerate(BACKEND_REGISTRY)}
BK_ID_TO_ENTRY = {i: entry for i, entry in enumerate(BACKEND_REGISTRY)}


# ============================================================================
# COMPRESS ENGINE
# ============================================================================

def compress_best(data):
    """Try all transform x backend combinations, return (compressed, tx_id, bk_id)."""
    if len(data) <= 2:
        return data, 0, 0  # raw + zlib1

    best_size = len(data) + 100
    best_data = data
    best_tx = 0
    best_bk = 0

    transforms = TRANSFORM_REGISTRY
    # Skip BWT-based transforms for large data
    if len(data) > 32768:
        transforms = [(n, f, i) for n, f, i in transforms
                      if 'bwt' not in n and 'mtf' not in n]

    for tx_id, (tx_name, tx_fwd, tx_inv) in enumerate(transforms):
        try:
            transformed = tx_fwd(data)
        except:
            continue
        if transformed is None or len(transformed) > len(data) * 2:
            continue

        for bk_id, (bk_name, bk_fwd, bk_inv) in enumerate(BACKEND_REGISTRY):
            try:
                compressed = bk_fwd(transformed)
            except:
                continue
            if compressed is not None and len(compressed) < best_size:
                best_size = len(compressed)
                best_data = compressed
                best_tx = tx_id
                best_bk = bk_id

    return best_data, best_tx, best_bk


def pack(data):
    """Compress data and wrap in .tp container."""
    orig_size = len(data)
    compressed, tx_id, bk_id = compress_best(data)

    # Header: MAGIC(4) + VERSION(1) + FLAGS(1) + ORIGSIZE(8) + TXID(1) + BKID(1) = 16 bytes
    header = MAGIC
    header += struct.pack('<B', 5)           # version
    header += struct.pack('<B', 0x01)        # flags: has original size
    header += struct.pack('<Q', orig_size)   # original size
    header += struct.pack('<B', tx_id)       # transform ID
    header += struct.pack('<B', bk_id)       # backend ID

    return header + compressed, tx_id, bk_id


def unpack(data):
    """Decompress .tp container to original data."""
    if len(data) < 16:
        raise ValueError("File too short to be a .tp file")

    magic = data[:4]
    if magic != MAGIC:
        raise ValueError(f"Invalid magic: expected {MAGIC!r}, got {magic!r}")

    version = data[4]
    flags = data[5]
    orig_size = struct.unpack('<Q', data[6:14])[0]
    tx_id = data[14]
    bk_id = data[15]
    payload = data[16:]

    # Decompress backend
    bk_name, bk_fwd, bk_inv = BK_ID_TO_ENTRY[bk_id]
    decompressed = bk_inv(payload)

    # Inverse transform
    tx_name, tx_fwd, tx_inv = TX_ID_TO_ENTRY[tx_id]
    original = tx_inv(decompressed)

    if len(original) != orig_size:
        raise ValueError(f"Size mismatch: expected {orig_size}, got {len(original)}")

    return original


# ============================================================================
# CLI
# ============================================================================

def cmd_compress(args):
    """Compress a file."""
    input_path = args.input
    output_path = args.output or (input_path + '.tp')

    with open(input_path, 'rb') as f:
        data = f.read()

    t0 = time.time()
    compressed, tx_id, bk_id = pack(data)
    elapsed = time.time() - t0

    with open(output_path, 'wb') as f:
        f.write(compressed)

    orig = len(data)
    comp = len(compressed)
    ratio = orig / comp if comp > 0 else float('inf')
    tx_name = TX_ID_TO_ENTRY[tx_id][0]
    bk_name = BK_ID_TO_ENTRY[bk_id][0]

    # Compare with industry
    zl = len(zlib.compress(data, 9))
    try: bz = len(bz2.compress(data, 9))
    except: bz = orig
    try: lz = len(lzma.compress(data))
    except: lz = orig
    best_ind = min(zl, bz, lz)

    print(f"tpress v{__version__}")
    print(f"  Input:     {input_path} ({orig:,} bytes)")
    print(f"  Output:    {output_path} ({comp:,} bytes)")
    print(f"  Ratio:     {ratio:.2f}x ({(1-comp/orig)*100:.1f}% reduction)")
    print(f"  Pipeline:  {tx_name} + {bk_name}")
    print(f"  Time:      {elapsed*1000:.0f}ms")
    print(f"  vs zlib9:  {zl:,}B  vs bz2: {bz:,}B  vs lzma: {lz:,}B  best: {best_ind:,}B")
    if comp - 16 <= best_ind:  # subtract header
        print(f"  RESULT:    BEATS or TIES best industry ({best_ind/(comp-16):.3f}x)")
    else:
        print(f"  RESULT:    {(comp-16)/best_ind:.3f}x of best industry")


def cmd_decompress(args):
    """Decompress a .tp file."""
    input_path = args.input
    if args.output:
        output_path = args.output
    elif input_path.endswith('.tp'):
        output_path = input_path[:-3]
    else:
        output_path = input_path + '.orig'

    with open(input_path, 'rb') as f:
        data = f.read()

    t0 = time.time()
    original = unpack(data)
    elapsed = time.time() - t0

    with open(output_path, 'wb') as f:
        f.write(original)

    print(f"tpress v{__version__}")
    print(f"  Input:     {input_path} ({len(data):,} bytes)")
    print(f"  Output:    {output_path} ({len(original):,} bytes)")
    print(f"  Time:      {elapsed*1000:.0f}ms")


def cmd_analyze(args):
    """Analyze a file and show optimal compression pipeline."""
    input_path = args.input

    with open(input_path, 'rb') as f:
        data = f.read()

    orig = len(data)
    print(f"tpress v{__version__} -- Analyzing {input_path} ({orig:,} bytes)")
    print("=" * 70)

    # Entropy
    counts = Counter(data)
    n = len(data)
    import math
    ent = -sum(c/n * math.log2(c/n) for c in counts.values() if c > 0)
    unique = len(counts)
    print(f"  Entropy:     {ent:.2f} bits/byte ({ent/8*100:.1f}% of max)")
    print(f"  Unique vals: {unique}/256")
    print(f"  Most common: {counts.most_common(3)}")

    # Try all pipelines, show top 10
    results = []
    for tx_id, (tx_name, tx_fwd, tx_inv) in enumerate(TRANSFORM_REGISTRY):
        try:
            transformed = tx_fwd(data)
        except:
            continue
        if transformed is None or len(transformed) > orig * 2:
            continue
        for bk_id, (bk_name, bk_fwd, bk_inv) in enumerate(BACKEND_REGISTRY):
            try:
                compressed = bk_fwd(transformed)
            except:
                continue
            if compressed is not None:
                results.append((len(compressed), tx_name, bk_name))

    results.sort()
    print(f"\n  Top 10 pipelines:")
    print(f"  {'Rank':>4} {'Size':>8} {'Ratio':>7} {'Transform':>15} {'Backend':>8}")
    for i, (size, tx, bk) in enumerate(results[:10]):
        print(f"  {i+1:4d} {size:7d}B {orig/size:6.2f}x {tx:>15} {bk:>8}")

    # Industry comparison
    zl = len(zlib.compress(data, 9))
    try: bz = len(bz2.compress(data, 9))
    except: bz = orig
    try: lz = len(lzma.compress(data))
    except: lz = orig
    print(f"\n  Industry:  zlib9={zl:,}B  bz2={bz:,}B  lzma={lz:,}B")
    best_ind = min(zl, bz, lz)
    if results:
        advantage = best_ind / results[0][0]
        print(f"  tpress:    {results[0][0]:,}B ({results[0][1]}+{results[0][2]})")
        print(f"  Advantage: {advantage:.3f}x {'WINS' if advantage > 1.005 else 'TIE' if advantage > 0.995 else 'LOSES'}")


def cmd_bench(args):
    """Run comprehensive benchmark."""
    import numpy as np
    np.random.seed(42)

    print(f"tpress v{__version__} -- Benchmark")
    print("=" * 100)

    tests = {}
    tests['counter_4K'] = bytes([i % 256 for i in range(4096)])
    tests['sine_4K'] = (128 + 100*np.sin(np.linspace(0, 50*np.pi, 4096))).clip(0, 255).astype(np.uint8).tobytes()
    N = 64
    tests['gradient_2d'] = np.array([[(i+j)%256 for j in range(N)] for i in range(N)], dtype=np.uint8).tobytes()
    tests['quadratic'] = bytes([int((i*i/16) % 256) for i in range(4096)])
    tests['english_4K'] = (b"the quick brown fox jumps over the lazy dog. " * 100)[:4096]
    tests['python_2K'] = (b"def compress(data):\n    return zlib.compress(data, 9)\n" * 40)[:2048]
    tests['json_4K'] = (b'{"id":1,"name":"test","values":[1,2,3]}\n' * 105)[:4096]
    tests['log_4K'] = (b"[2026-03-25] INFO: Request 42ms\n" * 130)[:4096]
    tests['html_4K'] = (b"<html><body><p>Hello world</p><div class='test'>Content</div></body></html>" * 50)[:4096]
    tests['csv_4K'] = (b"name,age,score\nAlice,30,95\nBob,25,87\nCarol,35,92\n" * 100)[:4096]

    for name, path in [('CLAUDE.md', 'CLAUDE.md'), ('formalrank.py', '04-computation/formalrank.py')]:
        full = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), path)
        if os.path.exists(full):
            with open(full, 'rb') as f:
                tests[name] = f.read()

    tests['random_4K'] = np.random.randint(0, 256, 4096, dtype=np.uint8).tobytes()
    tests['zeros_4K'] = bytes(4096)
    tests['binary_exe'] = bytes([(i*7+13)%256 for i in range(4096)])
    tests['low_entropy'] = bytes(np.random.choice([0,1,2,3], 4096).astype(np.uint8))
    tests['mostly_zeros'] = bytes([0]*3800 + list(range(256)) + [0]*40)
    tests['english_64K'] = (b"the quick brown fox jumps over the lazy dog. " * 1500)[:65536]

    hdr = f"  {'Data':>15} {'Raw':>8} {'tp5':>8} {'zlib9':>8} {'bz2':>8} {'lzma':>8} {'Best':>8} {'Adv':>8} {'Pipeline':>20}"
    print(f"\n{hdr}")

    wins = ties = losses = 0
    for name, data in tests.items():
        raw = len(data)
        t0 = time.time()
        packed, tx_id, bk_id = pack(data)
        elapsed = (time.time() - t0) * 1000
        tp_size = len(packed)

        # Verify roundtrip
        roundtrip = unpack(packed)
        assert roundtrip == data, f"ROUNDTRIP FAIL on {name}!"

        zl = len(zlib.compress(data, 9))
        try: bz = len(bz2.compress(data, 9))
        except: bz = raw
        try: lz = len(lzma.compress(data))
        except: lz = raw
        best_ind = min(zl, bz, lz)

        # Compare payload (excluding our 16-byte header) vs raw industry
        payload_size = tp_size - 16
        ratio = best_ind / payload_size if payload_size > 0 else 0

        if ratio > 1.005: wins += 1; tag = "WIN"
        elif ratio < 0.995: losses += 1; tag = "LOSE"
        else: ties += 1; tag = "TIE"

        tx_name = TX_ID_TO_ENTRY[tx_id][0]
        bk_name = BK_ID_TO_ENTRY[bk_id][0]
        print(f"  {name:>15} {raw:7d}B {tp_size:7d}B {zl:7d}B {bz:7d}B {lz:7d}B {best_ind:7d}B {ratio:7.3f}x {tx_name}+{bk_name:>8} {tag}")

    total = wins + ties + losses
    print(f"\n  SCORE: {wins}W {ties}T {losses}L / {total}")
    print(f"  Win rate: {wins/total*100:.0f}%, Never-worse: {(wins+ties)/total*100:.0f}%")
    print(f"  All roundtrips verified OK")


def cmd_verify(args):
    """Verify roundtrip on a file."""
    with open(args.input, 'rb') as f:
        data = f.read()

    packed, tx_id, bk_id = pack(data)
    original = unpack(packed)

    if original == data:
        tx_name = TX_ID_TO_ENTRY[tx_id][0]
        bk_name = BK_ID_TO_ENTRY[bk_id][0]
        print(f"OK: {len(data):,}B -> {len(packed):,}B -> {len(original):,}B "
              f"({tx_name}+{bk_name}, {len(data)/len(packed):.2f}x)")
    else:
        print(f"FAIL: roundtrip mismatch!")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=f'tpress v{__version__} -- Tournament-Powered Compressor',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s compress myfile.bin           # creates myfile.bin.tp
  %(prog)s decompress myfile.bin.tp      # restores myfile.bin
  %(prog)s analyze myfile.bin            # show optimal pipeline
  %(prog)s bench                         # full benchmark
  %(prog)s verify myfile.bin             # test roundtrip
        """)

    sub = parser.add_subparsers(dest='command')

    p_comp = sub.add_parser('compress', aliases=['c'], help='Compress a file')
    p_comp.add_argument('input', help='Input file path')
    p_comp.add_argument('-o', '--output', help='Output file path')

    p_decomp = sub.add_parser('decompress', aliases=['d'], help='Decompress a .tp file')
    p_decomp.add_argument('input', help='Input .tp file path')
    p_decomp.add_argument('-o', '--output', help='Output file path')

    p_analyze = sub.add_parser('analyze', aliases=['a'], help='Analyze file')
    p_analyze.add_argument('input', help='Input file path')

    p_bench = sub.add_parser('bench', aliases=['b'], help='Run benchmark')

    p_verify = sub.add_parser('verify', aliases=['v'], help='Verify roundtrip')
    p_verify.add_argument('input', help='Input file path')

    args = parser.parse_args()

    if args.command in ('compress', 'c'):
        cmd_compress(args)
    elif args.command in ('decompress', 'd'):
        cmd_decompress(args)
    elif args.command in ('analyze', 'a'):
        cmd_analyze(args)
    elif args.command in ('bench', 'b'):
        cmd_bench(args)
    elif args.command in ('verify', 'v'):
        cmd_verify(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
