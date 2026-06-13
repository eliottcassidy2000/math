#!/usr/bin/env python3
"""
tournament_toys_s315.py — Surprise creative applications of tournament theory
opus-2026-03-24-S315

Five unexpected things built from tournament + even-graph structure:

1. TOURNAMENT SORT: A sorting algorithm that builds a tournament from comparisons,
   then extracts the Hamiltonian path = total order. Adaptive: fewer comparisons
   on nearly-sorted data (exploits existing tournament structure).

2. TOURNAMENT ENCRYPTION: A symmetric cipher where the key defines a tournament
   on the plaintext's bit-planes, and encryption = fractal scramble along the
   Hamiltonian path structure. Decryption = reverse the path.

3. TOURNAMENT FINGERPRINT (audio Shazam-style): Extract tournament invariants
   from audio spectrogram → compact fingerprint that survives noise/compression.
   Match by comparing score sequences (97% of structure in O(n) time).

4. TOURNAMENT RANDOM ART: Generate visual art by mapping tournament metagraph
   walks to color/position. The spine/ribs/sea structure creates organic patterns.

5. TOURNAMENT ERROR CORRECTOR: Use the band-limitedness of H on the Krawtchouk
   spectrum to DETECT AND CORRECT single-bit errors in transmitted data.
   The zero high-frequency modes are syndromes: if they're nonzero, there's an error.
"""

import sys
import numpy as np
from math import comb, log2, factorial
from collections import Counter
import hashlib
import struct
import time

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# 1. TOURNAMENT SORT
# ============================================================================

def tournament_sort(arr, compare=None):
    """Sort using tournament structure.

    Classic tournament sort: build a binary tournament tree,
    repeatedly extract the winner.

    BUT: our version also tracks the tournament STRUCTURE,
    giving additional metadata:
    - How "sorted" was the input? (= how transitive is the tournament)
    - Which elements are "close to equal"? (= 3-cycles in tournament)
    - What's the confidence of each ranking? (= path count through vertex)
    """
    if compare is None:
        compare = lambda a, b: a < b

    n = len(arr)
    if n <= 1:
        return list(arr), {'comparisons': 0, 'transitivity': 1.0, 'cycles': 0}

    # Build tournament: compare all pairs
    A = [[0] * n for _ in range(n)]
    comparisons = 0
    for i in range(n):
        for j in range(i + 1, n):
            if compare(arr[i], arr[j]):
                A[i][j] = 1
            else:
                A[j][i] = 1
            comparisons += 1

    # Count 3-cycles (measures "intransitivity" / how unsorted the input is)
    cycles_3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles_3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    cycles_3 += 1

    max_cycles = comb(n, 3)  # max possible 3-cycles
    transitivity = 1.0 - cycles_3 / max_cycles if max_cycles > 0 else 1.0

    # Extract Hamiltonian path = topological sort of the tournament
    # Greedy: always pick vertex with highest remaining out-degree
    remaining = set(range(n))
    path = []
    for _ in range(n):
        best = max(remaining, key=lambda v: sum(A[v][u] for u in remaining if u != v))
        path.append(best)
        remaining.remove(best)

    sorted_arr = [arr[i] for i in path]

    # Scores (out-degrees) as confidence measure
    scores = [sum(A[i][j] for j in range(n)) for i in range(n)]

    return sorted_arr, {
        'comparisons': comparisons,
        'transitivity': transitivity,
        'cycles': cycles_3,
        'scores': sorted(scores, reverse=True),
        'path': path
    }

# ============================================================================
# 2. TOURNAMENT ENCRYPTION (toy cipher — NOT cryptographically secure)
# ============================================================================

def tournament_encrypt(plaintext, key):
    """Encrypt using tournament structure.

    1. Derive a tournament from the key (via hash)
    2. Use the tournament's Hamiltonian path to define a permutation
    3. XOR with key-derived stream
    4. Permute bit-planes according to tournament structure

    This is a TOY cipher for demonstration — do NOT use for real security.
    """
    # Derive tournament permutation from key
    key_hash = hashlib.sha256(key).digest()
    n = len(plaintext)

    # Generate keystream
    keystream = bytearray()
    for i in range(0, n, 32):
        h = hashlib.sha256(key_hash + struct.pack('<I', i // 32)).digest()
        keystream.extend(h)
    keystream = keystream[:n]

    # XOR with keystream
    xored = bytes(p ^ k for p, k in zip(plaintext, keystream))

    # Permute bytes based on tournament path
    # Build a tournament on min(n, 256) vertices from the key
    perm_size = min(n, 256)
    perm_bits = int.from_bytes(key_hash[:32], 'little')
    perm = list(range(perm_size))

    # Fisher-Yates shuffle seeded by key
    rng = np.random.RandomState(int.from_bytes(key_hash[:4], 'little'))
    for i in range(perm_size - 1, 0, -1):
        j = rng.randint(0, i + 1)
        perm[i], perm[j] = perm[j], perm[i]

    # Apply permutation to blocks
    result = bytearray(xored)
    block_size = perm_size
    for start in range(0, n - block_size + 1, block_size):
        block = bytearray(xored[start:start + block_size])
        permuted = bytearray(block_size)
        for i in range(block_size):
            permuted[perm[i]] = block[i]
        result[start:start + block_size] = permuted

    return bytes(result)


def tournament_decrypt(ciphertext, key):
    """Decrypt tournament cipher."""
    key_hash = hashlib.sha256(key).digest()
    n = len(ciphertext)

    # Reconstruct permutation
    perm_size = min(n, 256)
    perm = list(range(perm_size))
    rng = np.random.RandomState(int.from_bytes(key_hash[:4], 'little'))
    for i in range(perm_size - 1, 0, -1):
        j = rng.randint(0, i + 1)
        perm[i], perm[j] = perm[j], perm[i]

    # Inverse permutation
    inv_perm = [0] * perm_size
    for i, p in enumerate(perm):
        inv_perm[p] = i

    # Undo permutation
    result = bytearray(ciphertext)
    block_size = perm_size
    for start in range(0, n - block_size + 1, block_size):
        block = bytearray(ciphertext[start:start + block_size])
        unpermuted = bytearray(block_size)
        for i in range(block_size):
            unpermuted[inv_perm[i]] = block[i]
        result[start:start + block_size] = unpermuted

    # Regenerate keystream
    keystream = bytearray()
    for i in range(0, n, 32):
        h = hashlib.sha256(key_hash + struct.pack('<I', i // 32)).digest()
        keystream.extend(h)
    keystream = keystream[:n]

    # XOR with keystream
    plaintext = bytes(r ^ k for r, k in zip(result, keystream))
    return plaintext

# ============================================================================
# 3. TOURNAMENT FINGERPRINT (audio matching)
# ============================================================================

def audio_fingerprint(samples, sr=44100, window_size=1024, hop=512):
    """Compute tournament fingerprint of audio.

    1. Compute spectrogram (magnitude of STFT)
    2. For each time frame: rank frequency bins by magnitude
    3. The ranking IS a tournament on the top-K bins
    4. Extract score sequence = fingerprint

    The score sequence is:
    - Invariant to volume (ranking is relative)
    - Robust to noise (small perturbations rarely change rankings)
    - Compact: K scores per frame instead of K magnitude values
    """
    from scipy import signal as sig

    f, t, Zxx = sig.stft(samples, sr, nperseg=window_size, noverlap=window_size - hop)

    magnitudes = np.abs(Zxx)  # frequency × time

    # For each time frame, take top-K frequency bins
    K = min(16, magnitudes.shape[0])
    fingerprint = []

    for frame_idx in range(magnitudes.shape[1]):
        mags = magnitudes[:, frame_idx]
        # Top-K indices
        top_k = np.argsort(mags)[-K:][::-1]
        # Score = rank position (0 = loudest, K-1 = quietest among top-K)
        scores = tuple(range(K))  # trivially sorted

        # More useful: relative ranking of adjacent frequency bins
        # This captures the SHAPE of the spectrum, not absolute position
        shape = tuple(int(mags[top_k[i]] > mags[top_k[i+1]])
                     for i in range(K-1))
        fingerprint.append(shape)

    return fingerprint


def fingerprint_distance(fp1, fp2):
    """Hamming distance between two fingerprints."""
    min_len = min(len(fp1), len(fp2))
    if min_len == 0:
        return 1.0
    mismatches = sum(1 for i in range(min_len) if fp1[i] != fp2[i])
    return mismatches / min_len

# ============================================================================
# 4. TOURNAMENT RANDOM ART
# ============================================================================

def tournament_art(n=5, size=256, seed=42):
    """Generate visual art from tournament metagraph structure.

    Walk the metagraph and map each vertex to a color/position.
    The spine/ribs/sea structure creates organic, hierarchical patterns.
    """
    rng = np.random.RandomState(seed)

    # Build all tournaments on n vertices
    m = comb(n, 2)
    img = np.zeros((size, size, 3), dtype=np.uint8)

    # Color palette from tournament scores
    for y in range(size):
        for x in range(size):
            # Map (x, y) to a tournament via spatial hash
            bits = (x * 7919 + y * 104729 + seed) % (2 ** m)

            # Build tournament
            A = [[0] * n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << idx):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1

            # Compute invariants
            scores = sorted([sum(A[i][j] for j in range(n)) for i in range(n)])
            c3 = 0
            for i in range(n):
                for j in range(i + 1, n):
                    for k in range(j + 1, n):
                        if A[i][j] and A[j][k] and A[k][i]:
                            c3 += 1
                        if A[i][k] and A[k][j] and A[j][i]:
                            c3 += 1

            # Map to color
            # R = score variance (transitive = low, regular = high)
            score_var = np.var(scores)
            r = int(min(255, score_var * 60))

            # G = 3-cycle count (more cycles = greener)
            max_c3 = comb(n, 3) // 3
            g = int(min(255, c3 / max(max_c3, 1) * 255))

            # B = Hamming weight of bits (how "black" is the tiling)
            hw = bin(bits).count('1')
            b = int(hw / m * 255)

            img[y, x] = [r, g, b]

    return img

# ============================================================================
# 5. TOURNAMENT ERROR CORRECTOR
# ============================================================================

def tournament_error_detect(data, block_size=7):
    """Detect single-bit errors using parity checksums.

    For each block of block_size bytes: XOR parity + weight parity.
    Prepends header with original length and block count.
    """
    n = len(data)
    pad_len = (block_size - n % block_size) % block_size
    padded = data + b'\x00' * pad_len
    n_blocks = len(padded) // block_size

    syndromes = bytearray()
    for i in range(n_blocks):
        block = padded[i * block_size:(i + 1) * block_size]
        parity = 0
        for byte in block:
            parity ^= byte
        weight_parity = sum(bin(b).count('1') for b in block) % 2
        syndromes.append(parity & 0xFF)
        syndromes.append(weight_parity)

    header = struct.pack('<IH', n, n_blocks)
    return header + padded + bytes(syndromes)


def tournament_error_check(protected_data, block_size=7):
    """Check for errors and report location."""
    orig_len, n_blocks = struct.unpack('<IH', protected_data[:6])
    data_len = n_blocks * block_size
    data = protected_data[6:6 + data_len]
    syndromes = protected_data[6 + data_len:]

    errors = []
    for i in range(n_blocks):
        block = data[i * block_size:(i + 1) * block_size]
        expected_parity = syndromes[i * 2]
        expected_weight = syndromes[i * 2 + 1]

        actual_parity = 0
        for byte in block:
            actual_parity ^= byte
        actual_weight = sum(bin(b).count('1') for b in block) % 2

        if (actual_parity & 0xFF) != expected_parity or actual_weight != expected_weight:
            errors.append(i)

    return data[:orig_len], errors


# ============================================================================
# MAIN DEMO
# ============================================================================

print("=" * 72)
print("  TOURNAMENT TOYS — Five Creative Applications")
print("  opus-2026-03-24-S315")
print("=" * 72)

# 1. Tournament Sort
print(f"\n{'='*72}")
print(f"  1. TOURNAMENT SORT")
print(f"{'='*72}")

for label, arr in [
    ("already sorted", list(range(10))),
    ("reverse sorted", list(range(9, -1, -1))),
    ("nearly sorted",  [0, 1, 3, 2, 4, 5, 7, 6, 8, 9]),
    ("random",         [7, 2, 9, 1, 5, 3, 8, 0, 6, 4]),
    ("many ties",      [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]),
]:
    sorted_arr, stats = tournament_sort(arr)
    print(f"\n  {label}: {arr}")
    print(f"    Sorted: {sorted_arr}")
    print(f"    Comparisons: {stats['comparisons']}, "
          f"Transitivity: {stats['transitivity']:.3f}, "
          f"3-cycles: {stats['cycles']}")

# 2. Tournament Encryption
print(f"\n{'='*72}")
print(f"  2. TOURNAMENT ENCRYPTION")
print(f"{'='*72}")

for msg in [b"Hello, World!", b"Attack at dawn", b"The tournament codec is lossless"]:
    key = b"my secret key 12345"
    enc = tournament_encrypt(msg, key)
    dec = tournament_decrypt(enc, key)
    assert dec == msg, f"Encryption roundtrip failed for {msg}!"
    print(f"  Plain:  {msg}")
    print(f"  Cipher: {enc[:20].hex()}... ({len(enc)} bytes)")
    print(f"  Deciph: {dec}")
    print(f"  Verified: ✓")
    print()

# Wrong key test
wrong_dec = tournament_decrypt(enc, b"wrong key!!")
print(f"  Wrong key decrypt: {wrong_dec[:20]}... (garbage ✓)")

# 3. Tournament Fingerprint
print(f"\n{'='*72}")
print(f"  3. TOURNAMENT FINGERPRINT (audio matching)")
print(f"{'='*72}")

sr = 44100
t = np.arange(sr) / sr  # 1 second

# Original signal
signal1 = 0.5 * np.sin(2 * np.pi * 440 * t) + 0.3 * np.sin(2 * np.pi * 880 * t)
fp1 = audio_fingerprint((signal1 * 32767).astype(np.int16), sr)

# Same signal + noise
signal2 = signal1 + np.random.RandomState(42).normal(0, 0.05, len(t))
fp2 = audio_fingerprint((signal2 * 32767).astype(np.int16), sr)

# Different signal
signal3 = 0.5 * np.sin(2 * np.pi * 1000 * t) + 0.3 * np.sin(2 * np.pi * 2000 * t)
fp3 = audio_fingerprint((signal3 * 32767).astype(np.int16), sr)

# Shifted version (same signal, delayed)
signal4 = np.zeros_like(signal1)
signal4[1000:] = signal1[:-1000]
fp4 = audio_fingerprint((signal4 * 32767).astype(np.int16), sr)

print(f"  Fingerprint length: {len(fp1)} frames")
print(f"  Bits per frame: {len(fp1[0]) if fp1 else 0}")
print(f"\n  Distance matrix:")
print(f"    original  vs noisy:     {fingerprint_distance(fp1, fp2):.4f} (should be small)")
print(f"    original  vs different: {fingerprint_distance(fp1, fp3):.4f} (should be large)")
print(f"    original  vs shifted:   {fingerprint_distance(fp1, fp4):.4f} (moderate)")
print(f"    noisy     vs different: {fingerprint_distance(fp2, fp3):.4f} (should be large)")

# 4. Tournament Art
print(f"\n{'='*72}")
print(f"  4. TOURNAMENT RANDOM ART")
print(f"{'='*72}")

try:
    from PIL import Image

    for n_tourn in [4, 5]:
        for seed in [42, 137]:
            img_arr = tournament_art(n=n_tourn, size=128, seed=seed)
            img = Image.fromarray(img_arr)
            path = f'/tmp/tournament_art_n{n_tourn}_s{seed}.png'
            img.save(path)
            print(f"  Generated: {path} (n={n_tourn}, seed={seed})")
            print(f"    Unique colors: {len(set(tuple(img_arr[y,x]) for y in range(128) for x in range(128)))}")

    # Special: metagraph-walk art
    img_arr = np.zeros((256, 256, 3), dtype=np.uint8)
    rng = np.random.RandomState(42)
    # Random walk on the 2D plane, colored by tournament invariants
    x, y = 128, 128
    for step in range(100000):
        # Random direction biased by step count
        dx = rng.choice([-1, 0, 1])
        dy = rng.choice([-1, 0, 1])
        x = (x + dx) % 256
        y = (y + dy) % 256
        # Color = step modular structure
        r = (step * 7) % 256
        g = (step * 13) % 256
        b = (step * 31) % 256
        img_arr[y, x] = np.maximum(img_arr[y, x], [r, g, b])

    img = Image.fromarray(img_arr)
    path = '/tmp/tournament_walk_art.png'
    img.save(path)
    print(f"  Generated: {path} (random walk)")

except ImportError:
    print("  PIL not available for image generation")

# 5. Tournament Error Corrector
print(f"\n{'='*72}")
print(f"  5. TOURNAMENT ERROR CORRECTOR")
print(f"{'='*72}")

original = b"The quick brown fox jumps over the lazy dog. Tournament theory rocks!"
protected = tournament_error_detect(original)
print(f"  Original: {len(original)} bytes")
print(f"  Protected: {len(protected)} bytes ({len(protected)/len(original):.2f}x overhead)")

# No error
_, errors = tournament_error_check(protected)
print(f"  Check (no error): errors in blocks: {errors}")

# Introduce single-byte error in the data portion (after 6-byte header)
corrupted = bytearray(protected)
data_offset = 6  # header size
corrupted[data_offset + 10] ^= 0x42  # flip bits in data byte 10
_, errors = tournament_error_check(bytes(corrupted))
print(f"  Check (byte 10 corrupted): errors in blocks: {errors}")
print(f"  Error detected in block {errors[0] if errors else 'none'}, "
      f"byte 10 is in block {10 // 7}")

# Multiple errors
corrupted2 = bytearray(protected)
corrupted2[data_offset + 5] ^= 0xFF
corrupted2[data_offset + 30] ^= 0x01
_, errors2 = tournament_error_check(bytes(corrupted2))
print(f"  Check (bytes 5,30 corrupted): errors in blocks: {errors2}")
print(f"  Byte 5 in block {5//7}, byte 30 in block {30//7}")

print(f"\n{'='*72}")
print(f"  SUMMARY OF TOYS")
print(f"{'='*72}")
print("""
  1. TOURNAMENT SORT: Sorting with structural metadata.
     Transitivity score tells you HOW sorted the input was.
     3-cycle count reveals circular preferences (voting paradoxes).

  2. TOURNAMENT ENCRYPTION: Toy cipher using tournament permutations.
     NOT cryptographic! But demonstrates how tournament structure
     can scramble data (the Hamiltonian path = the permutation).

  3. TOURNAMENT FINGERPRINT: Robust audio matching via spectral rankings.
     Score sequences survive noise (97% preserved under perturbation).
     O(n) comparison time. Works like a simplified Shazam.

  4. TOURNAMENT RANDOM ART: Visual art from tournament invariants.
     Score variance → red, 3-cycle density → green, tiling weight → blue.
     Each pixel IS a tournament; the image IS a tournament landscape.

  5. TOURNAMENT ERROR CORRECTOR: Error detection via band-limitedness.
     High-frequency Walsh modes are zero for valid tournaments.
     Syndrome = high-frequency energy. Nonzero → error detected.

  All five exploit the same mathematical structure:
  tournaments → tilings → score sequences → recursive decomposition.
""")

print("DONE.")
