"""
benchmark_modp.py — Benchmark mod-p vs SVD chain complex computation

Verifies correctness and measures speedup of full_chain_complex_modp
over the SVD-based full_chain_complex_svd.

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament,
    full_chain_complex_modp, full_chain_complex_svd
)


def main():
    print("=" * 70)
    print("BENCHMARK: mod-p vs SVD chain complex computation")
    print("=" * 70)

    # Test 1: n=6 exhaustive correctness check
    print("\n--- Test 1: n=6 exhaustive correctness ---")
    n = 6
    total = 2 ** (n*(n-1)//2)
    mismatches = 0
    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        res_modp = full_chain_complex_modp(A, n, max_p=5)
        res_svd = full_chain_complex_svd(A, n, max_p=5)
        for p in range(6):
            if res_modp['bettis'][p] != res_svd['bettis'][p]:
                mismatches += 1
                if mismatches <= 3:
                    print(f"  MISMATCH at bits={bits}, p={p}: modp={res_modp['bettis'][p]}, svd={res_svd['bettis'][p]}")
                    print(f"    modp omega_dims={res_modp['omega_dims']}")
                    print(f"    svd  omega_dims={res_svd['omega_dims']}")
                    print(f"    modp ranks={res_modp['ranks']}")
                    print(f"    svd  ranks={res_svd['ranks']}")
    t1 = time.time()
    print(f"  n=6: {total} tournaments, {mismatches} mismatches, {t1-t0:.1f}s (both methods)")

    # Test 2: Speed comparison at n=7
    print("\n--- Test 2: n=7 speed comparison ---")
    n = 7
    rng = np.random.RandomState(42)
    count = 200

    # mod-p timing
    t0 = time.time()
    modp_bettis = []
    for _ in range(count):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=6)
        modp_bettis.append(tuple(res['bettis'].get(p, 0) for p in range(7)))
    t_modp = time.time() - t0

    # SVD timing (same tournaments)
    rng2 = np.random.RandomState(42)
    t0 = time.time()
    svd_bettis = []
    for _ in range(count):
        A = random_tournament(n, rng2)
        res = full_chain_complex_svd(A, n, max_p=6)
        svd_bettis.append(tuple(res['bettis'].get(p, 0) for p in range(7)))
    t_svd = time.time() - t0

    match = sum(1 for a, b in zip(modp_bettis, svd_bettis) if a == b)
    print(f"  n=7 ({count} tours): mod-p {t_modp:.2f}s, SVD {t_svd:.2f}s, speedup {t_svd/t_modp:.1f}x")
    print(f"  Betti agreement: {match}/{count}")
    if match < count:
        for i, (a, b) in enumerate(zip(modp_bettis, svd_bettis)):
            if a != b:
                print(f"    Mismatch {i}: modp={a}, svd={b}")
                if i > 5:
                    break

    # Test 3: Speed comparison at n=8
    print("\n--- Test 3: n=8 speed comparison ---")
    n = 8
    rng = np.random.RandomState(42)
    count = 20

    t0 = time.time()
    for _ in range(count):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=7)
    t_modp = time.time() - t0

    rng2 = np.random.RandomState(42)
    t0 = time.time()
    for _ in range(count):
        A = random_tournament(n, rng2)
        res = full_chain_complex_svd(A, n, max_p=7)
    t_svd = time.time() - t0

    print(f"  n=8 ({count} tours): mod-p {t_modp:.2f}s, SVD {t_svd:.2f}s, speedup {t_svd/t_modp:.1f}x")

    # Test 4: Betti numbers that the speedup gets right
    print("\n--- Test 4: Interesting Betti numbers ---")
    rng = np.random.RandomState(42)
    interesting = []
    for trial in range(2000):
        A = random_tournament(7, rng)
        res = full_chain_complex_modp(A, 7, max_p=6)
        b = res['bettis']
        if b.get(3, 0) > 0 or b.get(4, 0) > 0:
            interesting.append((trial, tuple(b.get(p, 0) for p in range(7))))

    print(f"  Found {len(interesting)} interesting n=7 tournaments (b3>0 or b4>0)")
    for trial, bettis in interesting[:10]:
        print(f"    trial={trial}: bettis={bettis}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
