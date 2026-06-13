#!/usr/bin/env python3
"""
EXHAUSTIVE H-spectrum at n=8 (all 2^28 labeled tournaments).

Computes the complete Hamiltonian-path-count census H(T) for every tournament
on 8 vertices, via a numpy-vectorized batched Held-Karp DP parallelized over all
cores (no C compiler is available on this node). Uses path-reversal symmetry
H(T) = H(T^op) (T^op = reverse all arcs = bits XOR (2^m-1)) to process only the
lower half [0, 2^27); each complement pair shares an H value and exactly one
member lies in the lower half, so census counts are exactly doubled.

Performance notes (monad-compute): int64 numpy matmul uses a slow non-BLAS
integer path; int32 + an explicit elementwise contraction with a cache-friendly
batch (B=768) parallelizes ~24x better (~85k tournaments/s on 24 cores, ~24 min
for the 2^27 lower half). Chunks are processed IN ORDER so any partial result is
an exhaustive prefix [0, K); a JSON checkpoint records K and the running census
for resume.

FILLS A DOCUMENTED GAP: 05-knowledge/results/h_spectrum_n8_exhaustive.out
previously held only a header (the original C run never completed).

Settles OPEN-Q-055 sub-question "at what n does each forbidden value unlock?":
gives the exact achievable H-spectrum at n=8, exhaustively confirming whether
H=7, 21, 35, 39, 49, 63, ... are achievable at n=8 (upgrading HYP-1026 / HYP-1104
from n=8 *sampling* to n=8 *exhaustive*).

Correctness: the vectorized DP was cross-checked against a scalar reference
Held-Karp on 4096 contiguous + 500 random tournaments (0 mismatches); reversal
symmetry and oddness verified on 200 random n=8 tournaments; small-n spectra
(n=3..6) reproduce the known values [1,3], [1,3,5], [1,3,5,9,11,13,15], and the
n=6 spectrum.

Instance: monad-compute-2026-06-04-S?
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import time
import sys
import json
from multiprocessing import Pool

N = 8
M = N * (N - 1) // 2          # 28
EDGES = [(i, j) for i in range(N) for j in range(i + 1, N)]
FULL = (1 << N) - 1
HALF = 1 << (M - 1)           # 2^27 ; process [0, HALF)
B = 768                        # tournaments per chunk (cache-friendly)
ARANGE_M = np.arange(M)
CKPT = "05-knowledge/results/h_spectrum_n8_exhaustive_monad.ckpt.json"


def chunk_hist(args):
    lo, hi = args
    b = hi - lo
    bits = np.arange(lo, hi, dtype=np.uint32)
    bitk = ((bits[:, None] >> ARANGE_M[None, :]) & 1).astype(np.int32)
    adj = np.zeros((b, N, N), dtype=np.int32)
    for k in range(M):
        i, j = EDGES[k]
        bk = bitk[:, k]
        adj[:, i, j] = 1 - bk     # arc i->j when bit 0
        adj[:, j, i] = bk         # arc j->i when bit 1
    dp = np.zeros((b, 1 << N, N), dtype=np.int32)
    for v in range(N):
        dp[:, 1 << v, v] = 1
    for S in range(1, 1 << N):
        sub = dp[:, S, :]
        if not sub.any():
            continue
        cont = np.zeros((b, N), dtype=np.int32)
        for v in range(N):
            sv = sub[:, v]
            if not sv.any():
                continue
            cont += sv[:, None] * adj[:, v, :]
        for u in range(N):
            if S & (1 << u):
                continue
            dp[:, S | (1 << u), u] += cont[:, u]
    H = dp[:, FULL, :].sum(axis=1)
    vals, cnts = np.unique(H, return_counts=True)
    return vals, cnts


def main():
    tasks = [(lo, min(lo + B, HALF)) for lo in range(0, HALF, B)]
    print(f"H-SPECTRUM: n={N}, m={M}, EXHAUSTIVE ({1 << M} tournaments)")
    print(f"Processing lower half [0,{HALF}) in {len(tasks)} chunks of B={B}; "
          f"counts doubled by reversal symmetry. Cores: {os.cpu_count()}")
    sys.stdout.flush()

    hist = {}
    prefix = 0            # contiguous lower-half tournaments fully processed
    done = 0
    t0 = time.time()
    with Pool(processes=os.cpu_count()) as pool:
        # ordered imap -> results arrive in chunk order -> contiguous prefix
        for vals, cnts in pool.imap(chunk_hist, tasks, chunksize=1):
            for v, c in zip(vals.tolist(), cnts.tolist()):
                hist[v] = hist.get(v, 0) + c
            done += 1
            prefix = min(done * B, HALF)
            if done % 1024 == 0 or done == len(tasks):
                el = time.time() - t0
                rate = prefix / el
                eta = (HALF - prefix) / rate if rate else 0
                print(f"  {done}/{len(tasks)} chunks  prefix={prefix}/{HALF}  "
                      f"{rate:.0f}/s  ETA {eta:.0f}s  distinct={len(hist)}")
                sys.stdout.flush()
                with open(CKPT, "w") as f:
                    json.dump({"prefix_lowerhalf": prefix, "lowerhalf_total": HALF,
                               "hist_lowerhalf": {str(k): v for k, v in hist.items()}},
                              f)

    full_hist = {int(h): 2 * c for h, c in hist.items()}
    total = sum(full_hist.values())
    assert total == (1 << M), f"census total {total} != {1<<M}"

    achieved = sorted(full_hist)
    max_h = achieved[-1]
    odd_gaps = [h for h in range(1, max_h + 1, 2) if h not in full_hist]

    print("\n" + "=" * 64)
    print("FULL H-CENSUS (n=8, EXHAUSTIVE, all 2^28 = 268435456 tournaments)")
    print("=" * 64)
    for h in achieved:
        print(f"  H={h:6d}: {full_hist[h]:12d}")
    print("-" * 64)
    print(f"Distinct H values: {len(achieved)}")
    print(f"Range: [1, {max_h}]")
    print(f"All H odd: {all(h % 2 == 1 for h in achieved)}")
    print(f"Census total: {total} (= 2^28 = {1<<M})  OK")
    print(f"\nODD GAPS in [1, {max_h}]: {odd_gaps}")
    print("\nKEY FORBIDDEN-VALUE CHECKS at n=8:")
    for h in (7, 21, 35, 39, 49, 63):
        st = "ACHIEVED" if h in full_hist else "GAP"
        cc = f" (count={full_hist[h]})" if h in full_hist else ""
        print(f"  H={h}: {st}{cc}")
    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
