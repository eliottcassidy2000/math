#!/usr/bin/env python3
"""
ramsey_hybrid.py — two-phase witness search for R(B_{n-1},B_n) >= 4n-1:
Phase 1: structured SA over two-level model (U0,U1 symmetric, D free) — fast
         exact FFT energies, gets within a few dozen violations.
Phase 2: general-graph tabu from the phase-1 adjacency — free of structural
         constraints, fixes the residual.
Phase 2 uses batched candidate moves with exact incremental C = A^2.

Usage: python3 ramsey_hybrid.py n [p1_iters] [p2_iters] [seed] [rounds]
"""
import numpy as np
import sys, time
from ramsey_book_sa2 import energy as sa2_energy, brute_verify, qr_set
from ramsey_tabu_general import energy_full, flip_delta, apply_flip, pair_energy

def phase1(n, iters, rng):
    N = 4*n - 2; m = 2*n - 1
    U0 = np.zeros(m); U1 = np.zeros(m); D = np.zeros(m)
    if rng.random() < 0.5:
        D = qr_set(m).astype(np.float64)
    else:
        D[rng.permutation(m)[:m//2]] = 1
    for a in range(1, (m+1)//2):
        if rng.random() < 0.5: U0[a] = U0[m-a] = 1
        if rng.random() < 0.5: U1[a] = U1[m-a] = 1
    e = sa2_energy(U0, U1, D, n, N)
    best = e; bU0, bU1, bD = U0.copy(), U1.copy(), D.copy()
    T0 = 2.5
    for it in range(iters):
        T = T0 * max(0.01, 1 - it/iters)
        mv = rng.random()
        if mv < 0.25:
            a = rng.integers(1, m); b = (m - a) % m
            U0[a] = 1 - U0[a]; U0[b] = U0[a]; undo = ("U0", a, b)
        elif mv < 0.5:
            a = rng.integers(1, m); b = (m - a) % m
            U1[a] = 1 - U1[a]; U1[b] = U1[a]; undo = ("U1", a, b)
        else:
            a = int(rng.integers(0, m)); D[a] = 1 - D[a]; undo = ("D", a, a)
        e2 = sa2_energy(U0, U1, D, n, N)
        if e2 <= e or rng.random() < np.exp((e - e2)/max(T, 1e-9)):
            e = e2
            if e < best:
                best = e; bU0, bU1, bD = U0.copy(), U1.copy(), D.copy()
            if e == 0:
                break
        else:
            typ, a, b = undo
            if typ == "U0": U0[a] = 1 - U0[a]; U0[b] = U0[a]
            elif typ == "U1": U1[a] = 1 - U1[a]; U1[b] = U1[a]
            else: D[a] = 1 - D[a]
    # build adjacency from best
    U0, U1, D = bU0, bU1, bD
    A = np.zeros((N, N), dtype=np.int8)
    for i in range(m):
        for j in range(m):
            dd = (j - i) % m
            if i != j:
                if U0[dd]: A[i, j] = 1
                if U1[dd]: A[m+i, m+j] = 1
            if D[dd]:
                A[i, m+j] = 1; A[m+j, i] = 1
    return A, best

def phase2(A, n, iters, rng, tag=""):
    N = A.shape[0]
    e, C, deg = energy_full(A, n, N)
    C = C.astype(np.int64); deg = deg.astype(np.int64)
    best = e
    tabu = {}
    t0 = time.time()
    last_improve = 0
    for it in range(iters):
        bsz = 32
        us = rng.integers(0, N, size=bsz); vs = rng.integers(0, N, size=bsz)
        bestdelta = None; bestuv = None
        for u, v in zip(us, vs):
            u = int(u); v = int(v)
            if u == v: continue
            if tabu.get((min(u,v), max(u,v)), 0) > it: continue
            dlt = flip_delta(A, C, deg, u, v, n, N)
            if bestdelta is None or dlt < bestdelta:
                bestdelta = dlt; bestuv = (u, v)
        if bestuv is None: continue
        u, v = bestuv
        if bestdelta <= 0 or rng.random() < 0.03:
            apply_flip(A, C, deg, u, v)
            e += bestdelta
            tabu[(min(u,v), max(u,v))] = it + int(rng.integers(10, 60))
            if e < best:
                best = e; last_improve = it
                if best < 6:
                    print(f"  {tag} it={it} e={e} ({time.time()-t0:.0f}s)", flush=True)
            if e == 0:
                e2, _, _ = energy_full(A, n, N)
                print(f"  {tag} ZERO at it={it}, recheck={e2}", flush=True)
                if e2 == 0:
                    return A, 0
                e = e2
        if it - last_improve > 60000:
            break
    return A, best

def main():
    n = int(sys.argv[1])
    p1 = int(sys.argv[2]) if len(sys.argv) > 2 else 150000
    p2 = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    rounds = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    N = 4*n - 2
    rng = np.random.default_rng(seed)
    globalbest = 10**9
    for r in range(rounds):
        A, e1 = phase1(n, p1, rng)
        efull, _, _ = energy_full(A, n, N)
        print(f"[round {r}] phase1 struct-energy={e1} full-energy={efull}", flush=True)
        A, e2 = phase2(A, n, p2, rng, tag=f"[r{r}]")
        print(f"[round {r}] phase2 final e={e2}", flush=True)
        globalbest = min(globalbest, e2)
        if e2 == 0:
            fn = f"ramsey_hybrid_witness_n{n}_seed{seed}_r{r}"
            np.save(fn + ".npy", A)
            with open(fn + ".txt", "w") as f:
                f.write("".join(str(A[i, j]) for j in range(N) for i in range(N)) + "\n")
            print("WITNESS SAVED " + fn, flush=True)
            # independent brute verify
            ok, won, woff, _ = brute_verify_from_A(A, n, N)
            print(f"final brute verify: ok={ok} worst_on={won}(<= {n-2}) worst_off={woff}(<= {n-1})", flush=True)
            return
    print(f"no witness, best={globalbest}", flush=True)

def brute_verify_from_A(A, n, N):
    M = A.astype(np.int64)
    C = M @ M
    comp = 1 - M
    np.fill_diagonal(comp, 0)
    Cc = comp @ comp
    ok = True; won = -1; woff = -1
    for i in range(N):
        for j in range(i+1, N):
            if M[i, j]:
                won = max(won, C[i, j])
                if C[i, j] > n-2: ok = False
            else:
                woff = max(woff, Cc[i, j])
                if Cc[i, j] > n-1: ok = False
    return ok, won, woff, A

if __name__ == "__main__":
    main()
