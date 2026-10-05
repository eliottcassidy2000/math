#!/usr/bin/env python3
"""
Syracuse-law partition functions to level 20 (opus, 2026-10-05): deciding equality in the spectrum bound.

Same recursion as collatz_syracuse_multifractal_20261005.py,
    mu_n(c) = sum_{a >= 1, 2^a c = 1 mod 3} 2^-a mu_(n-1)((2^a c - 1)/3 mod 3^(n-1)),
with a parallel numba kernel.  Levels <= LFULL are stored in float64; the level LFULL+1 is STREAMED (never stored):
its partition functions Z(q) = sum_c mu(c)^q, the atoms at c = -1, -5, -17, 1 and the maximum are accumulated on
the fly.  Output: per-level slopes s_n(q) = -log(Z_n(q)/Z_(n-1)(q))/log 3 against the bound
tau(q) = min(q-1, log_3(2^q-1)), and the deficit d_n(q) = bound - s_n(q), whose decay decides equality:
summable deficits (e.g. n^-3/2) mean equality (polynomial prefactor), a positive limit means strict inequality.
Usage: python collatz_syracuse_multifractal_deep_20261005.py [LFULL=19] [A=44]
"""
import sys, math, time, json
import numpy as np
import numba as nb
from numba import prange

QS = np.array([0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0, 4.0, 6.0, 8.0])

@nb.njit(parallel=True, cache=True)
def level_full(prev, n, A, pw2):
    """compute mu_n (float64, length 3^n) from prev = mu_(n-1); pw2[a] = 2^a mod 3^n (int64)."""
    M = 3 ** n; Mp = 3 ** (n - 1)
    out = np.zeros(M, dtype=np.float64)
    for c in prange(M):
        r = c % 3
        if r == 0:
            continue
        # 2^a c = 1 mod 3: a even iff c = 1 mod 3; a odd iff c = 2 mod 3
        a0 = 2 if r == 1 else 1
        s = 0.0
        w = 0.5 ** a0
        for a in range(a0, A + 1, 2):
            y = ((pw2[a] * c - 1) // 3) % Mp
            s += w * prev[y]
            w *= 0.25
        out[c] = s
    return out

@nb.njit(parallel=True, cache=True)
def level_stream(prev, n, A, pw2, qs, nchunks, atoms_idx):
    """stream mu_n from prev = mu_(n-1): returns (Z(q) array, max, atoms values)."""
    M = 3 ** n; Mp = 3 ** (n - 1)
    nq = qs.shape[0]
    Zpart = np.zeros((nchunks, nq))
    mxpart = np.zeros(nchunks)
    chunk = (M + nchunks - 1) // nchunks
    for ch in prange(nchunks):
        lo = ch * chunk; hi = min(M, lo + chunk)
        zloc = np.zeros(nq); mx = 0.0
        for c in range(lo, hi):
            r = c % 3
            if r == 0:
                continue
            a0 = 2 if r == 1 else 1
            s = 0.0
            w = 0.5 ** a0
            for a in range(a0, A + 1, 2):
                y = ((pw2[a] * c - 1) // 3) % Mp
                s += w * prev[y]
                w *= 0.25
            if s > 0.0:
                ls = math.log(s)
                for i in range(nq):
                    zloc[i] += math.exp(qs[i] * ls)
                if s > mx:
                    mx = s
        for i in range(nq):
            Zpart[ch, i] = zloc[i]
        mxpart[ch] = mx
    Z = np.zeros(nq)
    for i in range(nq):
        Z[i] = Zpart[:, i].sum()
    atoms = np.zeros(atoms_idx.shape[0])
    for k in range(atoms_idx.shape[0]):
        c = atoms_idx[k] % M
        r = c % 3
        if r != 0:
            a0 = 2 if r == 1 else 1
            s = 0.0
            w = 0.5 ** a0
            for a in range(a0, A + 1, 2):
                y = ((pw2[a] * c - 1) // 3) % Mp
                s += w * prev[y]
                w *= 0.25
            atoms[k] = s
    return Z, mxpart.max(), atoms

def main():
    LFULL = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 44
    out = []; rec = {}
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    def bound(q):
        return min(q - 1, math.log(2 ** q - 1) / math.log(3)) if q > 1 else (q - 1)
    t0 = time.time()
    P(f"deep Syracuse partition functions: full levels <= {LFULL}, streamed level {LFULL + 1}; valuations <= {A}; q = {list(QS)}")
    mu = np.array([1.0]); Zprev = None
    for n in range(1, LFULL + 2):
        M = 3 ** n
        pw2 = np.array([pow(2, a, M) for a in range(A + 1)], dtype=np.int64)
        atoms_idx = np.array([-1, -5, -17, 1, 5, 7], dtype=np.int64)
        if n <= LFULL:
            new = level_full(mu, n, A, pw2)
            mu = new
            pos = mu[mu > 0]
            Z = np.array([float(np.sum(pos ** q)) for q in QS])
            mx = float(mu.max()); atoms = np.array([mu[i % M] for i in atoms_idx])
            mass = float(mu.sum())
        else:
            Z, mx, atoms = level_stream(mu, n, A, pw2, QS, 1024, atoms_idx)
            mass = float("nan")
        line = f"n={n:2d} [{time.time()-t0:6.0f}s] mass {mass:.10f} max {mx:.4e} (x2^n {mx * 2.0 ** n:.4f}) atoms(-1,-5,-17,1,5,7) alpha_n: " + " ".join(f"{-math.log(v)/(n*math.log(3)):.4f}" if v > 0 else "nan" for v in atoms)
        rec[n] = {"Z": [float(z) for z in Z], "max": mx, "atoms": [float(v) for v in atoms]}
        if Zprev is not None:
            slopes = [-math.log(Z[i] / Zprev[i]) / math.log(3) for i in range(len(QS))]
            rec[n]["slopes"] = slopes
            line += " | deficits bound-slope: " + " ".join(f"q{QS[i]}:{bound(QS[i]) - slopes[i]:+.4f}" for i in range(len(QS)))
        P(line)
        Zprev = Z
        with open(__file__.replace(".py", ".json"), "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1)
        with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(out) + "\n")
    P("")
    P("deficit decay test: d_n * n^1.5 (constant => deficits summable => equality with a polynomial prefactor; growing => strict inequality)")
    for q_i, q in enumerate(QS):
        if q <= 1: continue
        vals = []
        for n in range(10, LFULL + 2):
            if n in rec and "slopes" in rec[n]:
                d = bound(q) - rec[n]["slopes"][q_i]
                vals.append(f"{n}:{d * n ** 1.5:.3f}")
        P(f"  q={q}: " + " ".join(vals))
    P(f"  [{time.time()-t0:.0f}s total]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
