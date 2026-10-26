#!/usr/bin/env python3
"""Independent numerical Ramanujan scan for G_q (weighted adjacency A on V-orbits) and for the twisted adjacency
A_chi of the branch-signed operator B_chi = P_U - P_A + P_D, primes QMIN <= q < QMAX.
Both operators are self-adjoint for the orbit weights |O|, so we diagonalise the symmetrised matrices
W^(1/2) A W^(-1/2) with eigvalsh (stable, real spectrum).  Records the extreme nontrivial eigenvalues (both edges)
and flags |mu| > 2 sqrt 3.  Borderline cases (|  |mu| - 2 sqrt 3 | < 1e-8) are listed for exact follow-up (PARI).
"""
import sys, time
import numpy as np
from sympy import primerange

R2 = 2 * np.sqrt(3.0)


def orbits_and_mats(q):
    INF = q
    inv = [0] * (q + 1)
    for t in range(1, q):
        inv[t] = pow(t, -1, q)
    def neg(t): return t if t == INF else (-t) % q
    def S(t): return 0 if t == INF else (INF if t == 0 else (-inv[t]) % q)
    def D(t): return INF if t == INF else (t + 2) % q
    V = [lambda t: t, neg, S, lambda t: S(neg(t))]
    chi = [1, -1, 1, -1]
    oid = [-1] * (q + 1); reps = []; size = []
    for t in range(q + 1):
        if oid[t] < 0:
            o = [v(t) for v in V]
            for p in o: oid[p] = len(reps)
            reps.append(min(o)); size.append(len(set(o)))
    n = len(reps)
    A = np.zeros((n, n))
    for a, s in enumerate(reps):
        for v in V:
            A[a, oid[D(v(s))]] += 1
    alive = []
    for a, s in enumerate(reps):
        if all(chi[k] == 1 for k, v in enumerate(V) if v(s) == s):
            alive.append(a)
    pos = {a: i for i, a in enumerate(alive)}
    Ac = np.zeros((len(alive), len(alive)))
    for i, a in enumerate(alive):
        s = reps[a]
        for k, v in enumerate(V):
            p = D(v(s)); b = oid[p]
            if b in pos:
                u = [j for j, w in enumerate(V) if w(reps[b]) == p][0]
                Ac[i, pos[b]] += chi[k] * chi[u]
    w = np.sqrt(np.array(size, dtype=float))
    As = A * w[:, None] / w[None, :]
    wa = w[alive]
    Acs = Ac * wa[:, None] / wa[None, :]
    return As, Acs


def main():
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    QMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    t0 = time.time()
    rows = []
    for q in primerange(QMIN, QMAX):
        As, Acs = orbits_and_mats(q)
        symerr = max(np.abs(As - As.T).max(), np.abs(Acs - Acs.T).max() if Acs.size else 0)
        ev = np.linalg.eigvalsh(As)
        assert abs(ev[-1] - 4) < 1e-9 and ev[-2] < 4 - 1e-6, (q, ev[-3:])
        lo, hi = ev[0], ev[-2]
        evc = np.linalg.eigvalsh(Acs) if Acs.size else np.array([0.0])
        rows.append((q, lo, hi, evc[0], evc[-1], symerr))
    nq = len(rows)
    exA = [(q, lo, hi) for q, lo, hi, *_ in rows if max(-lo, hi) > R2]
    exC = [(q, clo, chi_) for q, _, _, clo, chi_, _ in rows if max(-clo, chi_) > R2]
    border = [q for q, lo, hi, clo, chi_, _ in rows if min(abs(-lo - R2), abs(hi - R2), abs(-clo - R2), abs(chi_ - R2)) < 1e-8]
    wq, wr = max(((q, max(-lo, hi) / R2) for q, lo, hi, *_ in rows), key=lambda z: z[1])
    wqc, wrc = max(((q, max(-clo, chi_) / R2) for q, _, _, clo, chi_, _ in rows), key=lambda z: z[1])
    print(f"primes {QMIN} <= q < {QMAX}: {nq} primes; max symmetrisation error {max(r[5] for r in rows):.1e}")
    print(f"G_q (A): non-Ramanujan at {len(exA)} primes, fraction {len(exA)/nq:.4f}")
    print("   [q, lowest nontrivial mu, highest nontrivial mu, edge]:")
    for q, lo, hi in exA:
        edge = 'low (mu < -2 sqrt3)' if -lo > R2 else 'high (mu > 2 sqrt3)'
        if -lo > R2 and hi > R2: edge = 'both'
        print(f"   {q} {lo:.6f} {hi:.6f} {edge}")
    print(f"   worst ratio max|mu|/(2 sqrt 3) = {wr:.5f} at q = {wq}; largest |mu| = {max(max(-lo, hi) for _, lo, hi in exA) if exA else 0:.4f}")
    print(f"B_chi (A_chi): non-Ramanujan at {len(exC)} primes, fraction {len(exC)/nq:.4f}: "
          + ", ".join(f"{q} ({clo:.5f}, {chi_:.5f})" for q, clo, chi_ in exC))
    print(f"   worst ratio = {wrc:.5f} at q = {wqc}")
    print(f"borderline (within 1e-8 of 2 sqrt 3, need exact check): {border}")
    print(f"elapsed {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
