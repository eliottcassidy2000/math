#!/usr/bin/env python3
"""Independent mod-P singular-support engine for Vir_{p,q} (audit of HYP-9294 table / note 3.2).

For D <= Dmax: I_D = L_{-1} I_{D-1} + L_{-2} I_{D-2} (I_{2N} = <v>), dims checked against the
Rocha-Caridi character, grI by leading forms (PBW length), arc ideal J = (x_0^N)_d, and the
number of minimal differential-ideal generators of grI in each (weight, degree):
    dim grI_{D,r} - dim( sum_j x_j grI_{D-j-2,r-1} + d grI_{D-1,r} ).
Also dim ker(pi)_D = dim grI_D - dim J_D.
Usage: python vir_ss_audit.py p q Dmax [P]
"""
import sys
import time
import flint
from alin import rref_vectors, rank, nullspace, leading_forms
from virmod import VirVac, partitions_min2, poly_mul, poly_d

p, q, Dmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P = int(sys.argv[4]) if len(sys.argv) > 4 else 2147483647
one = flint.nmod(1, P)
c = one - one * 6 * (p - q) ** 2 / (p * q)
N = (p - 1) * (q - 1) // 2
W = 2 * N
t0 = time.time()
vir = VirVac(c, one)
bases = {D: partitions_min2(D) for D in range(Dmax + 1)}


def rocha_caridi(M):
    """char of L(c_{p,q},0) truncated: (1/(q)_inf) sum_k (q^{A_k} - q^{B_k})."""
    ser = [0] * (M + 1)
    for k in range(-50, 51):
        a = ((2 * p * q * k + q - p) ** 2 - (q - p) ** 2) // (4 * p * q)
        b = ((2 * p * q * k + q + p) ** 2 - (q - p) ** 2) // (4 * p * q)
        if 0 <= a <= M:
            ser[a] += 1
        if 0 <= b <= M:
            ser[b] -= 1
    for n in range(1, M + 1):  # divide by (q;q)_inf  -> multiply by 1/(1-q^n)
        for i in range(n, M + 1):
            ser[i] += ser[i - n]
    return ser


rc = rocha_caridi(Dmax)

# singular vector
eq = {}
for n in (1, 2):
    for b in bases[W]:
        for m, cc in vir.act(n, b).items():
            eq.setdefault((n, m), {})[b] = cc
ker = nullspace(list(eq.values()), bases[W], P)
assert len(ker) == 1, len(ker)
v = ker[0]
I = {W: [v]}
for D in range(W + 1, Dmax + 1):
    gens = [vir.actv(-1, w) for w in I[D - 1]]
    if D - 2 >= W:
        gens += [vir.actv(-2, w) for w in I[D - 2]]
    I[D] = [row for _, row in rref_vectors(gens, bases[D], P)]
ok_rc = all(len(bases[D]) - len(I.get(D, [])) == rc[D] for D in range(Dmax + 1))
print(f"Vir({p},{q}) N={N} mod {P}: dim V_D - dim I_D == Rocha-Caridi for D<={Dmax}: {ok_rc}", flush=True)

gr = {}
for D in range(W, Dmax + 1):
    lf = leading_forms(I[D], bases[D], degree=len, modulus=P)
    for g in lf:
        r = len(next(iter(g)))
        gr.setdefault((D, r), []).append(g)


def x_times(j, g):
    return poly_mul({(j + 2,): one}, g)


gens_found = []
kers = []
f = {(2,) * N: one}
dfs = [f]
for _ in range(Dmax):
    dfs.append(poly_d(dfs[-1]))
for D in range(W, Dmax + 1):
    degs = sorted({r for (DD, r) in gr if DD == D})
    for r in degs:
        target = gr[(D, r)]
        cols = [m for m in bases[D] if len(m) == r]
        sub = []
        for g in gr.get((D - 1, r), []):
            sub.append(poly_d(g))
        for j in range(0, D - W - 2 + 1):
            for g in gr.get((D - j - 2, r - 1), []):
                sub.append(x_times(j, g))
        rs = rank(sub, cols, P) if sub else 0
        ng = len(target) - rs
        if ng:
            gens_found.append((D, r, ng))
    # J and ker(pi)
    Jg = []
    for s in range(0, D - W + 1):
        for m in bases[D - W - s]:
            Jg.append(poly_mul({m: one}, dfs[s]))
    dJ = rank(Jg, bases[D], P)
    kers.append((D, len(I[D]) - dJ))
print("minimal generators (weight, degree, multiplicity):", gens_found)
print("dim ker(pi)_D:", [k for k in kers if k[1]])
print("first nonzero ker(pi) weight:", next((D for D, k in kers if k), None), " (2N+3 =", W + 3, ")")
print(f"time {time.time() - t0:.1f}s")
