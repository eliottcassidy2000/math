# Audit (HYP-9371 section 8), task 2b: the sign-network engine against the first-principles kernel on ALL 1076
# excess-4 cores (including the 7- and 8-vertex cores that decomp.c cannot enumerate, MAXE = 10) and on the
# excess-5 cores with V <= 5.
# For path lengths lambda_e, the engine's kernel is the projection of its 5-basis weight tensor W on the z^lambda
# coefficients of the five basis functions,  K_sn(lambda) = sum_b W[b] prod_e c_(lambda_e)(b_e),
#   c_l(R_+) = 1, c_l(R_-) = (-1)^l, c_l(Q_(s1 s2)) = [x^l] (-r_s1(x) r_s2(x) (1 - x^2)/(1 + x^2)),
# obtained by contracting the sign network with the projected edge tensors (label axis = path length 1..Lam), so the
# 5^E tensor is never formed.  a5_kernel.c computes K(lambda) from the THM-4757 definition.  Checked:
#   (i)  lambda = (1, ..., 1) on every core (the [z^E] coefficient: E = 5..12 at excess 4, E <= 10 at excess 5);
#   (ii) NRAND random lambda in {1,2,3}^E with sum(lambda - 1) <= 3 per core (junctions / m2 paths);
#   (iii) every length vector of the z^10 kernel table (a5_task5_K_L10.json, excess 1..5), all orders lambda_e <= 6.
# Usage: python3 a5_task2b_kernel_vs_engine.py KERNEL_BINARY [NRAND]
import os, sys, json, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, D); sys.path.insert(0, HERE)
import loop_engine5 as LE5
import a5_z10_lib as Z

binary = sys.argv[1]; NRAND = int(sys.argv[2]) if len(sys.argv) > 2 else 4
LAMMAX = 10
_LOADED = {}
def cores_of(e):
    if e not in _LOADED: _LOADED[e] = Z.load(Z.CORE_FILES[e])
    return _LOADED[e]

def series_coeffs(Lam):
    N = Lam + 1
    rp = [0] + [1] * (N - 1); rm = [0] + [(-1) ** k for k in range(1, N)]
    cm = [1] + [0] * (N - 1)
    for j in range(1, N):
        if 2 * j < N: cm[2 * j] = 2 * (-1) ** j
    def mul(a, b): return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(N)]
    Q = {(s1, s2): [-x for x in mul(mul(rp if s1 == 1 else rm, rp if s2 == 1 else rm), cm)] for s1 in (1, -1) for s2 in (1, -1)}
    C = np.zeros((Lam, 5), dtype=np.int64)            # C[l-1, b]
    for l in range(1, Lam + 1):
        C[l - 1] = [rp[l], rm[l], Q[(1, 1)][l], Q[(1, -1)][l], Q[(-1, -1)][l]]
    assert Q[(1, -1)] == Q[(-1, 1)]
    return C
CC = series_coeffs(LAMMAX)
ops5 = LE5.ops_basis5()
def ops_proj(Lam):
    C = CC[:Lam].T                                     # (5, Lam)
    return {'E': np.einsum('ijb,bl->ijl', ops5['E'], C), 'F+': ops5['F+'] @ C, 'F-': ops5['F-'] @ C}

def engine_K(V, edges, aut, lams):
    Lam = max(max(l) for l in lams)
    T = LE5.SignCore(V, edges, aut).contract('basis5', ops_proj(Lam))
    return [int(T[tuple(x - 1 for x in l)]) for l in lams]

rng = random.Random(20261022)
t0 = time.time()
tot = 0; bad = 0
for e in (4, 5):
    cores = cores_of(e); byV = {}
    for i, (V, edges, aut) in enumerate(cores):
        E = len(edges)
        lams = [tuple([1] * E)]
        if e == 4:
            for _ in range(NRAND):
                l = [1] * E
                for _ in range(rng.randint(1, 3)): l[rng.randrange(E)] += 1
                if max(l) <= 3: lams.append(tuple(l))
        mine = Z.kernels(binary, V, edges, lams)
        eng = engine_K(V, edges, aut, lams)
        nb = sum(a != b for a, b in zip(mine, eng)); bad += nb; tot += len(lams)
        s = byV.setdefault(V, [0, 0, 0, 0]); s[0] += 1; s[1] += len(lams); s[2] += nb; s[3] += sum(1 for x in mine if x)
        if nb: print(f"   excess {e} core {i} (V={V}): MISMATCH {[(l, a, b) for l, a, b in zip(lams, mine, eng) if a != b][:3]}", flush=True)
    for V, s in sorted(byV.items()):
        print(f"excess {e}, V = {V}: {s[0]} cores, {s[1]} length vectors ({s[3]} with K != 0), {s[2]} mismatches", flush=True)
    print(f"   ({time.time() - t0:.0f}s)", flush=True)

# (iii) the full z^10 table
tab = json.load(open(os.path.join(HERE, 'a5_task5_K_L10.json')))
nz = 0; nb3 = 0; nv = 0
for key, Ks in tab.items():
    L, e, i = map(int, key.split(':'))
    V, edges, aut = cores_of(e)[i]
    loops = {k for k, (a, b) in enumerate(edges) if a == b}
    lams = list(Z.compositions(L, len(edges), loops))
    if not lams: continue
    eng = engine_K(V, edges, aut, lams)
    nb3 += sum(a != b for a, b in zip(Ks, eng)); nv += len(lams); nz += sum(1 for x in Ks if x)
print(f"z^10 kernel table: {len(tab)} cores, {nv} length vectors ({nz} with K != 0), {nb3} mismatches against the engine")
print(f"TOTAL: {tot + nv} kernel values compared, {bad + nb3} mismatches ({time.time() - t0:.0f}s)")
