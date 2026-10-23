# Audit (HYP-9371 section 8), task 2c: variant of a5_task2b_kernel_vs_engine.py that runs ONLY part (ii') at the end
# (random length vectors with EVEN total on all 1076 excess-4 cores; 2b's random vectors mostly had odd totals, where
# every kernel vanishes identically).  Shared description of 2b: the sign-network engine against the first-principles kernel on ALL 1076
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

rng = random.Random(4757)
# (ii') every excess-4 core: NRAND random length vectors with EVEN total (odd totals vanish identically: every run of a
# contributing decomposition has even length), extras sum(lambda - 1) in 1..5, entries <= 4; plus lambda = 1 when E is even.
t0 = time.time(); tot = 0; bad = 0
cores = cores_of(4); byV = {}
for i, (V, edges, aut) in enumerate(cores):
    E = len(edges)
    lams = [tuple([1] * E)] if E % 2 == 0 else []
    tries = 0
    while len(lams) < NRAND + (E % 2 == 0) and tries < 1000:
        tries += 1
        l = [1] * E
        for _ in range(rng.randint(1, 5)): l[rng.randrange(E)] += 1
        if sum(l) % 2 == 0 and max(l) <= 4 and tuple(l) not in lams: lams.append(tuple(l))
    mine = Z.kernels(binary, V, edges, lams)
    eng = engine_K(V, edges, aut, lams)
    nb = sum(a != b for a, b in zip(mine, eng)); bad += nb; tot += len(lams)
    s_ = byV.setdefault(V, [0, 0, 0, 0, 0]); s_[0] += 1; s_[1] += len(lams); s_[2] += nb; s_[3] += sum(1 for x in mine if x)
    s_[4] = max(s_[4], max(abs(x) for x in mine))
    if nb: print(f"   core {i} (V={V}): MISMATCH {[(l, a, b) for l, a, b in zip(lams, mine, eng) if a != b][:3]}", flush=True)
for V, s_ in sorted(byV.items()):
    print(f"excess 4, V = {V}: {s_[0]} cores, {s_[1]} even-length vectors ({s_[3]} with K != 0, max |K| {s_[4]}), {s_[2]} mismatches", flush=True)
print(f"TOTAL: {tot} kernel values compared, {bad} mismatches ({time.time() - t0:.0f}s)")
