#!/usr/bin/env python3
"""procgen_tcpc_20261001_orchestrator_check.py -- the orchestrator's independent audit of the tcpc lane (2026-10-01).

Written from the definitions in the lane note; the lane's code was not read. Engines: the orchestrator's own
exact arc-HP engine (procgen_petersen_20261001_orchestrator_check.c, general digraphs, N <= 18 here) and its own
mod-2 bitset engine (procgen_selfie_20261001_orchestrator_check.c, mode par2, tournaments, N <= 26).

Checks:
  A. D_r (rotational sheets S = {1..m}, cross set Y by r mod 4) is a tournament, x -> x+1 on both sheets is an
     automorphism, and psi o tau (sheet swap composed with a shift) is an anti-automorphism that is one 2r-cycle
     (so D_r is anti-circulant, cf. THM-4529).
  B. D_r is all-odd for r = 3, 5, 7, 9 (exact counts) and r = 11, 13 (mod 2); H(D_7) = 24540117, H(D_9) = 116670839805.
  C. D_3 = QR_7 - v, D_5 = QR_11 - v, D_7 = QR_127[mu_14] (isomorphism); D_9 is not QR_19 - v.
  D. T1: in random two-sheet clocks with reversed sheets and a symmetric cross set (X = -X), every vertical arc lies on an
     odd number of Hamiltonian paths (r = 3, 5, 7).
  E. Control: the plain interval cross set X = {0..m-1} is all-odd for r = 3, 7 and not for r = 5, 9.
  F. Syracuse clock step law: for odd A with 3 not dividing A, log_2(S(A)) = 2 (A mod 3) - v (mod 6) in (Z/9)^x
     (2 is a primitive root mod 9), i.e. S(A) mod 9 depends only on (A mod 3, v mod 6).
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile
import time

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
OKS = []
ARCS = PAR = None


def ok(cond, msg):
    OKS.append(bool(cond))
    print(('[OK] ' if cond else '[FAIL] ') + msg, flush=True)


def two_sheet(r, S, X):
    """vertices a_j = j, b_j = r + j"""
    N = 2 * r
    out = [set() for _ in range(N)]
    for j in range(r):
        for k in range(r):
            if j != k:
                if (k - j) % r in S:
                    out[j].add(k)              # a_j -> a_k
                if (j - k) % r in S:
                    out[r + j].add(r + k)      # b_j -> b_k
            if (k - j) % r in X:
                out[j].add(r + k)              # a_j -> b_k
            else:
                out[r + k].add(j)              # b_k -> a_j
    return [sorted(s) for s in out]


def D(r):
    m = (r - 1) // 2
    S = set(range(1, m + 1))
    if m % 2 == 1:
        Y = {0} | {d % r for t in range(1, (m - 1) // 2 + 1) for d in (t, -t)}
    else:
        Y = {d % r for t in range(1, m // 2 + 1) for d in (t, -t)}
    return two_sheet(r, S, Y)


def is_tournament(out):
    N = len(out)
    return all((j in out[i]) != (i in out[j]) for i in range(N) for j in range(N) if i != j)


def exact_arcs(out):
    N = len(out)
    masks = ['%x' % sum(1 << j for j in out[i]) for i in range(N)]
    r = subprocess.run([ARCS, str(N)] + masks, capture_output=True, text=True, check=True).stdout.split('\n')
    H = int(r[0].split()[1])
    c = {}
    for line in r[1:]:
        if line.strip():
            u, v, x = map(int, line.split())
            c[(u, v)] = x
    return H, c


def par2_all_odd(out):
    N = len(out)
    s = ''.join('1' if j in out[i] else '0' for i in range(N) for j in range(i + 1, N))
    r = subprocess.run([PAR, 'par2', s], capture_output=True, text=True, check=True).stdout
    # "par2 N=..: H mod 2 = h, #odd arcs = k of M"
    h = int(r.split('H mod 2 = ')[1].split(',')[0])
    k = int(r.split('#odd arcs = ')[1].split(' ')[0])
    M = int(r.split(' of ')[1].split()[0])
    return h, k, M


def qr_restricted(p, verts):
    Q = {(x * x) % p for x in range(1, p)}
    idx = {v: i for i, v in enumerate(verts)}
    return [sorted(idx[w] for w in verts if w != v and (w - v) % p in Q) for v in verts]


def G(out):
    return nx.DiGraph([(i, j) for i in range(len(out)) for j in out[i]])


def check_A_B_C():
    for r in (3, 5, 7, 9, 11, 13):
        out = D(r)
        N = 2 * r
        tour = is_tournament(out)
        tau = [(x + 1) % r if x < r else r + (x - r + 1) % r for x in range(N)]
        aut = all(((tau[y] in out[tau[x]]) == (y in out[x])) for x in range(N) for y in range(N) if x != y)
        # psi: a_j <-> b_{j}; look for a shift t with sigma = psi o tau^t an anti-automorphism that is a 2r-cycle
        found = None
        for t in range(r):
            sig = [r + (x + t) % r if x < r else (x - r + t) % r for x in range(N)]
            anti = all(((sig[x] in out[sig[y]]) == (y in out[x])) for x in range(N) for y in range(N) if x != y)
            if anti:
                cyc, y = 1, sig[0]
                while y != 0:
                    y = sig[y]
                    cyc += 1
                if cyc == N:
                    found = t
                    break
            # also psi with reflection a_j <-> b_{-j}
            sig = [r + (-x + t) % r if x < r else (-(x - r) + t) % r for x in range(N)]
            anti = all(((sig[x] in out[sig[y]]) == (y in out[x])) for x in range(N) for y in range(N) if x != y)
            if anti:
                cyc, y = 1, sig[0]
                while y != 0:
                    y = sig[y]
                    cyc += 1
                if cyc == N:
                    found = ('refl', t)
                    break
        ok(tour and aut and found is not None,
           f'A: D_{r} (N = {N}) is a tournament with translation automorphism and a single-2r-cycle anti-automorphism ({found})')
        if r <= 9:
            H, c = exact_arcs(out)
            allodd = all(v % 2 == 1 for v in c.values())
            extra = {7: 24540117, 9: 116670839805}.get(r)
            ok(allodd and (extra is None or H == extra), f'B: D_{r}: H = {H}, all {len(c)} arc counts odd (exact)')
        else:
            h, k, M = par2_all_odd(out)
            ok(h == 1 and k == M, f'B: D_{r} (N = {N}): H odd and {k} of {M} arcs odd (mod-2 bitset engine)')
    # isomorphisms
    iso3 = nx.is_isomorphic(G(D(3)), G(qr_restricted(7, list(range(1, 7)))))
    iso5 = nx.is_isomorphic(G(D(5)), G(qr_restricted(11, list(range(1, 11)))))
    mu14 = sorted({pow(2, k, 127) for k in range(7)} | {(-pow(2, k, 127)) % 127 for k in range(7)})
    iso7 = nx.is_isomorphic(G(D(7)), G(qr_restricted(127, mu14)))
    not9 = not nx.is_isomorphic(G(D(9)), G(qr_restricted(19, list(range(1, 19)))))
    ok(iso3 and iso5 and iso7 and not9, 'C: D_3 = QR_7 - v, D_5 = QR_11 - v, D_7 = QR_127[mu_14]; D_9 is not QR_19 - v')


def check_D():
    rng = random.Random(4532)
    total = bad = 0
    for r in (3, 5, 7):
        m = (r - 1) // 2
        for trial in range(25):
            S = set()
            for d in range(1, m + 1):
                S.add(d if rng.random() < 0.5 else (r - d) % r)
            X = set()
            if rng.random() < 0.5:
                X.add(0)
            for d in range(1, m + 1):
                if rng.random() < 0.5:
                    X.add(d)
                    X.add((r - d) % r)
            out = two_sheet(r, S, X)
            assert is_tournament(out)
            H, c = exact_arcs(out)
            for j in range(r):
                arc = (j, r + j) if (r + j) in out[j] else (r + j, j)
                total += 1
                if c[arc] % 2 == 0:
                    bad += 1
    ok(bad == 0 and total == 375, f'D: T1 on 75 random symmetric two-sheet clocks (r = 3, 5, 7): {total} vertical arcs, {bad} even')


def check_E():
    res = {}
    for r in (3, 5, 7, 9):
        m = (r - 1) // 2
        out = two_sheet(r, set(range(1, m + 1)), set(range(0, m)))
        H, c = exact_arcs(out)
        res[r] = all(v % 2 == 1 for v in c.values())
    ok(res == {3: True, 5: False, 7: True, 9: False}, f'E: interval cross set all-odd: {res}')


def check_F():
    log2 = {pow(2, k, 9): k for k in range(6)}
    bad = 0
    n = 0
    for A in range(1, 400001, 2):
        if A % 3 == 0:
            continue
        x = 3 * A + 1
        v = 0
        while x % 2 == 0:
            x //= 2
            v += 1
        n += 1
        if log2[x % 9] != (2 * (A % 3) - v) % 6:
            bad += 1
    ok(bad == 0, f'F: log_2(S(A)) = 2(A mod 3) - v (mod 6) for all {n} odd A < 4e5 prime to 3')


def main():
    global ARCS, PAR
    t0 = time.time()
    with tempfile.TemporaryDirectory() as d:
        ARCS = os.path.join(d, 'arcs')
        PAR = os.path.join(d, 'sa')
        subprocess.run(['cc', '-O2', '-o', ARCS, os.path.join(HERE, 'procgen_petersen_20261001_orchestrator_check.c')], check=True)
        subprocess.run(['cc', '-O2', '-o', PAR, os.path.join(HERE, 'procgen_selfie_20261001_orchestrator_check.c')], check=True)
        print('==== A/B/C ====', flush=True); check_A_B_C()
        print('==== D ====', flush=True); check_D()
        print('==== E ====', flush=True); check_E()
        print('==== F ====', flush=True); check_F()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
