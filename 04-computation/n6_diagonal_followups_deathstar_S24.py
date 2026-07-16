#!/usr/bin/env python3
"""death-star-2026-07-16-S24 (HYP-7051 part 2): census follow-up questions.
Q1 separation power: does adding d(W) to (|W|, H, c3, Aut, SC) refine the atlas partition?
Q2 affine dimension a(W) (rank of W - w0 over F2): the companion invariant; exceptionality
   (d(W) > counting-LB) vs affine structure.
Q3 d(W) of structured unions: all-SC tilings, c3-parity sets, H-level unions.
Q4 THE F7 BRIDGE TEST: for small LRC movies, the state set S in Z7^6 with movie adjacency;
   compute the F7 diagonal degree d7(S) (delta functions by degree-<=d polys over F7,
   evaluation rank) and compare avg movie degree vs candidate bounds c*d7 — evidence for a
   Z7 one-inclusion lemma (the 2607.14068 Lemma 2.3 analog) => movie rate bounds.
"""
from itertools import permutations, combinations, product as iproduct
from fractions import Fraction as Fr
from math import comb
from collections import defaultdict
import sys, time
exec(open('04-computation/n6_diagonal_census_deathstar_S24.py').read().split('if __name__')[0])

def affine_dim(W):
    w0 = W[0]
    basis = []
    for w in W[1:]:
        x = w ^ w0
        for b in basis:
            hb = b.bit_length() - 1
            if (x >> hb) & 1:
                x ^= b
        if x:
            basis.append(x); basis.sort(key=lambda z: -z.bit_length())
    return len(basis)

def f7_diag_degree(pts, m):
    """min d: delta functions on pts (in Z7^m) realizable by deg<=d polys over F7.
    Evaluation rank over F7 of monomials x^a (a_i <= 6, sum 'degrees' with x_i^k as deg k)."""
    import numpy as np
    P = 7
    npts = len(pts)
    # build monomial exponent lists by total degree, each exponent 0..6
    for d in range(0, 6 * m + 1):
        rows = []
        # exponents with sum = dd <= d, each 0..6 — generate incrementally
        def gen(prefix, rem, slots):
            if slots == 0:
                if rem == 0:
                    rows.append(tuple(prefix))
                return
            for e in range(0, min(6, rem) + 1):
                gen(prefix + [e], rem - e, slots - 1)
        allrows = []
        for dd in range(0, d + 1):
            gen([], dd, m)
        allrows = rows
        # evaluation matrix mod 7
        A = np.zeros((len(allrows), npts), dtype=np.int64)
        for ri, expo in enumerate(allrows):
            for ci, p in enumerate(pts):
                v = 1
                for e, x in zip(expo, p):
                    if e:
                        v = (v * pow(int(x), e, P)) % P
                A[ri, ci] = v
        # rank mod 7 via gaussian elimination
        Amat = A % P
        r = 0
        rows_ = Amat.copy()
        ncols = npts
        for c in range(ncols):
            piv = None
            for rr in range(r, rows_.shape[0]):
                if rows_[rr, c] % P:
                    piv = rr; break
            if piv is None: continue
            rows_[[r, piv]] = rows_[[piv, r]]
            inv = pow(int(rows_[r, c]), P - 2, P)
            rows_[r] = (rows_[r] * inv) % P
            for rr in range(rows_.shape[0]):
                if rr != r and rows_[rr, c] % P:
                    rows_[rr] = (rows_[rr] - rows_[rr, c] * rows_[r]) % P
            r += 1
            if r == min(rows_.shape[0], ncols): break
        if r == npts:
            return d
        rows = []
    return 6 * m

def movie_states_graph(E):
    mov = [f for f in E if f > 0]
    evs = []
    for i, f in enumerate(mov):
        for w in range(7 * f):
            evs.append((Fr(w, 7 * f), i))
    evs.sort()
    pos = []
    i = 0
    while i < len(evs):
        p = evs[i][0]
        while i < len(evs) and evs[i][0] == p: i += 1
        pos.append(p)
    n = len(pos)
    seq = []
    for j in range(n):
        a = pos[j]; b = pos[j + 1] if j + 1 < n else Fr(1)
        mid = (a + b) / 2
        seq.append(tuple(int((f * mid % 1) * 7) for f in mov))
    states = sorted(set(seq))
    stid = {s: i for i, s in enumerate(states)}
    edges = set()
    deg = defaultdict(int)
    for j in range(n):
        u, v = stid[seq[j]], stid[seq[(j + 1) % n]]
        edges.add((min(u, v), max(u, v)))
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    avg = 2 * len(edges) / len(states)
    return states, avg, n

if __name__ == "__main__":
    t0 = time.time()
    classes = {}
    for bits in range(1 << M):
        T = tourn(bits)
        classes.setdefault(canon_only(T), []).append(bits)
    rows = []
    for key, W in classes.items():
        T = tourn(W[0])
        _, aut = canon_and_aut(T)
        sc = (canon_only(rev(T)) == key)
        rows.append(dict(W=W, sz=len(W), H=H_count(T), c3=c3_of(T), aut=aut, sc=sc,
                         dW=diag_degree(W), clb=counting_lb(len(W)), aff=affine_dim(W)))
    print("Q1: separation power of the invariant tuples")
    def parts(keyf):
        d = defaultdict(int)
        for r in rows: d[keyf(r)] += 1
        return len(d), max(d.values())
    for name, kf in [("(sz,H,c3,Aut,SC)", lambda r: (r['sz'], r['H'], r['c3'], r['aut'], r['sc'])),
                     ("+d(W)", lambda r: (r['sz'], r['H'], r['c3'], r['aut'], r['sc'], r['dW'])),
                     ("+d(W)+aff", lambda r: (r['sz'], r['H'], r['c3'], r['aut'], r['sc'], r['dW'], r['aff']))]:
        ncl, mx = parts(kf)
        print(f"  {name}: {ncl} cells (max cell {mx}) of 56")
    print("\nQ2: affine dimension vs exceptionality")
    for r in sorted(rows, key=lambda r: (r['sz'],)):
        if r['dW'] > r['clb']:
            print(f"  EXC: sz={r['sz']} H={r['H']} SC={r['sc']} d={r['dW']} clb={r['clb']} aff={r['aff']} (m=10)")
    affs = [(r['aff'], r['sz']) for r in rows]
    full = sum(1 for a, s in affs if a == min(10, s - 1) or a == 10)
    print(f"  classes with full-possible affine dim: {full}/56")
    print("\nQ3: structured unions")
    scW = [w for r in rows if r['sc'] for w in r['W']]
    print(f"  all-SC tilings: |W|={len(scW)} d(W)={diag_degree(scW)} clb={counting_lb(len(scW))}")
    for par in [0, 1]:
        W = [w for r in rows if r['c3'] % 2 == par for w in r['W']]
        print(f"  c3-parity {par}: |W|={len(W)} d(W)={diag_degree(W)} clb={counting_lb(len(W))}")
    print(f"  [{time.time()-t0:.0f}s]")
    sys.stdout.flush()
    print("\nQ4: THE F7 BRIDGE TEST — movie state sets, avg degree vs F7 diagonal degree")
    for E, lab in [([0, 3, 4, 5, 6, 7], "[3..7] consec"),
                   ([0, 3, 5, 8, 11, 14], "small generic"),
                   ([0, 2, 4, 6, 8, 10], "dilate 2x[1..5]")]:
        states, avg, nev = movie_states_graph(E)
        if len(states) > 260:
            states_s = states[:260]
            note = f"(F7-deg on first 260 of {len(states)})"
        else:
            states_s = states; note = ""
        d7 = f7_diag_degree(states_s, len(E) - 1)
        print(f"  {lab}: V={len(states)} events={nev} avg-deg={avg:.2f} d7(S)={d7} {note} "
              f"ratio avg/d7={avg/d7 if d7 else float('nan'):.2f}")
        sys.stdout.flush()
    print(f"[total {time.time()-t0:.1f}s]")
