"""Audit S50, second pass - b4: Steps 1 and 2 of the proof of THM-4745 (5), checked with own code.
  (S1a) configuration F: A_loc = NT_s (ending at v) u NT_r (starting at v), B_loc = C_3, arcs A_left' -> B' -> A_right';
        for s, r in {3,4,5} and EVERY back-arc path L_X v L_Y of A_left (L_Y starting at x_0) and EVERY path R_X v R_Y
        of A_right with R_X ending at the in-neighbour of v: L_X b1 R_X v L_Y b2 R_Y is a Hamiltonian path whose
        restriction to B is b1 v b2 (not a path).  F' (both blocks C_3, A_right' -> B' -> A_left'): b1 x0 x1 v z1 z2 b2.
  (S1b) the case analysis: v a junction of A_loc = NT_s u NT_r (s, r in {3,4,5}); B_loc any chain with parts in
        {2,3,4} and at most 6 vertices, v at any position of B_loc.  If T_loc is exact, both sub-gluings
        A_left u B_loc and A_right u B_loc are exact, hence (induction) one-way: so T_loc is one of the 4 combinations
        of one-way sub-gluings.  Check: in every combination where both sub-gluings ARE exact, T_loc is NOT exact.
  (S2)  Step 2: (a) both interior (s, r in {4,5}): the forced arcs of (F1),(B1) and (F2),(B2) contradict each other;
        (b) interior in A = NT_s (s = 4, 5), B = NT_r (r = 3..5) at an end (C_3: any vertex): with the forced arcs the
        pair (x_1..x_t x_N x_0 x_(t+1)..x_(N-1), Q_1 v b) has >= 2 merges (all Hamiltonian paths enumerated)."""
import itertools, time
from collections import Counter
import numpy as np
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def nt(s):
    return L.chain((s - 1,))


def is_path(adj, seq):
    return all((adj[seq[t]] >> seq[t + 1]) & 1 for t in range(len(seq) - 1))


def all_hp(adj, verts):
    verts = list(verts); S = 0
    for x in verts:
        S |= 1 << x
    out = []

    def dfs(p, used):
        if used == S:
            out.append(tuple(p)); return
        nxt = adj[p[-1]] & S & ~used
        while nxt:
            low = nxt & -nxt
            w = low.bit_length() - 1
            p.append(w); dfs(p, used | low); p.pop()
            nxt ^= low
    for s in verts:
        dfs([s], 1 << s)
    return out


# ------------------------------------------------------------------ (S1a)
okF = True; nw = 0
for s in (3, 4, 5):
    for r in (3, 4, 5):
        # vertices: A_left = 0..s-1 with v = s-1; A_right = v, s..s+r-2 ; B = v, b1, b2
        v = s - 1
        n = s + r - 1 + 2
        b1, b2 = n - 2, n - 1
        adj = [0] * n
        AL = nt(s)
        for x in range(s):
            for y in range(s):
                if (AL[x] >> y) & 1:
                    adj[x] |= 1 << y
        mapR = [v] + list(range(s, s + r - 1))
        AR = nt(r)
        for x in range(r):
            for y in range(r):
                if (AR[x] >> y) & 1:
                    adj[mapR[x]] |= 1 << mapR[y]
        for x in range(s - 1):                 # A_left' -> A_right' (chain arcs, forward)
            for z in range(s, s + r - 1):
                adj[x] |= 1 << z
        adj[v] |= 1 << b1; adj[b1] |= 1 << b2; adj[b2] |= 1 << v
        for x in range(s - 1):
            adj[x] |= (1 << b1) | (1 << b2)
        for z in range(s, s + r - 1):
            adj[b1] |= 1 << z; adj[b2] |= 1 << z
        Lpaths = [p for p in all_hp(adj, range(s)) if p.index(v) + 1 < len(p) and p[p.index(v) + 1] == 0]
        zin = mapR[r - 1]
        Rpaths = [p for p in all_hp(adj, mapR) if p.index(v) > 0 and p[p.index(v) - 1] == zin]
        for Lp in Lpaths:
            iv = Lp.index(v); LX, LY = list(Lp[:iv]), list(Lp[iv + 1:])
            for Rp in Rpaths:
                jv = Rp.index(v); RX, RY = list(Rp[:jv]), list(Rp[jv + 1:])
                W = LX + [b1] + RX + [v] + LY + [b2] + RY
                okF &= sorted(W) == list(range(n)) and is_path(adj, W) and not is_path(adj, [b1, v, b2])
                nw += 1
        okF &= len(Lpaths) == 2 ** (s - 2) and len(Rpaths) == 2 ** (r - 2)
# F'
adj = [0] * 7          # x0=0, x1=1, v=2, z1=3, z2=4, b1=5, b2=6
for u, w in [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2), (2, 5), (5, 6), (6, 2)]:
    adj[u] |= 1 << w
for x in (0, 1):
    for z in (3, 4):
        adj[x] |= 1 << z
for z in (3, 4):
    for b in (5, 6):
        adj[z] |= 1 << b
for b in (5, 6):
    for x in (0, 1):
        adj[b] |= 1 << x
Wp = [5, 0, 1, 2, 3, 4, 6]
okFp = is_path(adj, Wp) and not is_path(adj, [5, 2, 6]) and L.H(adj) > 27
check(okF, f'(S1a) configuration F, s, r in {{3,4,5}}: all {nw} witness choices are Hamiltonian paths restricting to b1 v b2')
check(okFp, "(S1a) configuration F': the witness b1 x0 x1 v z1 z2 b2 is valid and T is not exact")


# ------------------------------------------------------------------ (S1b)
def build(blocksA_left, blocksA_right, pb, jpos):
    """A_loc = IC(left parts + right parts) with v the junction; B_loc = IC(pb) with v at position jpos.
    Returns adjacency on A' u {v} u B', index lists, and helpers."""
    pa = blocksA_left + blocksA_right
    i = sum(blocksA_left)
    A = L.chain(pa); B = L.chain(pb)
    a, b = len(A), len(B)
    n = a + b - 1
    mA = {}; k = 0
    for x in range(a):
        if x == i:
            mA[x] = a - 1
        else:
            mA[x] = k; k += 1
    mB = {}; k = a
    for y in range(b):
        if y == jpos:
            mB[y] = a - 1
        else:
            mB[y] = k; k += 1
    base = [0] * n
    for x in range(a):
        for x2 in range(a):
            if (A[x] >> x2) & 1:
                base[mA[x]] |= 1 << mA[x2]
    for y in range(b):
        for y2 in range(b):
            if (B[y] >> y2) & 1:
                base[mB[y]] |= 1 << mB[y2]
    left = [mA[x] for x in range(0, i)]
    right = [mA[x] for x in range(i + 1, a)]
    Bp = [mB[y] for y in range(b) if y != jpos]
    v = a - 1
    return base, left, right, Bp, v, n, L.H(L.chain(blocksA_left)), L.H(L.chain(blocksA_right)), L.H(B)


def sub_H(adj, verts):
    verts = sorted(verts)
    idx = {x: r for r, x in enumerate(verts)}
    sub = [0] * len(verts)
    for x in verts:
        for y in verts:
            if (adj[x] >> y) & 1:
                sub[idx[x]] |= 1 << idx[y]
    return L.H(sub)


bad = []; tested = 0; ex_combo = 0
Bchains = [c for N in range(2, 6) for c in L.compositions(N, 2, 4)]
for sl in (2, 3, 4):
    for sr in (2, 3, 4):
        for pb in Bchains:
            for jpos in range(sum(pb) + 1):
                base, left, right, Bp, v, n, hL, hR, hB = build((sl,), (sr,), pb, jpos)
                for oL in ('fwd', 'bwd'):
                    for oR in ('fwd', 'bwd'):
                        adj = list(base)
                        for x, o in [(x, oL) for x in left] + [(x, oR) for x in right]:
                            for y in Bp:
                                if o == 'fwd':
                                    adj[x] |= 1 << y
                                else:
                                    adj[y] |= 1 << x
                        e1 = sub_H(adj, left + [v] + Bp) == hL * hB
                        e2 = sub_H(adj, right + [v] + Bp) == hR * hB
                        tested += 1
                        if e1 and e2:
                            ex_combo += 1
                            if L.H(adj) == hL * hR * hB:
                                bad.append(((sl,), (sr,), pb, jpos, oL, oR))
check(not bad, f'(S1b) junction of A_loc: in all {ex_combo} of {tested} one-way combinations where both sub-gluings are exact, '
      f'T_loc is not exact (A blocks of 3-5 vertices, B_loc any chain of <= 6 vertices, every position)' + (f' BAD {bad[:4]}' if bad else ''))


# ------------------------------------------------------------------ (S2)
ok_a = True
for s in (4, 5):
    for r in (4, 5):
        for t in range(1, s - 1):
            for u in range(1, r - 1):
                inA = set(range(t))                 # in-neighbours of x_t in NT_s (interior t)
                inB = set(range(u))
                # forced: every a in inA beats all of B'; every b in inB beats all of A'
                a0, b0 = min(inA), min(inB)
                ok_a &= True if (a0 in inA and b0 in inB) else False   # a0 -> b0 and b0 -> a0 both forced
check(ok_a, '(S2a) both interior: an in-neighbour of v in A and one in B are forced to beat each other (contradiction)')
ok_b = True; nb_ = 0
for s in (4, 5):
    N = s - 1
    for t in range(1, N):
        for r in (3, 4, 5):
            for end in ((0, r - 1) if r > 3 else (0, 1, 2)):
                A = nt(s); B = nt(r)
                n = s + r - 1
                mB = {}; k = s
                for y in range(r):
                    if y == end:
                        mB[y] = t
                    else:
                        mB[y] = k; k += 1
                adj = list(A) + [0] * (r - 1)
                for y in range(r):
                    for y2 in range(r):
                        if (B[y] >> y2) & 1:
                            adj[mB[y]] |= 1 << mB[y2]
                Ap = [x for x in range(s) if x != t]; Bp = list(range(s, n))
                for x in Ap:
                    for y in Bp:
                        if x < t:
                            adj[x] |= 1 << y
                        else:
                            adj[y] |= 1 << x
                P = list(range(1, t + 1)) + [N, 0] + list(range(t + 1, N))
                assert is_path(adj, P)
                Bset = set(Bp) | {t}
                Qs = [q for q in all_hp(adj, sorted(Bset)) if q.index(t) == len(q) - 2]
                Q = Qs[0]
                cnt = sum(1 for p in all_hp(adj, range(n)) if [x for x in p if x < s] == P and [x for x in p if x in Bset] == list(Q))
                ok_b &= cnt >= 2 and L.H(adj) > L.H(A) * L.H(B)
                nb_ += 1
check(ok_b, f'(S2b) interior in A, end of B: the stated pair has >= 2 merges and T_loc is not exact ({nb_} configurations)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
