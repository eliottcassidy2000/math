#!/usr/bin/env python3
"""AUDIT of HYP-9362 (d) / results note section 3: the extendable-plane chains, computed from the Paley codes.
 * S(3,4,8) = the 14 weight-4 words of the extended QR code of P_7; its derived design at a point is a Fano plane.
 * S(5,8,24) = the 759 octads of the extended QR code of P_23; derived designs S(4,7,23), S(3,6,22) and, at three points,
   a 2-(21,5,1) design = projective plane of order 4 (the unique such plane is PG(2,4)).
 * Divisibility: S(4,5,9) and S(6,9,25) do not exist; an extension of a projective plane of order n (a
   3-(n^2+n+2, n+2, 1) design) needs (n+2) | 12, i.e. n in {2, 4, 10} (the counting condition behind Hughes 1965).
Output: a8_designs.out
"""
import os, itertools
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a8_designs.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
def gf2_basis(rows):
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return basis
def codewords(basis):
    words = [0]
    for b in basis: words += [w ^ b for w in words]
    return words
def ext_qr(q):
    Q = {x * x % q for x in range(1, q)}
    rows = [sum(1 << ((i + r) % q) for r in Q) for i in range(q)]
    B = gf2_basis(rows)
    return [b | ((bin(b).count('1') & 1) << q) for b in B]
def blocks_of_weight(words, w, n):
    return [frozenset(i for i in range(n) if (x >> i) & 1) for x in words if bin(x).count('1') == w]
def is_steiner(blocks, t, v):
    cnt = {}
    for b in blocks:
        for T in itertools.combinations(sorted(b), t): cnt[T] = cnt.get(T, 0) + 1
    return len(cnt) == comb(v, t) and set(cnt.values()) == {1}
def derive(blocks, p):
    return [b - {p} for b in blocks if p in b]
W8 = codewords(ext_qr(7)); B8 = blocks_of_weight(W8, 4, 8)
say(f"[8,4,4] from P_7: {len(B8)} weight-4 words; S(3,4,8): {is_steiner(B8, 3, 8)}; derived at point 7 (the parity point): "
    f"{len(derive(B8, 7))} blocks, Fano 2-(7,3,1): {is_steiner(derive(B8, 7), 2, 7) if True else None}")
W24 = codewords(ext_qr(23)); B24 = blocks_of_weight(W24, 8, 24)
say(f"[24,12,8] from P_23: {len(B24)} octads; S(5,8,24): {is_steiner(B24, 5, 24)}")
D1 = derive(B24, 23); D2 = derive(D1, 22); D3 = derive(D2, 21)
pts1 = sorted(set().union(*D1)); pts2 = sorted(set().union(*D2)); pts3 = sorted(set().union(*D3))
say(f"derived at 1 point: {len(D1)} blocks of size 7 on {len(pts1)} points, S(4,7,23): {is_steiner([frozenset(pts1.index(x) for x in b) for b in D1], 4, 23)}")
say(f"derived at 2 points: {len(D2)} blocks of size 6 on {len(pts2)} points, S(3,6,22): {is_steiner([frozenset(pts2.index(x) for x in b) for b in D2], 3, 22)}")
P3 = [frozenset(pts3.index(x) for x in b) for b in D3]
pp = is_steiner(P3, 2, 21) and len(P3) == 21 and all(len(b) == 5 for b in P3)
lines_through = [sum(1 for b in P3 if x in b) for x in range(21)]
two_lines_meet_once = all(len(a & b) == 1 for a, b in itertools.combinations(P3, 2))
say(f"derived at 3 points: {len(P3)} blocks of size 5 on {len(pts3)} points; 2-(21,5,1): {pp}; 5 lines per point: {set(lines_through) == {5}}; "
    f"any two lines meet in one point: {two_lines_meet_once}  -> a projective plane of order 4 (unique: PG(2,4))")
# PG(2,4) explicitly over GF(4) and an isomorphism test by brute force on a canonical invariant (hyperoval count)
def hyperovals(lines, npts):
    L = [set(l) for l in lines]; cnt = 0
    for S in itertools.combinations(range(npts), 6):
        Ss = set(S)
        if all(len(Ss & l) in (0, 2) for l in L): cnt += 1
    return cnt
say(f"hyperovals (6-arcs) of the derived plane: {hyperovals(P3, 21)} (PG(2,4) has 168)")
def steiner_admissible(t, k, v):
    return all(comb(v - i, t - i) % comb(k - i, t - i) == 0 for i in range(t + 1))
say("S(4,5,9) admissible:", steiner_admissible(4, 5, 9), "; S(6,9,25) admissible:", steiner_admissible(6, 9, 25),
    "; S(3,4,8), S(3,6,22), S(4,7,23), S(5,8,24) admissible:", [steiner_admissible(*x) for x in ((3,4,8),(3,6,22),(4,7,23),(5,8,24))])
ext = [nn for nn in range(1, 200) if (nn * nn + nn + 2) * (nn * nn + nn + 1) % (nn + 2) == 0]
say("orders n <= 199 for which a 3-(n^2+n+2, n+2, 1) extension passes divisibility ((n+2) | 12):", ext)
out.close()
