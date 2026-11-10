#!/usr/bin/env python3
"""HYP-9657: brute-force check of the F1 family in psi-coordinates (f1_mass.py uses this formula).

J = A => B with the single arc between the two apexes reversed (apex of B beats apex of A), A = one-apex tournament
with word u, B = one-apex tournament with word v (u_i = 1 iff the apex beats chain vertex i).
(1) join law:   H(J) = g(u) g(v) - e s + (1 + n(u)) (1 + U(v)),  n(u) = value of u read with u_1 leading,
    U(v) = sum over zeros of v at position j of 2^(j-1), e = [u ends in 0 or is empty], s = [v starts with 1 or is empty];
(2) psi form:   H(J) = (2 psi(h1) + 1)(2 psi(h2) + 1) + 2^s (2 h1 + 1)(2 h2 + 1) - [s = 0], where 2 h1 is the value of
    the reduced core of u, 2 h2 the value of the reverse-complement of the reduced core of v, and s counts the
    trailing ones of u plus the leading zeros of v.
Universe: all pairs of words of lengths 0..4 (961 joins), counts by subset DP.
Run:  python3 f1_formula_check.py > f1_formula_check.out
"""
import itertools
def ham(n, beats):
    size = 1 << n
    cnt = [[0]*n for _ in range(size)]
    for v in range(n): cnt[1 << v][v] = 1
    for mask in range(1, size):
        row = cnt[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length()-1] += c
                    nxt ^= w
    return sum(cnt[size-1])
def one_apex_beats(u, offset):
    # chain vertices offset..offset+m-1 (increasing), apex offset+m ; u[i] = 1 iff apex -> chain_i
    m = len(u); arcs = []
    for i in range(m):
        for j in range(i+1, m): arcs.append((offset+i, offset+j))
        if u[i]: arcs.append((offset+m, offset+i))
        else: arcs.append((offset+i, offset+m))
    return arcs, offset + m
def Hword(u):
    arcs, apex = one_apex_beats(u, 0)
    n = len(u) + 1
    beats = [0]*n
    for a, b in arcs: beats[a] |= 1 << b
    return ham(n, beats)
def psi(h):
    return 0 if h == 0 else psi(h >> 1) + ((h >> 1) if h % 2 == 0 else ((h >> 1)//2 + 1))
bad = 0; tested = 0; badpsi = 0
for mu in range(0, 5):
    for mv in range(0, 5):
        for u in itertools.product([0,1], repeat=mu):
            for v in itertools.product([0,1], repeat=mv):
                arcsA, x = one_apex_beats(u, 0)
                arcsB, y = one_apex_beats(v, mu + 1)
                n = mu + mv + 2
                beats = [0]*n
                for a, b in arcsA + arcsB: beats[a] |= 1 << b
                for a in range(mu + 1):
                    for b in range(mu + 1, n):
                        if (a, b) == (x, y): beats[b] |= 1 << a          # back arc y -> x
                        else: beats[a] |= 1 << b
                H = ham(n, beats)
                gu, gv = Hword(u), Hword(v)
                nu = int("".join(map(str, u)), 2) if mu else 0            # MSB-first value of u
                Uv = sum(1 << j for j in range(mv) if v[j] == 0)          # zeros of v, position j (0-based) -> 2^j
                eu = 1 if (mu == 0 or u[-1] == 0) else 0
                sv = 1 if (mv == 0 or v[0] == 1) else 0
                pred = gu * gv - eu * sv + (1 + nu) * (1 + Uv)
                tested += 1
                if pred != H: bad += 1
                # psi-form: strip trailing ones of u and leading zeros of v (paddings s), reduced cores have even value 2h
                uu = list(u); s = 0
                while uu and uu[-1] == 1: uu.pop(); s += 1
                while uu and uu[0] == 0: uu.pop(0)
                vv = list(v)
                while vv and vv[0] == 0: vv.pop(0); s += 1
                while vv and vv[-1] == 1: vv.pop()
                h1 = (int("".join(map(str, uu)), 2) // 2) if uu else 0
                # reduced v: starts with 1, ends with 0; its dual coordinate: 1 + U(v) = 1 + (reverse-complement value) = 2 h2 + 1 with 2 h2 = value of rc(v)
                rc = [1 - b for b in reversed(vv)]
                h2 = (int("".join(map(str, rc)), 2) // 2) if rc else 0
                predpsi = (2*psi(h1) + 1) * (2*psi(h2) + 1) + (1 << s) * (2*h1 + 1) * (2*h2 + 1) - (1 if s == 0 else 0)
                if predpsi != H: badpsi += 1
print("tested joins:", tested, " join-law mismatches:", bad, " psi-form mismatches:", badpsi)
import sys
print("ALL CHECKS PASSED" if bad == 0 and badpsi == 0 else "SOME CHECK FAILED")
sys.exit(0 if bad == 0 and badpsi == 0 else 1)
