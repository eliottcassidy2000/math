#!/usr/bin/env python3
"""collatz_ck_fixed_points_20260930_klein_farey.py -- two follow-ups of the thirteenth note.
 (1) D42: locally-C_7 Cayley graphs of S_4 of degree 7 (Klein's {3,7} map on the genus-3 surface is one: PSL(2,7) contains S_4
     acting regularly on its 24 vertices); their induced 7-cycles, the rhombus graph, and whether C_7(G) iso G.
 (2) D40: the denominator law for mediants: with beta(u,v) = S_(uv) - S_(vu) = D_u S_v - D_v S_u, one has
     D_u S_(uv) = -2^(A_u) beta(u,v) mod D_(uv), so when gcd(D_u, D_v) = 1 the denominator of x_(uv) is
     D_(uv)/gcd(beta(u,v), D_(uv)); checked on all pairs of words with A_u, A_v <= 7.
"""
import itertools, time
from math import gcd
from collatz_ck_fixed_points_20260930 import small_groups, symmetric_subsets, cayley, is_locally_ck, induced_cycles, cycle_operator, is_isomorphic, rhombus_graph
from collatz_argument_styles_20260930 import carry, clock

T0 = time.time()
print("== (1) locally-C_7 Cayley graphs of S_4 of degree 7 ==")
els, mul, inv, ident = small_groups()["S_4"]
cnt = 0; loc = []
for S in symmetric_subsets(els, mul, inv, ident, 7):
    if len(S) != 7:
        continue
    cnt += 1
    G = cayley(els, mul, inv, S)
    if is_locally_ck(G, 7):
        loc.append((S, G))
print(" degree-7 connection sets: %d; locally C_7: %d (%.0fs)" % (cnt, len(loc), time.time() - T0))
seen_iso = []
for S, G in loc:
    if any(is_isomorphic(G, H) for H in seen_iso):
        continue
    seen_iso.append(G)
    sev = induced_cycles(G, 7, limit=5 * 24)
    R = rhombus_graph(G)
    line = " locally-C_7 graph #%d: induced 7-cycles %s" % (len(seen_iso), len(sev) if len(sev) <= 5 * 24 else ">120")
    if len(sev) == 24:
        C7, _ = cycle_operator(G, 7)
        line += "; C_7(G) iso G: %s; C_7(G) iso R(G): %s; R(G) degree %s, R iso G: %s" % (is_isomorphic(C7, G), is_isomorphic(C7, R), sorted(set(len(a) for a in R)), is_isomorphic(R, G))
    else:
        line += " (extra non-contractible 7-cycles; not a fixed point)"
    # triangles per vertex and girth-type data
    tri = sum(1 for u in G[0] for w in G[0] if u < w and w in G[u])
    line += "; triangles at a vertex %d" % tri
    print(line + " (%.0fs)" % (time.time() - T0), flush=True)

print("== (2) the denominator law for weighted mediants ==")
words = []
for A in range(1, 8):
    for cuts in range(1 << (A - 1)):
        w = []; run = 1
        for i in range(A - 1):
            if (cuts >> i) & 1:
                w.append(run); run = 1
            else:
                run += 1
        w.append(run); words.append(tuple(w))
ok_cong = True; ok_law = True; n_coprime = 0; n_pairs = 0; counter = None
for u in words:
    Su, Du, Au = carry(u), clock(u), sum(u)
    for v in words:
        Sv, Dv = carry(v), clock(v)
        Suv, Duv = carry(u + v), clock(u + v)
        beta = Su * Dv - Sv * Du
        n_pairs += 1
        ok_cong &= (Du * Suv + (2 ** Au) * beta) % Duv == 0
        if gcd(Du, Dv) == 1:
            n_coprime += 1
            den = abs(Duv) // gcd(Suv, abs(Duv))
            law = abs(Duv) // gcd(beta, abs(Duv))
            if den != law:
                ok_law = False; counter = counter or (u, v, den, law)
print(" congruence D_u S_(uv) = -2^(A_u) beta mod D_(uv) on all %d pairs: %s" % (n_pairs, ok_cong))
print(" denominator law den(x_(uv)) = |D_(uv)|/gcd(beta, D_(uv)) on the %d coprime-clock pairs: %s %s" % (n_coprime, ok_law, "" if ok_law else "counterexample %s" % (counter,)))
print(" so in the coprime-clock case the mediant is integral iff the clock of the whole divides the Farey determinant of the parts")
