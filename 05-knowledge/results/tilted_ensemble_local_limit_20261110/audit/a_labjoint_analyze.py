#!/usr/bin/env python3
"""AUDIT: analysis of the auditor's own labelled joint distributions (a_labjoint.c output, files labjoint_N.txt in
/tmp/claude-audit-ensemble or the directory given as argv[1]).  Checks: Theorem B at n = 3..7 as polynomial identities,
second-moment polynomials in w, residues mod 8, exact mod-4 discrepancies, and a file-by-file comparison with the
packet's out/lab_N.txt."""
import sys, os
from fractions import Fraction as Fr
from math import comb, factorial
D = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-audit-ensemble"
here = os.path.dirname(os.path.abspath(__file__))
def load(path): return [tuple(map(int, l.split())) for l in open(path) if l.strip()]
def pmul(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B): C[i + j] += x * y
    return C
def ppow(A, e):
    R = [1]
    for _ in range(e): R = pmul(R, A)
    return R
def padd(A, B):
    n = max(len(A), len(B)); return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(n)]
def strip(A):
    A = list(A)
    while len(A) > 1 and A[-1] == 0: A.pop()
    return A
def in_w(A):
    A = strip(A); out = {}
    while A != [0]:
        deg = len(A) - 1; assert deg % 2 == 0, "not symmetric"
        k = deg // 2; c = A[-1] * (-1) ** k; out[k] = c
        A = strip(padd(A, [-c * x for x in ppow([0, 1, -1], k)]))
    return [out.get(k, 0) for k in range(max(out) + 1)]
claimed_g = {3: [1, 2], 4: [1, 8], 5: [1, 22, 16], 6: [1, 52, 136], 7: [1, 114, 720, 272]}
claimed_sq = {3: [1, 8], 4: [1, 40, 16], 5: [1, 152, 544, 56, 48], 6: [1, 552, 6288, 6944, 3072], 7: [1, 2040, 55760, 189848, 199464, 72640, 4128]}
claimed_disc = {3: Fr(1, 2), 4: Fr(1, 2), 5: Fr(13, 32), 6: Fr(1, 4), 7: Fr(853, 8192)}
claimed_mod8 = {5: [.352, .234, .352, .063], 6: [.303, .176, .322, .199], 7: [.295, .246, .257, .202]}
Wn = {3: 8, 4: 32, 5: 158, 6: 928, 7: 6350}
ok = True
for n in range(3, 8):
    rows = load(os.path.join(D, f"labjoint_{n}.txt")); m = n * (n - 1) // 2
    assert sum(c for _, _, c in rows) == 2 ** m
    mean = [0]; sq = [0]
    for b, h, c in rows:
        basis = pmul(ppow([0, 1], m - b), ppow([1, -1], b))
        mean = padd(mean, [c * h * x for x in basis]); sq = padd(sq, [c * h * h * x for x in basis])
    gw = in_w(mean); sw = in_w(sq)
    d4 = Fr(abs(sum(c if h % 4 == 1 else -c for _, h, c in rows)), 2 ** m)
    r8 = [sum(c for _, h, c in rows if h % 8 == r) / 2 ** m for r in (1, 3, 5, 7)]
    allodd = all(h % 2 == 1 for _, h, _ in rows)
    e2 = sum(Fr(c, 4 ** j) for j, c in enumerate(sw)) == Fr(factorial(n) * Wn[n], 4 ** (n - 1))
    same = None
    pk = os.path.join(here, "..", "out", f"lab_{n}.txt")
    if os.path.exists(pk): same = sorted(load(pk)) == sorted(rows)
    good = (gw == claimed_g[n] and sw == claimed_sq[n] and d4 == claimed_disc[n] and allodd and e2 and same is not False)
    if n in claimed_mod8: good &= all(abs(r8[i] - claimed_mod8[n][i]) < 6e-4 for i in range(4))
    ok &= good
    print(f"n={n}: E_p[H] in w = {gw} (README {claimed_g[n]}); E_p[H^2] in w = {sw}; H always odd: {allodd}")
    print(f"      mod-4 discrepancy |P(1)-P(3)| = {d4} = {float(d4):.6f} (README {claimed_disc[n]}); residues mod 8 (1,3,5,7) = {[round(x, 4) for x in r8]}")
    print(f"      E[H^2] at w=1/4 = n! W(n)/4^(n-1) with W(n)={Wn[n]}: {e2}; identical to packet out/lab_{n}.txt: {same}; row OK: {good}")
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
