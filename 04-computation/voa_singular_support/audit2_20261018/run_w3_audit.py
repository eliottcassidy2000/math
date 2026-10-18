#!/usr/bin/env python3
"""Audit of HYP-9297 / HYP-9298 (vertex-algebra side).  usage: python run_w3_audit.py p q Dmax [P]
W_3 minimal model W_3(p,q), c = 2(1 - 12(p-q)^2/(pq)).  Maximal submodule = radical of the contravariant form.
Prints dim L, compares with (a) chi^{(5,6)}_0 + chi^{(5,6)}_3 for (4,5), (b) the A_2 product for p = 3;
prints the Li-degree-0 generators of gr I as polynomials in x = sigma(L_{-2}), y = sigma(W_{-3}) (mod P and as
reconstructed rationals), the arc-space Hilbert series, ker pi and the minimal generators of gr I."""
import sys, time
from fractions import Fraction
from aud_algebras import W3, LW, WW
from vamod import GradedAnalysis

p, q, Dmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P = int(sys.argv[4]) if len(sys.argv) > 4 else 2147483647
c = 2 * (1 - Fraction(12 * (p - q) ** 2, p * q))
va = W3(c, P)
W2max = 2 * Dmax
t0 = time.time()
G = va.gram_all(W2max)
I = va.radical(G)
print(f'W3({p},{q}) c={c} P={P}: gram+radical to weight {Dmax} in {time.time()-t0:.1f}s', flush=True)
dimV = [len(va.basis(2 * w)) for w in range(Dmax + 1)]
dimI = [len(I.get((2 * w, 0), [])) for w in range(Dmax + 1)]
dimL = [a - b for a, b in zip(dimV, dimI)]
print('dim V:', dimV)
print('dim L:', dimL)


def ratrec(a, P):
    """rational reconstruction of a mod P"""
    a %= P
    r0, r1, s0, s1 = P, a, 0, 1
    import math
    B = int(math.isqrt(P // 2))
    while r1 > B:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def vir_char(pp, qq_, r, s, N):
    """Rocha-Caridi character of Vir(pp,qq_) module (r,s), normalised to start at q^0, coefficients to q^N"""
    h = Fraction((pp * r - qq_ * s) ** 2 - (pp - qq_) ** 2, 4 * pp * qq_)
    num = [0] * (N + 1)
    for k in range(-50, 51):
        e1 = Fraction((2 * pp * qq_ * k + pp * r - qq_ * s) ** 2 - (pp * r - qq_ * s) ** 2, 4 * pp * qq_)
        e2 = Fraction((2 * pp * qq_ * k + pp * r + qq_ * s) ** 2 - (pp * r - qq_ * s) ** 2, 4 * pp * qq_)
        for e, sg in ((e1, 1), (e2, -1)):
            assert e.denominator == 1
            e = int(e)
            if 0 <= e <= N:
                num[e] += sg
    # divide by (q)_inf
    pn = [0] * (N + 1); pn[0] = 1
    for n in range(1, N + 1):
        for i in range(n, N + 1):
            pn[i] += pn[i - n]
    out = [sum(num[i] * pn[n - i] for i in range(n + 1)) for n in range(N + 1)]
    return h, out


if (p, q) == (4, 5):
    # 3-state Potts: Vir(5,6) modules h=0 (r,s)=(1,1) and h=3 (r,s)=(1,5) [Kac table of (5,6): h_{r,s}]
    N = Dmax
    h0, ch0 = vir_char(5, 6, 1, 1, N)
    # find (r,s) with h=3
    found = None
    for r in range(1, 5):
        for s in range(1, 6):
            h, _ = vir_char(5, 6, r, s, 0)
            if h == 3:
                found = (r, s)
    h3, ch3 = vir_char(5, 6, found[0], found[1], N)
    pred = [ch0[n] + (ch3[n - 3] if n >= 3 else 0) for n in range(N + 1)]
    print('chi^(5,6)_0 + q^3 chi^(5,6)_3 (h=3 at (r,s)=%s):' % (found,), pred)
    print('   MATCH' if pred == dimL else '   MISMATCH at %s' % [n for n in range(N + 1) if pred[n] != dimL[n]])
if p == 3:
    a = [0] * (Dmax + 1)
    for n in range(1, Dmax + 1):
        r = n % q
        a[n] = (0 if r in (0, 1, q - 1) else 1) + (0 if r in (0, 1, q - 1, 2, q - 2) else 1)
    ser = [1] + [0] * Dmax
    for n in range(1, Dmax + 1):
        for _ in range(a[n]):
            for i in range(n, Dmax + 1):
                ser[i] += ser[i - n]
    print('A_2 product mod', q, ':', ser)
    print('   MATCH' if ser == dimL else '   MISMATCH at %s' % [n for n in range(Dmax + 1) if ser[n] != dimL[n]])

t0 = time.time()
GA = GradedAnalysis(va, I, W2max)
GA.run()
summ = GA.summary_by_weight()
print(f'graded analysis in {time.time()-t0:.1f}s')
# Li-0 generators
print('Li-0 part of gr I (C_2 ideal) generators, as polynomials in x=L_{-2}, y=W_{-3}:')
for (D, nV, nI, nJ, ng) in GA.table:
    if D[2] == 0 and ng:
        # minimal generators: print basis of grI_D modulo Dec_D -- just print grI basis when Dec empty
        for v in GA.grI.get(D, []):
            terms = []
            for m, cc in sorted(v.items()):
                nx = sum(1 for X in m if X[0] == LW); ny = sum(1 for X in m if X[0] == WW)
                terms.append((ratrec(cc, P), f'x^{nx} y^{ny}'))
            print(f'   weight {D[0]//2}: (dim grI={nI}, dim decomposable={nI-ng}) ', terms)
print(' weight  dimV  dimI  dimL  arcHS  ker   minimal generators of grI (Li, -, #)')
for W2 in range(0, W2max + 1, 2):
    if W2 not in summ:
        continue
    o = summ[W2]
    gens = o['gens']
    flag = '  <-- NON-ARC' if any(g[0] > 0 for g in gens) else ''
    print(f" {W2//2:6d} {o['V']:5d} {o['I']:5d} {o['V']-o['I']:5d} {o['V']-o['J']:6d} {o['ker']:4d}   {gens}{flag}")
