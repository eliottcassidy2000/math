#!/usr/bin/env python3
"""
How small can a GMC counterexample be?   (mac-mini-2026-07-20-S133)
==================================================================
The published object is advertised as "cubic, 6 terms" in 4 real Gaussians.  The
mechanism run showed it is NOT term-minimal: dropping one of the six terms leaves a
counterexample in 2 of the 6 cases.  This finds the smaller ones and pushes further.

Criterion used throughout (the CORRECT one -- a Mathieu-Zhao subspace only requires
vanishing for m >> 0, so a counterexample must fail for infinitely many m):
    E[P^m] = 0 for all m = 1..M,  AND  E[Q P^m] != 0 at BOTH of the top two m.

Letters: Z, Zbar, W, Wbar  =  two independent standard complex Gaussians = 4 real.
    E[Z^a Zbar^b W^c Wbar^d] = delta_ab delta_cd a! c!
"""
from math import factorial
from itertools import product, combinations

def pmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2], k1[3]+k2[3])
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def expect(p):
    return sum(c * factorial(a) * factorial(cc)
               for (a, b, cc, d), c in p.items() if a == b and cc == d)

ONE = {(0,0,0,0): 1}
NAMES = ["Z", "Zb", "W", "Wb"]

def show(p):
    if not p: return "0"
    out = []
    for k, c in sorted(p.items()):
        s = "".join(NAMES[i] + ("" if k[i] == 1 else f"^{k[i]}")
                    for i in range(4) if k[i])
        out.append(("+" if c > 0 else "-") + (f"{abs(c)}" if abs(c) != 1 or not s else "") + s)
    return "".join(out).lstrip("+")

def test(P, Qs, M=9):
    """returns (Q, tail) for the first Q making P a counterexample, else None."""
    if not P: return None
    Pm = ONE; pw = []
    for m in range(1, M + 1):
        Pm = pmul(Pm, P)
        if expect(Pm) != 0: return None
        pw.append(Pm)
    for Q in Qs:
        vals = [expect(pmul(Q, x)) for x in pw]
        if vals[-1] != 0 and vals[-2] != 0: return (Q, vals)
    return None

def monos(dmax):
    out = []
    for d in range(1, dmax + 1):
        for k in product(range(d + 1), repeat=4):
            if sum(k) == d: out.append(k)
    return out

QS = [{k: 1} for k in monos(2)]

# ---------------------------------------------------------------- PART A
print("=" * 78)
print("PART A -- the 5-term counterexamples hiding inside the published 6-term one")
print("=" * 78)
Z = {(1,0,0,0):1}; ZB = {(0,1,0,0):1}; W = {(0,0,1,0):1}; WB = {(0,0,0,1):1}
def padd(*ps):
    out = {}
    for p in ps:
        for k, c in p.items(): out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}
P6 = pmul(padd(ONE, W), padd(pmul(ZB, padd(ONE, {k: -c for k, c in Z.items()})), WB))
print(f"  published P (6 terms) = {show(P6)}")
r = test(P6, QS)
print(f"    counterexample with Q = {show(r[0])},  E[QP^m] = {r[1][:6]}...")
print()
for k in list(P6):
    Pd = {kk: c for kk, c in P6.items() if kk != k}
    rr = test(Pd, QS)
    if rr:
        print(f"  DROP {show({k: P6[k]})}  ->  P = {show(Pd)}   ({len(Pd)} terms)")
        print(f"      Q = {show(rr[0])},  E[QP^m] = {rr[1][:6]}...")

# ---------------------------------------------------------------- PART B
print()
print("=" * 78)
print("PART B -- exhaustive: is there a QUADRATIC counterexample?  (deg <= 2)")
print("=" * 78)
ms2 = monos(2)
print(f"  {len(ms2)} monomials of degree <= 2; searching all supports of size <= 5, "
      f"coefficients in {{-1,+1}}")
best2 = []
for ksz in range(1, 6):
    hits = 0
    for supp in combinations(range(len(ms2)), ksz):
        for signs in product((-1, 1), repeat=ksz):
            P = {ms2[i]: s for i, s in zip(supp, signs)}
            rr = test(P, QS, M=8)
            if rr:
                hits += 1
                if len(best2) < 3: best2.append((P, rr))
    print(f"    support size {ksz}: {hits} counterexamples")
    if hits: break
if best2:
    for P, rr in best2:
        print(f"      P = {show(P)}   Q = {show(rr[0])}   E[QP^m] = {rr[1][:6]}...")
else:
    print("    => NO quadratic counterexample with <= 5 terms and unit coefficients.")
    print("       CUBIC appears to be necessary (search is bounded, not a proof).")

# ---------------------------------------------------------------- PART C
print()
print("=" * 78)
print("PART C -- exhaustive: fewest CUBIC terms?  (deg <= 3, supports of size <= 4)")
print("=" * 78)
ms3 = monos(3)
print(f"  {len(ms3)} monomials of degree <= 3")
found = {}
for ksz in range(1, 5):
    hits = []
    for supp in combinations(range(len(ms3)), ksz):
        for signs in product((-1, 1), repeat=ksz):
            P = {ms3[i]: s for i, s in zip(supp, signs)}
            rr = test(P, QS, M=8)
            if rr: hits.append((P, rr))
        if len(hits) > 40: break
    print(f"    support size {ksz}: {'>' if len(hits) > 40 else ''}{len(hits)} counterexamples")
    if hits:
        found[ksz] = hits
        for P, rr in hits[:4]:
            d = max(sum(k) for k in P)
            print(f"      P = {show(P):<28} (deg {d})  Q = {show(rr[0]):<6} "
                  f"E[QP^m] = {rr[1][:5]}...")
        break

print()
print("SUMMARY")
if found:
    k = min(found)
    print(f"  Smallest counterexample support found: {k} terms"
          f" (published claim was 6).")
else:
    print("  No counterexample with <= 4 terms and unit coefficients in this search box.")
print("  All searches are BOUNDED (unit coefficients, degree <= 3, support <= 5) --")
print("  absence of a hit is evidence, not proof.")
