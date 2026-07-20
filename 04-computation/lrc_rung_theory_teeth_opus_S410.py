"""
opus-2026-07-19-S410 (HYP-8020): three in-frame verifications for the RUNG THEORY
manifesto (the manifesto ships only with same-day teeth).

T1  BOUNDARY DEGENERATION PROBE (the rung stack's moduli picture): the lift family
    V_L = {i + L*k_i} (k = (1,2,1,2,...,1,2), 12 speeds) should degenerate, as
    L grows, toward gridmax-land: its maximizer value should approach a limit and
    its ACTIVE STRUCTURE should stabilize to a pattern of the k-vector -- the
    concrete face of "families degenerate to subgroups" (G-K S*_0 boundary).
T2  THETA-FLOW TREE STRUCTURE: CF expansions of the AP13 maximizer path
    (1,14),(3,14),(23,107),(14,65),(11,51),(23,78),(31,112),(21,55),(13,34) --
    is the flow Stern-Brocot-coherent (shared deepening prefixes) or wild?
T3  DUTY-LEDGER SATURATION: for each corpus family with ghost structure, the
    budget slack 1 - [(m+1)M + M/K]; extremality = slack 0 (duty saturation).
"""
from math import gcd
from fractions import Fraction

def scan0(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q: bg, bq, wit = m, q, (a, q)
    return Fraction(bg, bq), wit

def cf(a, q):
    num, den, out = q, a, []
    while den: out.append(num//den); num, den = den, num % den
    return out  # CF of a/q as [a0; ...] with a0 = 0 implied by a<q? here num/den = q/a
                # NOTE: this returns CF of q/a; the CF of a/q is [0] + this.

print("T1: boundary degeneration of the lift family V_L = {i + L*k_i}, k = (1,2)^6:")
k = [1,2]*6
for L in (30, 60, 120, 240):
    V = [i+1 + L*k[i] for i in range(12)]
    M, (a, q) = scan0(V, 3*L + 60)
    act = [i+1 for i, v in enumerate(V)
           if Fraction(min((v*a) % q, q-(v*a) % q), q) == M]
    print(f"  L={L}: M = {M} = {float(M):.5f} at t*={a}/{q} (q/L = {q/L:.3f}); "
          f"active slots {act}")

print("\nT2: CF of the AP13 theta-path maximizers (CF of a/q = [0; ...]):")
path = [(1,14),(3,14),(23,107),(14,65),(11,51),(23,78),(31,112),(21,55),(13,34)]
for a, q in path:
    print(f"  {a}/{q}: [0; {cf(a % q, q)}]")

print("\nT3: duty-ledger slack 1 - [(m+1)M + M/K] (extremality = 0):")
fams = [("DW", list(range(1,13))+[182], 12, 14),
        ("ladder3", list(range(1,12))+[13,36], 11, 3),
        ("ladder4", list(range(1,12))+[13,48], 11, 4),
        ("GW", list(range(1,12))+[13,24], 11, 2),
        ("F4(31)", list(range(1,30))+[31,120], 29, 4)]
for name, V, m, K in fams:
    M, wit = scan0(V, 300)
    slack = 1 - ((m+1)*M + M*Fraction(1, K))
    print(f"  {name}: M={M}; budget slack = {slack} = {float(slack):.5f}"
          f"{'  SATURATED' if slack == 0 else ''}")
