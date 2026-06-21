#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_S_formula_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

Pin the closed form of S(P7,Q7) = sum_j |e_j|  (e_j = 7 r_j - pq, residue-only).
And PROVE the residue-only law structurally.

Residue-only law (to be proven):  e_j depends only on (p mod 7, q mod 7).
Mechanism: e_j = 7 r_j - pq.  r_j counts the pq points  x=7frac(pk/q)+(2t+1)/(2q) in bin j.
Claim: 7 r_j - pq = G_j(p mod7, q mod7) for an explicit G.

We tabulate S as a 7x7 matrix over (P7,Q7), give its symmetry, and a closed form.
Conjecture from data: with a=P7 (=p mod7) and the slope s = a * (q^{-1} mod 7),
S = sum_j |e_j| equals 4 * (number of j with e_j>0 weighting)... let's just FIND it.

Observed values: S in {0,12,20,24,32,36,44}.  These are 4*{0,3,5,6,8,9,11}.
So S/4 in {0,3,5,6,8,9,11}.  Note 7 - {0,3,5,6,8,9,11}?? Let's relate to residues.
"""
from math import gcd

P = 7

def evec(p, q):
    r = [0] * P
    for k in range(q):
        base = 14 * ((p * k) % q)
        for t in range(p):
            j = ((base + 2 * t + 1) // (2 * q)) % P
            r[j] += 1
    return [7 * x - p * q for x in r]

def S_table():
    """S(a,b) for a=p%7,b=q%7 using a representative (smallest valid p,q in window)."""
    tab = {}
    # find a representative coprime (p,q) in window for each residue (a,b)
    reps = {}
    for q in range(1, 400):
        for p in range(q + 1, (43 * q) // 20 + 1):
            if 20 * p > 43 * q or gcd(p, q) != 1:
                continue
            key = (p % 7, q % 7)
            if key not in reps:
                reps[key] = (p, q)
    for key, (p, q) in reps.items():
        e = evec(p, q)
        tab[key] = sum(abs(x) for x in e)
    return tab, reps

def main():
    print("THREAD C: closed form of S(p mod 7, q mod 7)")
    print("=" * 72)
    tab, reps = S_table()
    print("S matrix (rows a=p%7=0..6, cols b=q%7=0..6); '.' = no window rep:")
    print("      b=0   1   2   3   4   5   6")
    for a in range(7):
        row = []
        for b in range(7):
            v = tab.get((a, b))
            row.append(f"{v:3d}" if v is not None else "  .")
        print(f"  a={a}: " + " ".join(row))

    # symmetry checks
    print("\nSymmetry S(a,b) =? S(b,a):",
          all(tab.get((a, b)) == tab.get((b, a)) for a in range(7) for b in range(7)
              if (a, b) in tab and (b, a) in tab))
    print("Symmetry S(a,b) =? S(-a,-b) (mod7):",
          all(tab.get((a, b)) == tab.get(((-a) % 7, (-b) % 7)) for (a, b) in tab
              if ((-a) % 7, (-b) % 7) in tab))
    print("S(a,b)=? S(7-a? ...) check S(a,b)=S(a,7-b):",
          all(tab.get((a, b)) == tab.get((a, (7 - b) % 7)) for (a, b) in tab
              if (a, (7 - b) % 7) in tab))

    # closed-form hunt: S in 4*{0,3,5,6,8,9,11}. Relate to slope s=a*binv mod7.
    print("\nClosed-form via slope s = (p mod7)*(q^{-1} mod 7) mod 7:")
    print("  (a,b) -> s, S, S/4")
    for b in range(1, 7):
        binv = pow(b, -1, 7)
        for a in range(1, 7):
            if (a, b) not in tab:
                continue
            s = (a * binv) % 7
            S = tab[(a, b)]
            print(f"   a={a} b={b}: s={s} S={S} S/4={S//4 if S%4==0 else S/4}")

    # The e-vector is a SHIFT of row-0 of the discrepancy. Let's see e directly as function of s.
    # group by s and see if S depends only on s
    print("\nDoes S depend only on slope s=a b^{-1} mod 7? group:")
    from collections import defaultdict
    byslope = defaultdict(set)
    for b in range(1, 7):
        binv = pow(b, -1, 7)
        for a in range(1, 7):
            if (a, b) in tab:
                byslope[(a * binv) % 7].add(tab[(a, b)])
    for s in sorted(byslope):
        print(f"   s={s}: S values {sorted(byslope[s])}")

if __name__ == "__main__":
    main()
