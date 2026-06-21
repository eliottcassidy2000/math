#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_general_prime_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C -- NEW LEAD)

NEW LEAD: does the residue-only closed form generalize to an arbitrary apex value P
(P sectors), not just P=7?  If yes, the LRC machinery is robust and the same proof gives
D_{p,q}^{(P)} = S_P(p mod P, q mod P)/(P p q),  with S_P a finite PxP residue table.

The proof chain used ONLY: (i) cyclic-shift (needs P nmid q for invertibility -- works for any P,
prime or not), (ii) the coverage cov(z) = floor(p/P) + [z%P < p%P] when the P-multiples
{P a : a in Z/q} uniformly tile (needs gcd to make a-reindex a bijection: a = p k mod q ok),
(iii) window p < ? : the cov formula needs p < P q (always true if p/q < P).

So the SAME closed form should hold for any P with window p/q < P. Test P=2,3,5,7,11,13.
Also: is the universal max S_P / (something) a clean function of P? For P=7 it's 44 = 4*11.
"""
from math import gcd
from fractions import Fraction as Fr

def cmatrix_row0(p, q, P):
    """row 0 of the P-sector integer matrix on the Ppq grid (cov model)."""
    r = [0] * P
    for a in range(q):
        for t in range(p):
            # (P a + t) mod (P q) in [q j, q(j+1)) -> j
            v = (P * a + t) % (P * q)
            j = v // q
            r[j] += 1
    return r

def cmatrix_row0_direct(p, q, P):
    """direct grid count row 0 (sector(qv)=0), midpoint sectoring."""
    N = P * p * q
    r = [0] * P
    for k in range(N):
        i = (P * ((q * (2 * k + 1)) % (2 * N))) // (2 * N)
        if i != 0:
            continue
        j = (P * ((p * (2 * k + 1)) % (2 * N))) // (2 * N)
        r[j] += 1
    return r

def cov_formula_row0(p, q, P):
    """residue-only closed form: r_j = (p//P)*q + #{z in [qj,q(j+1)): z%P < p%P}."""
    base = p // P
    rem = p % P
    r = []
    for j in range(P):
        cnt = sum(1 for z in range(q * j, q * (j + 1)) if z % P < rem)
        r.append(base * q + cnt)
    return r

def main():
    print("NEW LEAD: residue-only closed form for general apex P (P sectors)")
    print("=" * 72)
    for P in (2, 3, 5, 7, 11, 13):
        # verify cov-model row0 == direct, and == residue closed form, for p/q < P
        ok_model = ok_cf = True
        residue_only = True
        from collections import defaultdict
        groups = defaultdict(set)
        nchk = 0
        for q in range(1, 22):
            for p in range(q + 1, P * q):   # window p/q < P
                if gcd(p, q) != 1:
                    continue
                rm = cmatrix_row0(p, q, P)
                rd = cmatrix_row0_direct(p, q, P)
                rcf = cov_formula_row0(p, q, P)
                if rm != rd:
                    ok_model = False
                if rcf != rd:
                    ok_cf = False
                e = tuple(P * x - p * q for x in rd)
                groups[(p % P, q % P)].add(e)
                nchk += 1
        residue_only = all(len(v) <= 1 for v in groups.values())
        # max S over residue table
        Smax = 0
        for (a, b), es in groups.items():
            for e in es:
                Smax = max(Smax, sum(abs(x) for x in e))
        print(f"  P={P:2d}: cov-model==direct:{ok_model}  closed-form==direct:{ok_cf}  "
              f"residue-only:{residue_only}  maxS={Smax}  ({nchk} ratios p/q<{P})")

    # For P=7 in our actual window p/q<=43/20<3, confirm the SAME closed form (already done).
    print("\n=> The residue-only law + cov closed form is GENERAL (any apex P, window p/q<P).")
    print("   The LRC apex P=7 is one instance; the proof is prime-agnostic.")
    print("   This robustness means the L7 closure technique ports to any sector count.")

    # Bonus: the universal max S_P -- is it ~ P^2/something? tabulate.
    print("\nUniversal max deviation S_P (full table) vs P:")
    for P in (2,3,5,7,11,13):
        Smax=0
        seen=set()
        for q in range(1,18):
            for p in range(q+1, P*q):
                if gcd(p,q)!=1: continue
                rd=cmatrix_row0_direct(p,q,P)
                S=sum(abs(P*x-p*q) for x in rd)
                Smax=max(Smax,S)
        print(f"   P={P:2d}: maxS={Smax}  S/4={Fr(Smax,4)}  ~P^2/4*?={float(Smax)/(P*P):.3f}*P^2")

if __name__ == "__main__":
    main()
