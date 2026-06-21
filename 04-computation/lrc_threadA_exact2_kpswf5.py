#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_threadA_exact2_kpswf5.py  (kind-pasteur 2026-06-21, THREAD A part 4b -- CLOSED FORM, fixed scale)

f_P(a,b) = rowdef.  Observed (part 3):
   f_P(a,b) - 2 a b (P-a)(P-b)/P  =  g(||ab||_P)/P,   g(t)=2 t (P-t).
Equivalently:
   P * f_P(a,b) = 2 a b (P-a)(P-b) + 2 ||ab||_P (P - ||ab||_P).
i.e.  f_P(a,b) = [ 2 a b (P-a)(P-b) + 2 ||ab|| (P-||ab||) ] / P,   with ||x||=min(x%P,P-x%P).

Note 2ab(P-a)(P-b)/P is generally NOT an integer, nor is g/P; their SUM is.  Test this.
"""
from math import gcd
from fractions import Fraction as Fr


def rowdef_residue(a, b, P):
    for Qm in range(1, 80):
        q = b + P * Qm
        for Pm in range(1, 80):
            p = a + P * Pm
            if p <= q or p >= P * q:
                continue
            if gcd(p, q) != 1:
                continue
            rr = p % P
            rd = 0
            for j in range(P):
                Nj = sum(1 for z in range(q * j, q * (j + 1)) if z % P < rr)
                rd += abs(P * Nj - rr * q)
            return rd
    return None


def normP(x, P):
    m = x % P
    return min(m, P - m)


def main():
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    print("CLOSED FORM (fixed scale):")
    print("   f_P(a,b) = [ 2 a b (P-a)(P-b) + 2*||ab||*(P-||ab||) ] / P     (||x||=min(x%P,P-x%P))")
    print("=" * 80)
    allok = True
    for P in primes:
        h = (P - 1) // 2
        ok = True
        ex = []
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                rd = rowdef_residue(a, b, P)
                t = normP(a * b, P)
                num = 2 * a * b * (P - a) * (P - b) + 2 * t * (P - t)
                if Fr(num, P) != rd:
                    ok = False
                    if len(ex) < 3:
                        ex.append((a, b, rd, Fr(num, P)))
        print(f"   P={P:2d}: {'MATCH' if ok else 'NO '+str(ex)}")
        allok = allok and ok
    print(f"\n==> {'UNIFORM CLOSED FORM CONFIRMED (all primes 3..37)' if allok else 'FAILS'}")

    # Full S_P = P*f_P and the apex/sharp constants.
    print("\nFull discrepancy:  D_{p,q} = S_P(p%P,q%P)/(P p q),  S_P = P*f_P = 2ab(P-a)(P-b)+2||ab||(P-||ab||).")
    print("Sharp window constants for general P (window 1<p/q<P):")
    for P in primes:
        h = (P - 1) // 2
        # the LRC window is p/q in (1, ~2.15) for P=7; generally the sharp sup D*q etc.
        # max S_P / p over in-window residue reps -- but here just report maxS_P (full table max).
        maxS = 0
        argmax = None
        for a in range(1, P):
            for b in range(1, P):
                t = normP(a * b, P)
                aa, bb = normP(a, P), normP(b, P)
                S = 2 * aa * bb * (P - aa) * (P - bb) + 2 * t * (P - t)
                # careful: a,b enter via ||.|| in the area too
                if S > maxS:
                    maxS = S
                    argmax = (aa, bb, t)
        print(f"   P={P:2d}: maxS_P={maxS}  at (||a||,||b||,||ab||)={argmax}  (S/4={Fr(maxS,4)})")

    # Confirm against published general_prime maxS: 2,4,16,44,168,276 for P=2,3,5,7,11,13.
    print("\n   (published general_prime maxS: P=3->4, 5->16, 7->44, 11->168, 13->276)")


if __name__ == "__main__":
    main()
