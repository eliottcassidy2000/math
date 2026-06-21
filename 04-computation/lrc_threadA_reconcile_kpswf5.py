#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_threadA_reconcile_kpswf5.py  (kind-pasteur 2026-06-21, THREAD A -- naming reconciliation)

Make ABSOLUTELY sure our closed form reproduces the canon HYP-2739 numbers, with the
exact same meaning of every symbol.  Canon (HYP-2739 / HYP-2742):

  D_{p,q} = S(p%7,q%7)/(7pq),  S = 4 f(||p||_7,||q||_7),
  ||x||_7 in {0,1,2,3},
  f(a,b) = 0 if ab=0 / ab+3 if a!=b / ab+4-2|a-2| if a=b.
  HYP-2742 S-table (this is the '4f' object):
    S=[[0,0,0,0,0,0,0],[0,12,20,24,24,20,12],[0,20,32,36,36,32,20],
       [0,24,36,44,44,36,24],[0,24,36,44,44,36,24],[0,20,32,36,36,32,20],
       [0,12,20,24,24,20,12]]   (indexed by p%7,q%7 = 0..6)

So in CANON: the numerator of D is 'S = 4f', S(3,3)=44 is its max, D=S/(7pq).

Our rowdef object = sum_j |7 c_0j - pq| = the SAME 4f (HYP-2739 step says
'S = sum_j |e_j| = rowdef, full-matrix defect = 7S').  So:
   rowdef = 4f = canon S (max 44).   GOOD.
And our 'closed_form_num' = 2ab(P-a)(P-b)+2||ab||(P-||ab||) over /P gives rowdef? Let's
check: part4b matched rowdef = [2ab(P-a)(P-b)+2||ab||(P-||ab||)]/P.  For P=7,(3,3):
   2*9*16 + 2*||9||*(7-||9||) = 288 + 2*2*5 = 288+20=308; /7 = 44.  YES = canon S(3,3).

So the FINAL canon-aligned statement is:
   canon-S(a,b) = 4 f = [ 2 a b (P-a)(P-b) + 2 ||ab||_P (P-||ab||_P) ] / P,
                  with a=||p||_P, b=||q||_P  (a,b in 0..(P-1)/2),
   and D_{p,q} = canon-S / (P p q).
This is the uniform general-P law.  Here we VERIFY it reproduces the literal HYP-2742
7x7 table over the FULL residue range p%7,q%7 in 0..6 (using a=||p||,b=||q||).
"""
from fractions import Fraction as Fr


def normP(x, P):
    m = x % P
    return min(m, P - m)


def canonS(pres, qres, P):
    """canon S = 4f as a function of residues p%P,q%P via a=||p||,b=||q||."""
    a = normP(pres, P)
    b = normP(qres, P)
    if a == 0 or b == 0:
        return 0
    t = normP(a * b, P)   # ||ab||; note ||a||*||b|| then norm again
    num = 2 * a * b * (P - a) * (P - b) + 2 * t * (P - t)
    assert num % P == 0, (pres, qres, P, num)
    return num // P


def main():
    P = 7
    canon = [[0,0,0,0,0,0,0],
             [0,12,20,24,24,20,12],
             [0,20,32,36,36,32,20],
             [0,24,36,44,44,36,24],
             [0,24,36,44,44,36,24],
             [0,20,32,36,36,32,20],
             [0,12,20,24,24,20,12]]
    print("Reproduce literal HYP-2742 7x7 S-table from the uniform closed form:")
    ok = True
    for pr in range(7):
        row = []
        for qr in range(7):
            v = canonS(pr, qr, P)
            row.append(v)
            if v != canon[pr][qr]:
                ok = False
        print("   " + " ".join(f"{x:3d}" for x in row))
    print(f"   identical to canon HYP-2742 table: {ok}")

    # also reproduce the published P=7 f (rowdef/4): f = S/4.
    print("\n  published f = S/4 (a,b in 0..3):")
    for a in range(4):
        row = []
        for b in range(4):
            S = canonS(a, b, 7)
            row.append(Fr(S, 4))
        print("   " + " ".join(f"{str(x):>4s}" for x in row))
    print("   should be f(a,b): 0;  3,5,6; 5,8,9; 6,9,11 (the ab+3 / ab+4-2|a-2| values)")

    # general-P sharp constants in the ACTUAL LRC-type window.
    # For LRC the sharp face was D*q <= 12/7 at p/q=3/2 (P=7).  The closed form gives, for the
    # in-window minimal-denominator ratios, an exact face.  Report sup over coprime in-window.
    print("\nGeneral-P: the apex law and max are clean:")
    print("   canon-S(a,b)=0  <=>  P|p or P|q   (a=0 or b=0): the APEX law, all P.")
    for P in (3,5,7,11,13,17,19,23,29,31,37,41,43):
        h=(P-1)//2
        m=max(canonS(a,b,P) for a in range(P) for b in range(P))
        print(f"   P={P:2d}: max canon-S = {m}  (S/4={Fr(m,4)}); = [2h^2(P-h)^2 + 2||h^2||(P-||h^2||)]/P")


if __name__ == "__main__":
    main()
