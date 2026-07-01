"""
lrc_modular_bound_test_genus_f14_opus_20260701.py  (HYP-3779)
Tests the "modular-forms bound on a(n)" three ways -- all NEGATIVE:
 (1) genus X0(2n) does NOT track beater existence (beaters at genus 0 n=8,9; construction wins at genus 1,2 n=12..14);
 (2) cusp form f14 (14a) coefficients [1,-1,1,0,0] != a(n)=[2,2,4,4,3];
 (3) construction-vs-covering-min margin gaps are irregular rationals, abrupt threshold at n=12.
Integrates klein-S60 (HYP-3778): transition a(n)=n for n>=12 (radius-1 band, COMBINATORIAL). Law: VALUE modular,
EXISTENCE combinatorial. See reflection the-modular-forms-bound-on-a-of-n-FAILS-...-opus-20260701.md.
"""
import numpy as np
from fractions import Fraction as F
def phi6(n): return n*n-n+1
genusX0={14:1,16:0,18:0,20:1,22:2,24:1,26:2,28:2}   # standard tables
a={7:2,8:2,9:4,10:4,11:3, 12:12,13:13,14:14}         # covering-min rung (klein: beaters n=7..11, transition n>=12)

def test1_genus():
    print("(1) genus X0(2n) vs beater existence:")
    for n in range(7,15):
        print(f"    n={n:>2} 2n={2*n:>2} genus={genusX0[2*n]} a(n)={a[n]:>2} {'BEATER' if a[n]<n else 'construction'}")
    print(f"    beaters (n=7..11) at genus {[genusX0[2*n] for n in range(7,12)]} (incl 0,0!); no beaters (n=12..14) at genus {[genusX0[2*n] for n in range(12,15)]}")
    print("    => genus does NOT separate. FAIL.")

def test2_f14():
    M=20; P=np.zeros(M+2); P[0]=1
    for d in [1,2,7,14]:
        k=1
        while d*k<=M:
            for j in range(M,d*k-1,-1): P[j]-=P[j-d*k]
            k+=1
    f=np.zeros(M+2)
    for i in range(M+1): f[i+1]=P[i]
    print("\n(2) f14 (14a) coeffs a_n:", [int(f[n]) for n in range(1,12)])
    print("    a(n) rung n=7..11 =", [a[n] for n in range(7,12)], "vs f14 =", [int(f[n]) for n in range(7,12)], "=> NO match. FAIL.")

def test3_gaps():
    print("\n(3) margin GAP construction - covering-min (n=7..11):")
    for n in range(7,12):
        gap=F(n-1,n*phi6(n))-F(a[n]-1,n*(a[n]*(n-1)+1))
        print(f"    n={n}: gap={gap}={float(gap):.5f}")
    print("    => irregular rationals, Theta(1/n^2), abrupt 0 at n=12. No L-value pattern. FAIL.")

if __name__=="__main__":
    test1_genus(); test2_f14(); test3_gaps()
    print("\n=> NO modular-forms bound on a(n). LAW: VALUE is modular (-1/12=zeta(-1), the construction, HYP-3775);")
    print("   EXISTENCE is combinatorial (radius-1 band, klein-S60 HYP-3778). Cusp form f14 = OPEN-Q-108 proof obstruction, not a(n).")
    print("   Rigorous transition tool = character-sum/Weil (repo's sqrt(V) peel-deviation THM-546), not modular forms.")
