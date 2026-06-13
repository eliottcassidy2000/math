#!/usr/bin/env python3
"""
lrc_repeated_addition_hyperoperation_tower_s548.py    oracle-2026-06-01-S548

Use 'multiplication = repeated addition' creatively + recursively.

Multiplication is repeated addition; its inverse face is log (turns x into +). That
single bridge unifies two LRC objects I built separately:
  * the MULTIPLICATIVE cascade  |SAFE| = prod_i c_i           (S545)
  * the ADDITIVE entropy / free energy  F = -sum_i log c_i    (S543)
  =>  |SAFE| = prod c_i = exp(-F),  i.e. the cascade IS the exponential of a REPEATED
      ADDITION of log-clearances. Multiplication-as-repeated-addition, made literal.

And the whole LRC thread climbs a HYPEROPERATION TOWER (each level = the previous
REPEATED):
  H0 succ (+1)      : the center = shift / odometer (S541)
  H1 add (+)        : runner positions s*t = t+...+t ; Goldbach even=p+q ; holdback sum|v_i-v_j|
  H2 mult (x = rep+): the cascade prod c_i ; doubling 2q=q+q ; channels mod (n/2) ; speed s = s*t
  H3 exp (^ = rep x): #tournaments 2^C(n,2) ; Burnside Fix=2^e ; covering (1-2/n)^{n-1} ; entropy 2^d
  H4 tetration      : metagraph (tournaments of tournaments) ; Cayley-Dickson dims 2^(2^..)
The DOUBLED PRIME 2q=q+q is the SEED of H2 (the simplest non-trivial repeated addition);
it propagates up (the (q,q) Burnside gcd=q at H3; the rank-one tree).

We (1) verify the bridge |SAFE| = prod c_i = exp(-F) and the additive build of F;
(2) tabulate the tower magnitudes; (3) show the doubled-prime seed.
"""
from fractions import Fraction
from math import log, exp
from functools import reduce
from math import gcd

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def safe_measure(runset, n):
    thr = Fraction(1, n); sp = (0,) + tuple(runset)
    W = {Fraction(0), Fraction(1)}
    for s in sp[1:]:
        if s == 0: continue
        for m in range(0, abs(s) + 1):
            for sgn in (1, -1):
                t = Fraction(n * m + sgn, n * abs(s))
                if 0 <= t <= 1: W.add(t)
    Wl = sorted(W); tot = Fraction(0)
    for a, b in zip(Wl, Wl[1:]):
        mid = (a + b) / 2
        if all(d0(Fraction(s) * mid) >= thr for s in sp[1:]): tot += (b - a)
    return tot

def cascade(speeds, n):
    runs = speeds[1:]; prev = Fraction(1); cs = []; Fpart = []; Facc = 0.0
    for i in range(len(runs)):
        cur = safe_measure(runs[:i + 1], n)
        c = cur / prev if prev > 0 else Fraction(0); cs.append(c); prev = cur
        Facc += (-log(float(c)) if c > 0 else float('inf')); Fpart.append(Facc)
    return cs, Fpart

def holdback(speeds):
    return sum(abs(a - b) for i, a in enumerate(speeds[1:]) for b in speeds[1 + i + 1:])

def main():
    print("Multiplication = repeated addition: the cascade<->free-energy bridge + the tower (oracle-S548)\n")

    print("(1) BRIDGE: |SAFE| = prod c_i = exp(-F),  F = -sum log c_i (repeated addition of log-clearances)")
    for name, sp, n in [("generic n=5",(0,1,3,5,7),5),("generic n=6",(0,1,2,4,7,8),6),("AP n=5 (tight)",(0,1,2,3,4),5)]:
        cs, Fp = cascade(sp, n)
        prod = reduce(lambda a,b: a*b, cs, Fraction(1))
        meas = safe_measure(sp[1:], n)
        Ftot = Fp[-1]
        chk = (prod == meas)
        recon = exp(-Ftot) if Ftot != float('inf') else 0.0
        print(f"   {name}: prod c_i = {prod} = |SAFE|? {chk};  F=-sum log c = {Ftot if Ftot!=float('inf') else 'inf'};  exp(-F)={recon:.5f} ~ {float(meas):.5f}")
        print(f"      additive build F_k (repeated addition): {[round(f,3) if f!=float('inf') else 'inf' for f in Fp]}")
    print("   => the multiplicative cascade and the additive free energy are ONE object via log/exp.\n")

    print("(2) THE HYPEROPERATION TOWER (each level = the previous REPEATED), magnitudes for n=6")
    sp,n=(0,1,2,4,7,8),6
    hb=holdback(sp); wind=sum(sp[1:]); meas=safe_measure(sp[1:],n)
    cov=(1-2/n)**(n-1); total=2**(n*(n-1)//2)
    print(f"   H1 add (+):    holdback sum|v_i-v_j| = {hb};  winding sum v_i = {wind}   (runner s*t = t added s times)")
    print(f"   H2 mult (x):   cascade prod c_i = |SAFE| = {meas} = {float(meas):.4f}   (= exp(-F), repeated add of logs)")
    print(f"   H3 exp (^):    #tournaments 2^C(n,2) = {total};  covering main (1-2/n)^(n-1) = {cov:.4f}")
    print(f"   H4 tetration:  metagraph (tournaments OF tournaments) / Cayley-Dickson dims 2^(2^..) -- the iterated count")
    print()

    print("(3) THE DOUBLED PRIME 2q=q+q SEEDS H2 and propagates up; n=2q Burnside (q,q): Fix=2^e, e has gcd(q,q)=q")
    for q in (3,5,7):
        n=2*q; lam=(q,q)
        e = sum((l-1)//2 for l in lam) + gcd(q,q)   # the all-odd Burnside exponent for (q,q)
        print(f"   n={n}=2*{q}: doubled prime; (q,q) Burnside exponent e = {(q-1)//2}*2 + gcd({q},{q})={(q-1)+q} = {e}; Fix=2^{e}; apex=q (S547)")
    print()
    print("READING: 'multiplication = repeated addition' is not a side identity here -- it is the")
    print("SPINE. The runners are additive (speed s = the unit step t repeated s times). Their")
    print("clearances MULTIPLY into the cascade (= exp of a repeated ADDITION of log-clearances,")
    print("so S545's product and S543's free energy are one object). The cascades EXPONENTIATE into")
    print("the counts (Burnside, covering, 2^C(n,2)). The counts TETRATE into the metagraph /")
    print("Cayley-Dickson tower. LRC is the cross-level claim: the ADDITIVE orbit (runners reaching)")
    print("is never fully blocked by the MULTIPLICATIVE cascade, measured by the EXPONENTIAL count.")
    print("The doubled prime 2q=q+q is the seed of the multiplicative rung, climbing as (q,q).")

if __name__ == "__main__":
    main()
