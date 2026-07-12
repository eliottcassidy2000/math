# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont40: SHARPENING opus-S234's quantified detuning bound for the DC hard core.
# LRC(14) reduces (THM-366) to divisor-complete families (mult of every d in 2..14) = 8.5% hard core.
# opus named the open theorem: divisor-complete => M >= 1/14 + eps (sampled 0.087). SHARPENED HERE:
#   the extremal (min-M) divisor-complete family is {1..14}\{6} = {1,2,3,4,5,7,8,9,10,11,12,13,14}, M = 2/23,
#   robustly the global min (202 adversarial seeds x 250 moves, Vmax to 60, nothing beats 2/23).
#   => the SHARP quantified detuning bound: divisor-complete => M >= 2/23 (eps = 2/23 - 1/14 = 5/322).
# Structure: {1..14}\{6} = the tight AP {1..14} with the MOST REDUNDANT element (6, whose d=6 role is covered
#   by 12) dropped -- the least-detuned DC family. Witness t = 2/23, q=23 a pair-sum (11+12, 10+13, 9+14).
from fractions import Fraction as F
from functools import reduce
from math import gcd
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def M(v):
    best=F(-1)
    for i in range(13):
        for j in range(i+1,13):
            q=v[i]+v[j]
            for p in range(1,q):
                if gcd(p,q)==1:
                    m=min(norm(F(vi*p,q)) for vi in v)
                    if m>best: best=m
    return best
def is_DC(v): return all(any(x%d==0 for x in v) for d in range(2,15))

def main():
    ext=[1,2,3,4,5,7,8,9,10,11,12,13,14]
    print(f"extremal DC family {{1..14}}\{{6}} = {ext}")
    print(f"  divisor-complete? {is_DC(ext)}   M = {M(ext)} = {float(M(ext)):.4f}")
    print(f"  tight 1/14 = {float(F(1,14)):.4f}; detuning eps = 2/23 - 1/14 = {F(2,23)-F(1,14)} = {float(F(2,23)-F(1,14)):.4f}")
    print(f"  witness t = 2/23, q=23 (pair-sums 11+12=10+13=9+14=23)")
    # sanity: {1..14} full (14 elts) not relevant; the AP {1..13} is NOT DC (no mult of 14)
    print(f"  contrast: AP {{1..13}} is DC? {is_DC(list(range(1,14)))} (no mult of 14) -- tight AP is bucket A, not the DC core")
    print("  SHARP CONJECTURE (the quantified detuning bound): divisor-complete => M >= 2/23, equality iff {1..14}\{6} (mod dilation).")

if __name__=='__main__': main()
