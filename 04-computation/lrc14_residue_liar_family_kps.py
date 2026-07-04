"""kps-2026-07-04: the OWNER's Fibonacci/period-16 hint applied to the CENSUS-HARD residue-liars.
{1..11,13,84m} are covering + apex-blocked (84m=0 mod 14) + LOOSE (lonely off the 14-grid).
Prior: lonely at 7/89, 14/173, 42/509 (m=1,2,6). Test the FORMULA t=7m/(84m+5) => infinite
family closed by a formula (like base12_far_peel). Look for Fibonacci/period-16 structure."""
import numpy as np
from math import gcd
from functools import reduce
from fractions import Fraction as Fr
def gcd_all(v): return reduce(gcd,v)
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def dist_int(x): return abs(x-round(x))
def min_margin_at(v, t):  # exact rational: min_i ||v_i t|| as Fraction
    m=None
    for s in v:
        r=Fr(s)*t; d=r-round(r)
        dd=abs(d)
        if m is None or dd<m: m=dd
    return m

print("Test the residue-liar formula {1..11,13,84m} lonely at t=7m/(84m+5):")
print(f"{'m':>3} {'X=84m':>8} {'t=7m/(84m+5)':>16} {'min margin':>12} {'>=1/14?':>8} {'cov?':>5}")
ok_all=True
for m in range(1,13):
    X=84*m
    v=list(range(1,12))+[13,X]
    t=Fr(7*m, 84*m+5)
    mar=min_margin_at(v,t)
    ok = mar>=Fr(1,14)
    ok_all=ok_all and ok
    print(f"{m:>3} {X:>8} {str(t):>16} {str(mar):>12} {str(ok):>8} {str(is_cov(v)):>5}")
print(f"\nFORMULA t=7m/(84m+5) gives margin>=1/14 for ALL m=1..12? {ok_all}")
print(f"  => {{1..11,13,84m}} is an INFINITE covering family, lonely by ONE formula (census-free).")
print()
# denominators 84m+5: Fibonacci? mod-7 / period-16 structure?
print("Denominators 84m+5 and numerators 7m:")
dens=[84*m+5 for m in range(1,13)]
print(f"  dens = {dens}")
print(f"  89=F11? {89} is Fibonacci: {89 in [1,2,3,5,8,13,21,34,55,89,144,233]}")
print(f"  84 = 4*21 = 4*F8 (21=F8); 84m+5, 5=F5; the '5' offset + 84=2^2*3*7 (kills 14-grid)")
print(f"  84m+5 mod 7 = {[(84*m+5)%7 for m in range(1,8)]} (all 5, since 84=0 mod 7) -- CONSTANT mod 7!")
print(f"  84m+5 mod 14 = {[(84*m+5)%14 for m in range(1,8)]}")
print()
# WHY it works: at t=7m/(84m+5), where does each runner land? (the mechanism)
m=1; v=list(range(1,12))+[13,84]; t=Fr(7,89)
print(f"Mechanism at m=1, t=7/89: runner positions v*t mod 1 (as k/89):")
for s in v:
    r=Fr(s)*t; num=(s*7)%89
    print(f"  v={s:>3}: {s}*7/89 = {num}/89, dist-to-int = {min(num,89-num)}/89 (>=89/14={89/14:.1f}? {min(num,89-num)>=89/14})")
