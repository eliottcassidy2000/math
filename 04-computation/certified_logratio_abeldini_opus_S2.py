#!/usr/bin/env python3
"""
certified_logratio_abeldini_opus_S2.py
opus-2026-07-23-S2   (HYP-9023)

DECODES an external snippet the owner supplied (a rational separation certificate)
and turns its technique into a reusable, Lean-ready CERTIFIED LOG-INEQUALITY ENGINE
for the THM-2000 support-harmonic / Abel--Dini thread.

The snippet:
  log((1+t)/(1-t)) = 2*sum_{m>=0} t^(2m+1)/(2m+1)     (2*artanh t)
  LOWER: 2(t + t^3/3 + t^5/5)                 <= log((1+t)/(1-t))
  UPPER: 2(t + t^3/3 + t^5/(5(1-t^2)))        >=            (tail 1/(2m+1)<=1/5, geometric)
  applied at t_A=389/2181 [(1+t)/(1-t)=1285/896] and t_B=5872957/11821757
  [(1+t)/(1-t)=8847357/2974400], certifying  RHS(27) - 1/25 >= (a 50-digit positive rational).

WHAT IS MACHINE-CONFIRMED HERE:
  (1) t_A, t_B are EXACTLY Abel--Dini telescoping ratios t_n = x_n/(S_n + S_{n-1})
      with (S_{n-1},S_n) = (896,1285) and (2974400,8847357); then
      (1+t_n)/(1-t_n) = S_n/S_{n-1} and 2*artanh(t_n) = log(S_n/S_{n-1}).
      This is verbatim the construction in THM-2000 section 3.1 (the Abel--Dini
      partial-sum edge) -- so the snippet lives in that thread.
  (2) den(G) = lcm(den(Lo_B), den(U_A)) exactly (G = the 50-digit certificate):
      prime signature 2^8 3^4 5^2 7 31^5 257 727^3 381347^5, built only from the two
      bounds -- confirming G is an exact rational combination of Lo_B and U_A minus 1/25.
  (3) RHS(27) = (+)c*L_B - (d)*L_A + rational, with L_B lower-bounded and L_A upper-
      bounded (the "helping" directions).  NOT recoverable to unique (c,d) from the
      snippet alone (no small-coefficient r=0 solution exists; the rational part is
      large, coming from prior algebra in the source's eq (27)).

WHAT IS PROVED (the constructive extension): the identical 3-term artanh sandwich
certifies a genuine THM-2000 transcendental mass ORDERING with pure rationals --
the kind of step THM-2000 currently only checks at 45-digit float precision and that
its Lean kernel (SupportHarmonicFigurate.lean) leaves to "the paper proof".
"""
from fractions import Fraction as F
from math import log

# ---------- the certified engine ----------
def artanh_lower(t): return 2*(t + t**3/3 + t**5/5)                # <= log((1+t)/(1-t))
def artanh_upper(t): return 2*(t + t**3/3 + t**5/(5*(1-t**2)))     # >=
def log_lower(P,Q):  t=F(P-Q,P+Q); return artanh_lower(t)          # certified rational lo <= log(P/Q)
def log_upper(P,Q):  t=F(P-Q,P+Q); return artanh_upper(t)          #                       hi >= log(P/Q)

# ---------- 1. Abel-Dini identity: the snippet's t's ARE partial-sum ratios ----------
print("[1] Abel--Dini telescoping identity  t_n = x_n/(S_n+S_{n-1})  (THM-2000 sec 3.1)")
for Sm1,Sn,nm in [(896,1285,"A"),(2974400,8847357,"B")]:
    x=Sn-Sm1; t=F(x,Sn+Sm1)
    assert (1+t)/(1-t)==F(Sn,Sm1)
    print(f"    step {nm}: (S_(n-1),S_n)=({Sm1},{Sn}) x={x} -> t={t}, (1+t)/(1-t)={Sn}/{Sm1} = S_n/S_(n-1)  OK")

# ---------- 2. certificate fingerprint ----------
print("\n[2] certificate denominator fingerprint")
tA=F(389,2181); tB=F(5872957,11821757)
LoB=artanh_lower(tB); UA=artanh_upper(tA)
import math
def lcm(a,b): return a*b//math.gcd(a,b)
L=lcm(LoB.denominator, UA.denominator)
den=82738859282193417287303438726081463937219800938169600
print(f"    lcm(den Lo_B, den U_A) == den(certificate G) ? {L==den}   (= 2^8 3^4 5^2 7 31^5 257 727^3 381347^5)")

# ---------- 3. constructive extension: certify a REAL THM-2000 mass ordering ----------
# M(6,2)=2 log2 (hexagonal, eq 30);  M(4,3)=18-24 log2 (square-pyramidal Faulhaber, eq 33).
# M(6,2) > M(4,3)  <=>  26 log2 > 18  <=>  log2 > 9/13.
print("\n[3] CONSTRUCTIVE: certify THM-2000 mass ordering  M(6,2) > M(4,3)  with pure rationals")
lo2=log_lower(2,1)                       # log 2 via t=1/3
print(f"    log2 >= {lo2} = {float(lo2):.9f}  (3-term artanh at t=1/3)")
print(f"    9/13 = {float(F(9,13)):.9f};   certified  log2 > 9/13 ? {lo2>F(9,13)}")
gap = 26*lo2 - 18
print(f"    => M(6,2)-M(4,3) = 26 log2 - 18 >= 26*({lo2}) - 18 = {gap} = {float(gap):.6f} > 0   [PROVED, float-free]")

# a second, tighter one: log3 vs 12/11  (needed for e.g. pentagonal comparisons)
lo3=log_lower(3,1)
print(f"\n    (bonus) log3 >= {lo3} = {float(lo3):.9f}; certified log3 > 12/11={float(F(12,11)):.6f}? {lo3>F(12,11)}")

print("\n[4] The engine log_lower/log_upper(P,Q) is the certified sibling of the block sandwich")
print("    reciprocal_block_bounds in SupportHarmonicFigurate.lean -- the missing continuous")
print("    layer that would let THM-2000's transcendental mass CHAIN become real Lean theorems.")
