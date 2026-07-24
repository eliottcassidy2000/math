#!/usr/bin/env python3
"""
snippet_eq27_reconcile_opus_S4.py   opus-2026-07-23-S4

Reconcile opus-S3 ("c unpinnable") with mac-mini-S168 ("c=2457/6592 exact").
Independently verify the EXACT identity, test the uniqueness-of-low-height claim,
factor every key integer (smoothness = the diagnostic: Baker/log-prime-lattice vs
generic Abel-Dini), and locate the two logs relative to THM-2000 figurate masses.

Snippet object:  RHS(27) - 1/25 >= G  (G the 51/53-digit rational), certified by
  2(t+t^3/3+t^5/5) <= log((1+t)/(1-t)) <= 2(t+t^3/3+t^5/(5(1-t^2)))
at t_A=389/2181 (upper bound => U_A >= log_A) and t_B=5872957/11821757 (lower => L_B <= log_B).
"""
from fractions import Fraction as F
import math
from sympy import factorint
import mpmath as mp
mp.mp.dps = 90

tA = F(389, 2181)
tB = F(5872957, 11821757)

# the two ratios (1+t)/(1-t)
AA = (1 + tA) / (1 - tA)
BB = (1 + tB) / (1 - tB)
print("ratio_A (1+tA)/(1-tA) =", AA, " expected 1285/896  ->", AA == F(1285, 896))
print("ratio_B (1+tB)/(1-tB) =", BB, " expected 8847357/2974400 ->", BB == F(8847357, 2974400))

# Abel-Dini telescoping check: t = x/(Sn+Sn1), (Sn1,Sn) consecutive partial sums
for tag, Sn1, Sn in [("A", 896, 1285), ("B", 2974400, 8847357)]:
    x = Sn - Sn1
    t = F(x, Sn + Sn1)
    print(f"  [{tag}] x={x}, Sn-1+Sn={Sn+Sn1}, t=x/(Sn1+Sn)={t}  matches t_{tag}:",
          t == (tA if tag == "A" else tB))

# series bounds
def lower(t): return 2 * (t + t**3 / 3 + t**5 / 5)                    # <= log
def upper(t): return 2 * (t + t**3 / 3 + t**5 / (5 * (1 - t**2)))     # >= log
UA = upper(tA)   # >= log_A
LB = lower(tB)   # <= log_B

Gnum = 391926968594914200867482400554891567498742649630277
Gden = 82738859282193417287303438726081463937219800938169600
G = F(Gnum, Gden)

# ============ CLAIM 1: mac-mini exact identity ============
c = F(2457, 6592)
lhs = c * LB - UA - F(1, 25)
print("\n[CLAIM] (2457/6592)*L_B - U_A - 1/25 == G exactly :", lhs == G)
print("        residual (should be 0):", lhs - G)

# integer form
print("[CLAIM] integer form 2457*L_B - 6592*U_A >= 6592/25 :",
      2457 * LB - 6592 * UA, ">=", F(6592, 25), "->", 2457 * LB - 6592 * UA >= F(6592, 25))

# ============ CLAIM 2: uniqueness / low-height signal ============
# Target T = 1/25 + G = c*L_B - U_A  (assume r=0). For each integer d on U_A,
# c_d=(T + d*U_A)/L_B; the TRUE decode is the d giving anomalously low height.
T = F(1, 25) + G
print("\n--- height of c_d=(T+d*U_A)/L_B over integer d (r=0 model) ---")
for d in range(-2, 6):
    cd = (T + d * UA) / LB
    h = max(abs(cd.numerator), abs(cd.denominator))
    shown = str(cd) if h < 10**8 else "(height ~10^%.0f)" % math.log10(h)
    print(f"  d={d:+d}: log10(height)={math.log10(h):5.1f}   c={shown}")

# ============ factorizations (smoothness diagnostic) ============
print("\n--- factorizations ---")
keys = [("Sn-1(A)=896", 896), ("Sn(A)=1285", 1285), ("x(A)=389", 389),
        ("Sn1+Sn(A)=2181", 2181),
        ("Sn-1(B)=2974400", 2974400), ("Sn(B)=8847357", 8847357),
        ("x(B)=5872957", 5872957), ("Sn1+Sn(B)=11821757", 11821757),
        ("coef num 2457", 2457), ("coef den 6592", 6592),
        ("G num", Gnum), ("G den", Gden)]
for name, val in keys:
    print(f"  {name:22s} = {dict(factorint(val))}")

# ============ true values & margins ============
logA = mp.log(mp.mpf(1285) / 896)
logB = mp.log(mp.mpf(8847357) / 2974400)
val = mp.mpf(2457) / 6592 * logB - logA
print("\n--- true values ---")
print("  log_A =", mp.nstr(logA, 30))
print("  log_B =", mp.nstr(logB, 30))
print("  log_B / log_A =", mp.nstr(logB / logA, 20))
print("  (2457/6592)log_B - log_A =", mp.nstr(val, 30))
print("  1/25 =", mp.nstr(mp.mpf(1) / 25, 20), "   margin =", mp.nstr(val - mp.mpf(1) / 25, 20))
print("  certified G =", mp.nstr(mp.mpf(Gnum) / Gden, 20))

# ============ relate logs to THM-2000 masses / log primes ============
print("\n--- log-prime-lattice test (Baker reading) ---")
# log_A is EXACTLY an integer combo of log primes iff 1285,896 are smooth (they are):
# 1285=5*257, 896=2^7*7  => log_A = log5 + log257 - 7log2 - log7
la = mp.log(5) + mp.log(257) - 7 * mp.log(2) - mp.log(7)
print("  log_A - (log5+log257-7log2-log7) =", mp.nstr(logA - la, 10), "(0 => exact prime combo)")
# log_B exact prime combo depends on factoring 8847357:
fB = factorint(8847357)
print("  8847357 factored:", dict(fB), " -> largest prime =", max(fB))
fBd = factorint(2974400)
print("  2974400 factored:", dict(fBd))

# THM-2000 masses near these logs?
log2 = mp.log(2)
print("\n--- THM-2000 mass proximity ---")
print("  M(6,2)=2log2 =", mp.nstr(2 * log2, 15))
print("  M(4,3)=18-24log2 =", mp.nstr(18 - 24 * log2, 15))
print("  M(6,2)-M(4,3) = 26log2-18 =", mp.nstr(26 * log2 - 18, 15))
print("  log_A =", mp.nstr(logA, 15), "  ~ M(6,2)-M(4,3)? ratio:", mp.nstr(logA / (26 * log2 - 18), 8))

# PSLQ: is [log_A, log_B, 1] rationally dependent with small coeffs? (already know exact via primes)
rel = mp.pslq([logA, logB], maxcoeff=10**6, maxsteps=10**5)
print("\n  pslq[log_A, log_B] small integer relation:", rel, "(None => Q-independent at that height)")
