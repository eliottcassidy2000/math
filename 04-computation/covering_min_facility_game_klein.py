#!/usr/bin/env python3
"""
covering_min_facility_game_klein.py  --  klein-2026-07-01-S85

THE COVERING-MIN AS AN ADVERSARIAL FACILITY-LOCATION GAME (owner framing + CS6840 game theory, S84).
  - DEFENDER (observer) picks a time t; payoff = min_{v in S} ||vt|| = distance to the nearest runner (facility).
  - Defender maximizes: M(S) = max_t min_v ||vt|| = the covering-min (defender's best loneliness given S).
  - ADVERSARY picks the config S (|S|=n-1); value = min_S M(S). LRC: this >= 1/n.
It is adversarial facility location / the covering radius: facilities (runners) 'cover' the circle over time;
the observer seeks the least-covered point; the adversary places facilities to minimize the max-empty-radius.

Computes for n=14 (finite reduction mod Phi6, S68):
  (1) the game value M(construction) = max_k min_v dist(v*k mod Phi6);
  (2) the LP-DUAL certificate = the equioscillation 2-point {runner 1, killer} (S73/HYP-3813);
  (3) the DISCREPANCY POTENTIAL (Koksma-Hlawka): star-discrepancy of the runner cloud + the coverage potential
      (union of danger arcs) -- the potential function for the inf-measure / PoA argument;
  (4) the POCHHAMMER fiber fraction f(n)=(1/2)_{n-2}/(n-2)! ~ 1/sqrt(pi*n) -- the pi/EVEN/MEASURE side,
      complementary to the sqrt(p)/ODD/CERTIFICATE side (HYP-3818/HYP-3820).
"""
from fractions import Fraction as Fr
import math

def dist_mod(r, m): r %= m; return min(r, m-r)

# ---------------------------------------------------------------
n = 14; Phi6 = n*n - n + 1                    # 183 = 3*61
S = list(range(1, n-1)) + [n*(n-1)]           # construction {1..12, 182}
print("="*74)
print(f"ADVERSARIAL FACILITY-LOCATION GAME, n={n}, Phi6={Phi6}=3*61; S={S}")
print("="*74)

# (1) game value M(S) = (1/Phi6) max_k min_v dist(v*k mod Phi6)
best_k, best_val = 0, -1
for k in range(1, Phi6):
    mv = min(dist_mod(v*k, Phi6) for v in S)
    if mv > best_val: best_val, best_k = mv, k
print(f"(1) GAME VALUE M(S) = max_k min_v dist(vk mod Phi6) = {best_val}/{Phi6} = {best_val/Phi6:.5f}")
print(f"    achieved at k* = {best_k}  (t* = {best_k}/{Phi6}); 1/n = {1/n:.5f}; M(S) {'>=' if best_val/Phi6>=1/n else '<'} 1/n")
print(f"    M(S) = n/Phi6 = {n}/{Phi6} = {n/Phi6:.5f} (covering-min, THM-523)  [defender's best loneliness]")

# (2) LP-DUAL certificate: the binding runners at k* (equioscillation)
resid = sorted(((dist_mod(v*best_k, Phi6), v) for v in S))
binding = [v for d0, v in resid if d0 == best_val]
print(f"\n(2) LP-DUAL / equioscillation certificate at k*: binding runners (dist=M) = {binding}")
signed = sorted(((v*best_k) % Phi6, v) for v in S)   # phase-residue cloud
near = [(r if r <= Phi6//2 else r-Phi6, v) for r, v in signed if dist_mod(r, Phi6) == best_val]
print(f"    binding phase-residues (signed, /Phi6): {[(f'{r:+d}', v) for r, v in near]}  "
      f"=> the observer(0) is PINNED between {near[0][0]:+d} and {near[-1][0]:+d} = +-n (Chebyshev 2-point, S73)")
print(f"    dual weights (1/2,1/2) on {{runner 1, killer}}: moving t either way brings one binder closer => M certified.")

# (3) DISCREPANCY POTENTIAL (Koksma-Hlawka) at t*
cloud = sorted((Fr(v*best_k % Phi6, Phi6)) for v in S)          # runner positions on [0,1)
# star discrepancy D* = max_x | #{points < x}/N - x |  (evaluate at the points)
Nn = len(cloud); Dstar = 0.0
for j, x in enumerate(cloud + [Fr(1)]):
    Dstar = max(Dstar, abs(j/Nn - float(x)), abs((j)/Nn - float(x)))
# three-gap check
gaps = sorted({float(cloud[(i+1) % Nn] - cloud[i]) % 1 for i in range(Nn)})
# coverage potential: each runner v covers a danger arc of half-width r around 0 over time; at fixed t, the
# 'potential' PHI = total covered length if each of N runners owns an arc of width 2*M -> N*2*M vs 1
PHI = Nn * 2 * (best_val/Phi6)
print(f"\n(3) DISCREPANCY POTENTIAL (Koksma-Hlawka |avg f(vt) - int f| <= V(f) D*):")
print(f"    runner cloud at t* has star-discrepancy D* = {Dstar:.4f}; three-gap sizes (approx) = {[round(g,4) for g in gaps]}")
print(f"    coverage potential PHI = N*2*M = {Nn}*2*{best_val}/{Phi6} = {PHI:.4f}  (<1 => an uncovered gap MUST exist)")
print(f"    => the potential-function/PoA floor: if total danger 2(n-1)M < 1 the observer escapes; equality tunes M.")
print(f"    Koksma-Hlawka reads the lonely measure L_C = int 1_lonely as a discrepancy of the {{vt}} cloud.")

# (4) POCHHAMMER fiber fraction  f(n) = (1/2)_{n-2} / (n-2)!  ~ 1/sqrt(pi n)   [pi/even/measure side]
def pochhammer_half(m):     # (1/2)_m = prod_{j=0}^{m-1} (1/2 + j)
    r = Fr(1)
    for j in range(m): r *= Fr(1, 2) + j
    return r
def fact(m):
    r = 1
    for j in range(2, m+1): r *= j
    return r
m = n-2
f = pochhammer_half(m) / fact(m)
print(f"\n(4) POCHHAMMER fiber fraction f(n)=(1/2)_(n-2)/(n-2)!:  f({n}) = {float(f):.4f}  "
      f"~ 1/sqrt(pi*n) = {1/math.sqrt(math.pi*n):.4f}")
print(f"    the pi/EVEN/MEASURE side (Wallis/Gamma, iota-even) -- complementary to the sqrt(p)/ODD/CERTIFICATE")
print(f"    side (Gauss i*sqrt(p), HYP-3820). Two halves of the apex: measure ~ 1/sqrt(pi n), certificate ~ sqrt(p).")

print("\n" + "="*74)
print("NET: covering-min = adversarial facility location; the LP-dual is the 2-point equioscillation {1,killer};")
print("the discrepancy/coverage potential (Koksma-Hlawka) is the potential-function/PoA handle on inf-measure;")
print("and the two apex halves are f(n)~1/sqrt(pi n) (measure) and sqrt(p) (certificate).")
print("="*74)
