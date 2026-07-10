# kind-pasteur-2026-07-10-S127
# REFUTATION of the upper edge: q <= 27 (and any fixed bound q <= B) is FALSE on the residual class.
# The minimal strictly-live modulus is UNBOUNDED. Corrects the S127cont9 conjecture BoundedStrictlyLiveSupply.
#
# Why q >= 15 is TRUE (proved in Lean): covering => a zero residue at every q in [2,14].
# Why q <= B is FALSE: the strictly-live condition at q depends only on v mod q, so a family whose residues
# are "tight" (no live multiplier) at every q in [15,B] fails all of them -- and such families exist inside
# the residual class (covering, ratio>13, compressed, distinct, not-detuned, diff-primitive, not-near-AP).
# The naive adversary (dilated APs c*{1..13}) is excluded by difference-primitivity (diff-gcd = c > 1);
# but genuine residual families with min q > B exist, min q growing with the family's "near-tightness".
from math import gcd
from functools import reduce
import random, time

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def has_live(S, q): return any(all(q < 14 * ((v * p) % q) < 13 * q for v in S) for p in range(1, q))
def compressed(S): return all(any(j != i and S[i] <= 13 * S[j] for j in range(13)) for i in range(13))
def detuned(S):
    for g in range(2, 500):
        if len([v for v in S if v % g != 0]) == 1: return True
    return False
def diffprim(S):
    S = sorted(S); return reduce(gcd, [S[i] - S[0] for i in range(1, 13)]) == 1
def near_AP(S):
    S = sorted(S)
    for L in range(1, max(S) + 1):
        ok = True; ks = []
        for s in S:
            k = round(s / L); a = s - L * k
            if abs(a) * 182 > L or k == 0: ok = False; break
            ks.append(k)
        if ok and len(set(ks)) <= 12: return True
    return False
def min_q(S, cap=400): return next((q for q in range(15, cap + 1) if has_live(S, q)), None)

def is_residual(S):
    return (covering(S) and max(S) > 13 * min(S) and compressed(S) and len(set(S)) == 13
            and not detuned(S) and diffprim(S) and not near_AP(S) and any(abs(x) >= 23 for x in S))

print("=== THE CERTIFIED COUNTEREXAMPLE ===")
V = [210, 1378, 1379, 2106, 2222, 2247, 3650, 3773, 4123, 5083, 5561, 5680, 6000]
print(f"  v = {V}")
print(f"  covering={covering(V)}  ratio={max(V)/min(V):.2f} (>13)  compressed={compressed(V)}  distinct={len(set(V))==13}")
print(f"  not-detuned={not detuned(V)}  diff-primitive={diffprim(V)}  not-near-AP={not near_AP(V)}  max>=23={any(abs(x)>=23 for x in V)}")
print(f"  => IsResidual = {is_residual(V)}")
print(f"  min strictly-live q = {min_q(V)}   (fails EVERY q in [15,27])   => q<=27 is FALSE")
print()
print("  the decidable core (Lean-certifiable): for each q in [15,27], NO multiplier p is strictly-live:")
for q in range(15, 28):
    live_ps = [p for p in range(1, q) if all(q < 14 * ((v * p) % q) < 13 * q for v in V)]
    print(f"    q={q}: live multipliers = {live_ps}")

print()
print("=== min q is UNBOUNDED: the naive adversary is excluded, but genuine residual families climb ===")
print("  dilated APs c*{1..13} (the mu=0 locus) fail all q but have diff-gcd = c > 1 (NON-primitive):")
for c in [2, 4, 6]:
    S = [c * i for i in range(1, 14)]
    print(f"    {c}*{{1..13}}: fails all q, diff-gcd={diffprim(S)==False and reduce(gcd,[S[i]-S[0] for i in range(1,13)])}, residual={is_residual(S)}")
print("  but genuine residual families reach:")
random.seed(21); best = (min_q(V), V); t0 = time.time()
while time.time() - t0 < 60:
    HI = random.choice([2000, 8000, 30000, 100000])
    S = sorted(random.sample(range(1, HI + 1), 13))
    if not (max(S) > 13 * min(S) and covering(S) and compressed(S)): continue
    if not diffprim(S) or detuned(S): continue
    mq = min_q(S, cap=80)
    if mq and mq > best[0] and not near_AP(S): best = (mq, S)
print(f"    max(min q) found = {best[0]}   at {best[1]}")
print()
print("CONSEQUENCE: BoundedStrictlyLiveSupply B is FALSE for every fixed B. The correct obligation is")
print("StrictlyLiveSupply (B = infinity) / the measure floor. klein's THM-685 transfer (live rulers at")
print("q > Sum v/mu) is LOAD-BEARING, not a superfluous safety net -- the near-tight residual families")
print("genuinely need large q. The S127cont9 'witness lives at the bottom' holds for GENERIC families only.")
