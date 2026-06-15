#!/usr/bin/env python3
"""
lrc14_singular_series_adelic_macmini_0614s1.py  (mac-mini-2026-06-14-S1)

THE LRC(14) SINGULAR SERIES L(S), pushed on three creative fronts (THM-501 frontier):

  C'(14)  <==  inf_S L(S) > 0,   L(S) = lim_{q->inf} D(q,S)/q,
  D(q,S) = #{a in Z/q : v*a mod q not in B_q for all v in S},  B_q = +-{0..floor(q/14)}.

FRONT 1 (adversarial infimum): push far past THM-501's 120 samples. Structured search
  over dilated-AP cores d*{1..12} u {r} (the identified extremizers) + random high-energy
  configs, to pin inf L and confirm inf L > 0 with margin.
FRONT 2 (the 7-vanishing): the sinc s(t)=sin(pi t/7)/(pi t) VANISHES at t in 7Z. So in
  L = (6/7)^13 + sum_{exact relations} (6/7)^{13-|T|}(-1)^|T| prod s(t_v), every term whose
  relation has ANY coefficient t_v divisible by 7 drops out. Verify the structural effect:
  speed sets all ≡ 0 mod 7 vs mixed; relations forced coprime-to-7.
FRONT 3 (adelic / Euler-product test, HYP-2503): is L = beta_inf * prod_p beta_p?
  Compute D(p^e,S)/p^e along single primes p=2,3,5,7,11,13 as e grows; if L localizes,
  the product of single-prime limits should reconstruct L. Hypothesis to TEST (and likely
  refute): L is purely ARCHIMEDEAN (the exact-relation sum), p-adic data lives in the
  APPROACH (threshold q*), not in L -> single-prime limits each return the SAME (6/7)^13-ish
  value (no genuine factorization).
"""
import sys, math, random
from itertools import combinations

sys.stdout.reconfigure(line_buffering=True)
random.seed(140614)

def deficit(q, S):
    """D(q,S) = # a in [0,q) with v*a mod q outside B_q for ALL v. B_q radius = floor(q/14)."""
    rad = q // 14
    # a is 'dangerous-free' iff every v*a mod q has min(r,q-r) > rad
    cnt = 0
    for a in range(q):
        safe = True
        for v in S:
            r = (v * a) % q
            if r <= rad or r >= q - rad:
                safe = False
                break
        if safe:
            cnt += 1
    return cnt

def L_estimate(S, qs=(13999, 14000, 14001, 15013, 15015)):
    """Empirical singular series: average D(q,S)/q over several large q (window-averaged)."""
    vals = [deficit(q, S) / q for q in qs]
    return sum(vals) / len(vals), min(vals), max(vals)

def gcd_all(xs):
    from math import gcd
    g = 0
    for x in xs: g = gcd(g, x)
    return g

def is_Cprime_config(S):
    """primitive (gcd=1), 13 distinct positives, contains a multiple of 14."""
    return len(set(S)) == 13 and gcd_all(S) == 1 and any(v % 14 == 0 for v in S)

MAIN = (6/7) ** 13

print("=" * 74)
print(f"LRC(14) SINGULAR SERIES L(S).  Main term (6/7)^13 = {MAIN:.6f}")
print("=" * 74)

# ---- sanity: tight (non-primitive) L=0 vs primitive generic L ~ (6/7)^13 ----
tight = [14*i for i in range(1, 14)]   # 14*{1..13}: NON-primitive (gcd 14) -> excluded from C'
Lt, lo, hi = L_estimate(tight)
print(f"\nTIGHT 14*{{1..13}} (gcd=14, NOT a C' config): L ~ {Lt:.5f}  (D~0)  [{lo:.4f},{hi:.4f}]")
# primitive generic set with a multiple of 14, near-Sidon (geometric): 14 plus distinct geometrics
gen = sorted(set([14] + [2**k + 1 for k in range(1, 13)]))[:13]
assert is_Cprime_config(gen), gen
Lg, lo, hi = L_estimate(gen)
print(f"GENERIC primitive+mult14 (near-Sidon): L ~ {Lg:.5f}  (expect ~{MAIN:.4f})  [{lo:.4f},{hi:.4f}]  set={gen}")

# ---- FRONT 1: adversarial infimum over PRIMITIVE multiple-of-14 evader cores ----
print("\n" + "=" * 74)
print("FRONT 1: adversarial infimum over PRIMITIVE mult-of-14 configs (evader cores 7*{1..12}u{r})")
print("=" * 74)
results = []
# evader family: 7*{1..12} (contains 14=7*2) u {r}, r coprime to 7 -> primitive
core7 = [7*i for i in range(1, 13)]   # {7,14,...,84}, contains 14
for r in list(range(15, 300)):
    if r % 7 == 0 or r in core7:
        continue
    S = sorted(set(core7 + [r]))
    if not is_Cprime_config(S):
        continue
    L, lo, hi = L_estimate(S, qs=(14000, 15013))
    results.append((L, r, tuple(S)))
# dilated cores d*{1..12} u {14*m}: primitive iff gcd(d,14m)=1; contains 14*m
for d in [1, 3, 5, 9, 11, 13]:
    core = [d*i for i in range(1, 13)]
    for m in [13, 17, 19, 23, 29, 1]:
        rr = 14*m
        if rr in core: continue
        S = sorted(set(core + [rr]))
        if not is_Cprime_config(S): continue
        L, lo, hi = L_estimate(S, qs=(14000, 15013))
        results.append((L, f"d{d}m{m}", tuple(S)))
results.sort(key=lambda x: x[0])
print("lowest 14 L found (PRIMITIVE mult-of-14 configs):")
for L, tag, s in results[:14]:
    print(f"   r/tag={str(tag):8s}  L={L:.5f}   set={s[:4]}...{s[-1]}")
ap_min = results[0][0]
print(f"   --> min L over evader/dilated cores = {ap_min:.5f}")

# random PRIMITIVE small-entry configs containing a multiple of 14 (high additive energy)
print("\nrandom small-entry PRIMITIVE mult-of-14 configs (high additive energy):")
randmin = 1.0; randbest = None; tried = 0
while tried < 600:
    base = sorted(random.sample(range(1, 45), 12))
    mult = random.choice([14, 28, 42, 56, 70, 84])
    S = sorted(set(base + [mult]))
    if not is_Cprime_config(S):
        continue
    tried += 1
    L, lo, hi = L_estimate(S, qs=(14000, 15013))
    if L < randmin:
        randmin = L; randbest = S
print(f"   min L over {tried} random primitive configs = {randmin:.5f}   at {randbest}")
overall = min(ap_min, randmin)
print(f"\nOVERALL inf-L evidence (PRIMITIVE C' configs): min = {overall:.5f}  > 0 ?  {overall > 0}")

# ---- FRONT 2: the 7-vanishing s(7j)=0 ----
print("\n" + "=" * 74)
print("FRONT 2: the sinc 7-vanishing  s(t)=sin(pi t/7)/(pi t),  s(7j)=0")
print("=" * 74)
def s(t):
    return math.sin(math.pi * t / 7) / (math.pi * t) if t != 0 else 1/7
print("  s(t) for t=1..14:", [f"{s(t):+.4f}" for t in range(1, 15)])
print(f"  s(7) = {s(7):.2e}, s(14) = {s(14):.2e}  -> VANISH at multiples of 7")
print("  Consequence: any exact relation sum t_v v = 0 with a coefficient t_v in 7Z drops out.")
# empirical: a set that is all multiples of 7 (so 7 | every speed) -- its relations have 7|coeffs?
# Speeds 14*{1..13} are all multiples of 14=2*7. relation 1*v1 - ... ; coefficients t are integers,
# the relation lattice doesn't force t in 7Z, but the SPEEDS being 7-rich changes which q fire.
allmul7 = [14*7*i for i in range(1,14)]  # extra factor 7
L7, _, _ = L_estimate([x for x in allmul7], qs=(14000,15013))
print(f"  speeds 98*{{1..13}} (all div by 7^... ) L ~ {L7:.5f}  (vs tight {Lt:.4f})")

# ---- FRONT 3: adelic / single-prime localization test (HYP-2503) ----
print("\n" + "=" * 74)
print("FRONT 3: Euler-product / adelic test  L =? beta_inf * prod_p beta_p   (HYP-2503)")
print("=" * 74)
Stest = sorted(set([7*i for i in range(1,13)] + [13]))   # evader 7*{1..12} u {13}, primitive, contains 14
assert is_Cprime_config(Stest), Stest
print(f"test S = 7*{{1..12}} u {{13}} (primitive evader, high energy).  L (large q) ~ {L_estimate(Stest)[0]:.5f}")
print("single-prime sequences D(p^e,S)/p^e (does each converge to a NONTRIVIAL local factor?):")
for p in (2, 3, 5, 7, 11, 13):
    seq = []
    e = 1
    while p**e <= 60000:
        if p**e >= 14:
            seq.append(deficit(p**e, Stest) / p**e)
        e += 1
    print(f"   p={p:2d}: D(p^e)/p^e = {[f'{x:.4f}' for x in seq[-4:]]}  (last={seq[-1]:.4f} if any)" if seq else f"   p={p}: (no p^e in range)")
print("\nINTERPRETATION: if every single-prime limit ~ (6/7)^13 = main term (no genuine")
print("per-prime suppression), then L does NOT factor as an Euler product -- L is the")
print("ARCHIMEDEAN exact-relation sum, and the p-adic data lives in the APPROACH (threshold")
print("q*), not in L.  (Would REFINE/correct HYP-2503's Euler-product framing.)")

print("\nDONE.")
