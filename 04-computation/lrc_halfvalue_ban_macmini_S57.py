#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S57 -- HYP-4132: the HALF-VALUE BAN CRITERION at s = 2p.

Derivation (margin 2/25, s = 2p, p odd prime, mu = ceil(4p/25) = 3 for p in {13,17}):
  dilations a: odd units mod 2p <-> units mod p; margins symmetric under a -> -a.
  EVEN runner v = 2u: blocked <=> dist_p(u a) <= 1  (E_p = {0,+-2} mod 2p -- always).
  ODD  runner v:      blocked <=> v a == +-1 mod 2p (O_p = {+-1}, mu = 3).
  BAN of v = the +-class of rho(v)^{-1} mod p, where RHO(v) = v (v odd), v/2 (v even):
    odd-odd bans coincide  <=> v == +-w (mod p)        [odd residues mod 2p rigidity]
    even-even              <=> v/2 == +-w/2 (mod p)
    mixed                  <=> w == +-2v (mod p)       [the doubling fold]
  KILL-ALL: v == 0 mod 2p.  SAFE (no ban): v odd, p | v (antipode distance p).

  CRITERION: witness at 2p  <=>  no runner == 0 (mod 2p)  AND
             |{ +-rho(v) mod p : v in W, p not| rho(v) }| <= (p-1)/2 - 1.

Tests here (all small, exact):
  T1  the criterion vs brute-force witness search at s = 26, 34 on random sets
      + census-style survivors + the AP + lifts (exactness check).
  T2  the AP's rho-classes = perfect DOUBLE COVER at every p (the rigidity shape).
  T3  classify rho-FULL 12-sets at p = 13 (which families double-cover all 6):
      how much freedom beyond the AP?  (the coincidence-dividend question)
  T4  search for profile-passing families rho-full at ALL of p = 13, 17, 19, 23
      (mu(19) = mu(23) = 4: odd bans double -- criterion adjusted: odd runners
      ban TWO classes +-v^{-1}, +-3v^{-1}; evens still one) -- if none exist
      except AP-shapes, the 4-prime family alone discharges the template crux.
"""
from math import gcd
from functools import reduce
from itertools import combinations
import random, sys, time

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(57)

def dist(x, q):
    x %= q
    return min(x, q - x)

def brute_witness_2p(W, p):
    """exact: exists odd unit a mod 2p with all dist_{2p}(v a) >= mu."""
    s = 2 * p
    mu = -(-4 * p // 25)          # ceil(4p/25)
    for a in range(1, s, 2):
        if gcd(a, s) != 1:
            continue
        if all(dist(v * a, s) >= mu for v in W):
            return a
    return 0

def rho_classes(W, p):
    """the +-half-value classes mod p (None if some v == 0 mod 2p)."""
    cls = set()
    for v in W:
        if v % (2 * p) == 0:
            return None
        r = (v if v % 2 else v // 2) % p
        if r == 0:
            continue                       # odd p-multiple: safe, no ban
        cls.add(min(r, p - r))
    return cls

def criterion_mu3(W, p):
    cls = rho_classes(W, p)
    if cls is None:
        return False
    return len(cls) <= (p - 1) // 2 - 1

def criterion_mu4(W, p):
    """mu=4 (p = 19, 23): odd v bans +-v^{-1} AND +-3 v^{-1}; even v = 2u bans
    +-u^{-1} (E_p still {0,+-2}: dist_p(ua) <= 1).  Ban classes of v:
    odd: {+-v^{-1}, +-3v^{-1}} ~ classes of rho with rho = v and v/3;
    for the criterion count the union of banned +-dilation-classes."""
    s = 2 * p
    if any(v % s == 0 for v in W):
        return False
    banned = set()
    for v in W:
        if v % 2:
            if v % p == 0:
                continue
            vi = pow(v % p, -1, p)
            banned.add(min(vi % p, p - vi % p))
            b2 = (3 * vi) % p
            banned.add(min(b2, p - b2))
        else:
            u = (v // 2) % p
            if u == 0:
                return False               # even p-multiple: kill-all mod 2p? v=2u, p|u => 2p|v handled; p|u <=> 2p|v -- unreachable
            ui = pow(u, -1, p)
            banned.add(min(ui, p - ui))
    return len(banned) <= (p - 1) // 2 - 1

# ---------------- T1: exactness ----------------
log("T1: criterion vs brute force at s = 26 (p = 13) and s = 34 (p = 17)")
AP = list(range(1, 13))
tests = [AP]
for _ in range(4000):
    tests.append(sorted(random.sample(range(1, 200), 12)))
# lifts and doubled shapes
for _ in range(1000):
    C = random.sample(range(1, 13), random.randint(1, 4))
    W = [v for v in AP if v not in C] + [c + 13 * random.randint(1, 8) for c in C]
    tests.append(sorted(W))
bad = 0
for W in tests:
    for p in (13, 17):
        pred = criterion_mu3(W, p)
        act = brute_witness_2p(W, p) != 0
        if pred != act:
            bad += 1
            if bad <= 5:
                log(f"  MISMATCH p={p} W={W} pred={pred} act={act}")
log(f"  {len(tests)} sets x 2 primes: mismatches = {bad}")
log("T1b: mu=4 criterion at s = 38, 46 (p = 19, 23)")
bad4 = 0
for W in tests[:2000]:
    for p in (19, 23):
        pred = criterion_mu4(W, p)
        act = brute_witness_2p(W, p) != 0
        if pred != act:
            bad4 += 1
            if bad4 <= 5:
                log(f"  MISMATCH p={p} W={W} pred={pred} act={act}")
log(f"  2000 sets x 2 primes: mismatches = {bad4}")

# ---------------- T2: AP rho-structure ----------------
log("\nT2: the AP's rho multiset at p = 13, 17, 19, 23:")
for p in (13, 17, 19, 23):
    ms = {}
    for v in AP:
        r = (v if v % 2 else v // 2) % p
        c = min(r, p - r)
        ms[c] = ms.get(c, 0) + 1
    log(f"  p={p}: classes {dict(sorted(ms.items()))}  (need {(p-1)//2} classes for full)")

# ---------------- T3: rho-full 12-sets at p=13 ----------------
log("\nT3: rho-FULL at p=13 = perfect double cover (12 runners, 6 classes twice)?")
log("    structure: rho(v) in class c means v in {c, 2c, p-c, 2(p-c)} + 13Z x parity...")
# count rho-full 12-subsets of [1, 30] (small window, exact)
cnt_full = 0
cnt_all = 0
examples = []
for W in combinations(range(1, 31), 12):
    cnt_all += 1
    cls = rho_classes(W, 13)
    if cls is not None and len(cls) == 6:
        cnt_full += 1
        if len(examples) < 3 and sorted(W) != AP:
            examples.append(W)
log(f"  12-subsets of [1,30]: {cnt_all} total, rho-full = {cnt_full} ({100*cnt_full/cnt_all:.1f}%)")
log(f"  (rho-full = NO witness at 26; the census q=26 mass = the rho-deficient)")

# ---------------- T4: 4-prime simultaneous rho-fullness + profile ----------------
log("\nT4: profile-passing families failing the criterion at ALL of 26/34/38/46:")
def profile_pass(W):
    Ws = sorted(W)
    for m in range(2, 13):
        if not any(v % m == 0 for v in W):
            return False
    if 2 * Ws[-1] <= 23 * Ws[0]:
        return False
    if Ws[-1] > 24 * Ws[-2]:
        return False
    if 2 * Ws[-1] < 38:
        return False
    if reduce(gcd, W) != 1:
        return False
    for q in range(2, 26):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if not any((v * a) % q in (0, 1, q - 1) for v in W):
                return False
    return True
found = []
n_prof = 0
t0 = time.time()
# search census-style: random 12-subsets of [1,48] passing profile, test 4-prime failure
tries = 0
while time.time() - t0 < 120 and tries < 4_000_000:
    tries += 1
    W = sorted(random.sample(range(1, 49), 12))
    if max(W) < 25:
        continue
    if not profile_pass(W):
        continue
    n_prof += 1
    if (not criterion_mu3(W, 13) and not criterion_mu3(W, 17)
            and not criterion_mu4(W, 19) and not criterion_mu4(W, 23)):
        found.append(W)
        if len(found) <= 8:
            log(f"  4-PRIME-BLOCKED profile family: {W}")
log(f"  profile samples: {n_prof} (from {tries} tries); 4-prime-blocked: {len(found)}")
log(f"  (if such families exist, they must witness at some OTHER q <= 50 -- kps's (ii);")
log(f"   if none, the 2p-family {{26,34,38,46}} ALONE could carry the template dichotomy)")
log(f"\n[t = {time.time()-T0:.0f}s]")
