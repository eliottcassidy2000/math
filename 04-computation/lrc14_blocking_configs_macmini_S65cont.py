#!/usr/bin/env python3
"""
lrc14_blocking_configs_macmini_S65cont.py -- HYP-5730 phase 3: the BLOCKING-CONFIGURATION census

Small pair-sum rulers (q <= ~36) defeat counting certificates (union bound too weak); the exact
band check there is a BOUNDED, ABSOLUTE computation.  The theorem-shaped question becomes:

    Which 13-element subsets R of (Z/q) \\ {0} are BAND-BLOCKED
    (no multiplier p in [1,q-1] with all r*p mod q in [ceil(q/14), q - ceil(q/14)])?

If blocked configs are rare and rigid (AP/dilate-like), the live-ruler theorem factors:
[large rulers: counting certificates C0/C1/C2/C3] + [small rulers: blocked-config rigidity] +
[covering sets cannot be blocked-rigid on every small ruler AND evade counting on every large one].

Phase 3 additions:
  C0 (window certificate): all speeds in [q/14, 13q/14] for some pair sum q => p = 1 live.
      (= spread13_lonely realized on a pair-sum ruler; catches r <= 13 sets.)
  BLOCKING CENSUS: q = 15..26 exhaustive over 13-subsets of Z/q \\ {0}, up to UNIT DILATION;
      report count + structure (is every blocked config a dilate of a near-interval?).
  SLIVER ADVERSARY: r in (13, 17), Vmin >= 18 (all pair sums > 36, small-ruler exact vacuous),
      vs C0 v C1 v C2.
All exact integer arithmetic.
"""
from itertools import combinations
from math import gcd
from functools import reduce
import random

random.seed(68)

# ---------------------------------------------------------------- band machinery
def banded_p(R, q):
    """Exact: multiplier p in [1,q-1] with all r*p mod q in closed band, else None."""
    for p in range(1, q):
        ok = True
        for r in R:
            x = r * p % q
            if 14 * x < q or 14 * (q - x) < q:
                ok = False
                break
        if ok:
            return p
    return None

# ---------------------------------------------------------------- blocking census
print("=" * 100)
print("BLOCKING CENSUS: 13-subsets of (Z/q)\\{0}, band [ceil(q/14), q-ceil(q/14)], up to dilation")
print("=" * 100)
print("q : #subsets  #blocked  #blocked-dilation-classes  structure of blocked classes")

def dilation_class_rep(R, q):
    """Canonical representative of {u*R : u unit mod q} (lexicographically smallest sorted tuple)."""
    best = None
    for u in range(1, q):
        if gcd(u, q) != 1:
            continue
        T = tuple(sorted(u * r % q for r in R))
        if best is None or T < best:
            best = T
    return best

def longest_ap_mod(R, q):
    """Longest AP (common difference d, any) contained in R viewed in Z/q (as a SET)."""
    Rs = set(R)
    best = 1
    for d in range(1, q):
        for a in Rs:
            L = 1
            x = (a + d) % q
            while x in Rs and L < len(R):
                L += 1
                x = (x + d) % q
            best = max(best, L)
    return best

for q in range(15, 27):
    nsub = nblk = 0
    reps = {}
    for R in combinations(range(1, q), 13):
        nsub += 1
        if banded_p(R, q) is None:
            nblk += 1
            rep = dilation_class_rep(R, q)
            reps.setdefault(rep, 0)
            reps[rep] += 1
    # structure of blocked classes
    desc = []
    for rep in sorted(reps)[:4]:
        lap = longest_ap_mod(rep, q)
        holes = (q - 1) - 13
        desc.append(f"{list(rep)} LAP={lap}")
    print(f"q={q}: {nsub:8d} {nblk:8d} {len(reps):6d}   {'; '.join(desc) if desc else 'NONE BLOCKED'}")

# ---------------------------------------------------------------- certificates for the sliver test
def covering(S):
    return all(any(v % k == 0 for v in S) for k in range(2, 15))

def primitive(S):
    return reduce(gcd, S) == 1

def rulers_of(S):
    return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})

def c0_window(S):
    """Some pair sum q with all speeds in [q/14, 13q/14]: p=1 live."""
    for q in rulers_of(S):
        if all(14 * v >= q and 14 * (q - v) >= q for v in S):
            return q
    return None

def c1_fires(S, q):
    m = -(-q // 14) - 1
    classes = {}
    for v in S:
        r = v % q
        if r == 0:
            return False
        classes.setdefault(min(r, q - r), gcd(r, q))
    return sum(g * (2 * (m // g) + 1) - 1 for g in classes.values()) < q - 1

def c2_fires(S, q):
    if any(v % q == 0 for v in S):
        return False
    for k in sorted(kk for kk in range(15, q) if q % kk == 0):
        for s in range(1, k):
            if all(14 * (v * s % k) >= k and 14 * (k - v * s % k) >= k for v in S):
                return True
    return False

def certified(S):
    if c0_window(S) is not None:
        return True
    return any(c1_fires(S, q) or c2_fires(S, q) for q in rulers_of(S))

def live_somewhere(S):
    return any(banded_p(S, q) is not None for q in rulers_of(S))

# ---------------------------------------------------------------- sliver adversary
print()
print("=" * 100)
print("SLIVER ADVERSARY: r = Vmax/Vmin in (13, 17), Vmin >= 18 (all pair sums > 36) vs C0|C1|C2")
print("=" * 100)
def random_sliver(vmin_lo=18, vmin_hi=23):
    while True:
        vmin = random.randrange(vmin_lo, vmin_hi + 1)
        vmax = random.randrange(13 * vmin + 1, 17 * vmin)
        mids = random.sample(range(vmin + 1, vmax), 11)
        S = sorted([vmin] + mids + [vmax])
        if len(set(S)) == 13 and covering(S) and primitive(S):
            return S

found_defeater = None
n_cert = 0
N_SLIVER = 250
for i in range(N_SLIVER):
    S = random_sliver()
    if certified(S):
        n_cert += 1
    else:
        if live_somewhere(S):
            found_defeater = S
            print(f"  RESIDUAL (live, uncertified): {S}  r={max(S)/min(S):.2f}")
            if i > 30:
                break
        else:
            print(f"  !!! M<1/14 candidate: {S}")
print(f"random sliver sets: {n_cert}/{N_SLIVER} certified by C0|C1|C2")

# hill-climb inside the sliver
best = None
for restart in range(5):
    S = random_sliver()
    def score(T):
        c0 = c0_window(T) is not None
        n = sum(1 for q in rulers_of(T) if c1_fires(T, q) or c2_fires(T, q))
        return (1000 if c0 else 0) + n
    cur = score(S)
    for step in range(200):
        T = list(S)
        idx = random.randrange(13)
        T[idx] = max(1, T[idx] + random.choice([-5, -3, -2, -1, 1, 2, 3, 5]))
        T = sorted(set(T))
        if len(T) != 13 or not covering(T) or not primitive(T):
            continue
        if min(T) < 18 or max(T) <= 13 * min(T):     # stay in the sliver
            continue
        sc = score(T)
        if sc <= cur:
            S, cur = T, sc
        if cur == 0:
            break
    print(f"  sliver hill-climb restart {restart}: min score = {cur}  S={S}")
    if best is None or cur < best[0]:
        best = (cur, S)
print(f"SLIVER adversarial min (C0 blocked + #C1/C2-certified rulers) = {best[0]}")
if best[0] == 0:
    S = best[1]
    lv = [(q, banded_p(S, q)) for q in rulers_of(S) if banded_p(S, q)]
    print(f"  DEFEATER of C0|C1|C2 in the sliver: {S}")
    print(f"  its live rulers (exact): {lv[:8]}")
print()
print("Done.")
