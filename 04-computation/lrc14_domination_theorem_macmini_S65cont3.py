#!/usr/bin/env python3
"""
lrc14_domination_theorem_macmini_S65cont3.py -- THM-674 verification + the extended dodger

THE DOMINATION THEOREM (prime k > 14, danger radius j = ceil(k/14)-1, residues R nonzero):
  in G = Z_k^* / {+-1} (cyclic, order m = (k-1)/2), C = occupied classes, D = C^{-1},
  T = {classes of 1..j}:      BLOCKED  <=>  T*D = G.
  Proof: s bad <=> exists r in C, d in {+-1..+-j}: r*s = d  <=> class(s) in T*C^{-1}. QED
  (j=1: T={1}: blocked <=> D=G <=> ALL classes occupied = THM-672 T2.
   k=29: 2 primitive => log_2 iso, T={0,1}: blocked <=> indD u (indD+1) = Z/14
   <=> NO TWO CONSECUTIVE HOLES -- the cycle-domination statement.)
COROLLARY (j=2 primes): blocked => #classes >= sum of ceil(cycle_len/2) over the +t2-orbit
  cycles of Z/m, t2 = ind(2):  k=29: >=7/14; k=31 (ord2=5): >=9/15; k=37: >=9/18; k=41: >=10/20.
GENERAL LEDGER (all k, j >= 2): |A_r \\ 0| = (g-1) + 2g*floor(j/g), g = gcd(r,k)
  (non-units now reach danger elements divisible by g -- the [15,28] nesting breaks at j>=2).

Then: THE EXTENDED DODGER -- block ALL descents k in [15,42]; diagnose captures.
"""
from itertools import combinations
from math import gcd
from functools import reduce
import random

random.seed(74)

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))

def danger_j(k):
    return -(-k // 14) - 1          # ceil(k/14) - 1

def blocked(R, k):
    """Exact: no s in [1,k-1] with all r*s mod k in the closed band [ceil(k/14), k-ceil(k/14)]."""
    j = danger_j(k)
    bad = set(range(0, j + 1)) | set(range(k - j, k))
    for s in range(1, k):
        if all((r * s) % k not in bad for r in R):
            return False
    return True

def classes_of(R, k):
    return {min(r % k, k - r % k) for r in R if r % k}

def TD_covers(R, k):
    """Prime k: does T*D = G? T = classes of 1..j, D = inverse classes of occupied."""
    j = danger_j(k)
    C = classes_of(R, k)
    D = {pow(c, -1, k) for c in C}
    TD = set()
    for d in range(1, j + 1):
        for x in D:
            y = d * x % k
            TD.add(min(y, k - y))
    G = set(range(1, (k - 1) // 2 + 1))
    return TD >= G

# ---------------------------------------------------------------- verify the domination theorem
print("=" * 100)
print("VERIFY THM-674: blocked <=> T*D = G   (primes k = 29, 31, 37, 41 j=2; k = 43, 53 j=3)")
print("=" * 100)
for k in (29, 31, 37, 41, 43, 53):
    m = (k - 1) // 2
    viol = blk = 0
    min_classes = None
    pool = list(range(1, k))
    for _ in range(150000):
        R = tuple(sorted(random.sample(pool, 13)))
        b = blocked(R, k)
        t = TD_covers(R, k)
        if b != t:
            viol += 1
        if b:
            blk += 1
            nc = len(classes_of(R, k))
            min_classes = nc if min_classes is None else min(min_classes, nc)
    for _ in range(20000):                       # multisets too
        R = [random.randrange(1, k) for _ in range(13)]
        if blocked(R, k) != TD_covers(R, k):
            viol += 1
    print(f"k={k} (j={danger_j(k)}, m={m}): 170k samples, blocked found {blk}, "
          f"TD-characterization violations {viol}, min #classes among blocked: {min_classes}")

# ---------------------------------------------------------------- corollary: cycle bounds (j=2)
print()
print("=" * 100)
print("COROLLARY (j=2 primes): blocked => #classes >= per-cycle domination bound")
print("=" * 100)
def cycle_bound(k):
    """sum of ceil(len/2) over +ind(2)-orbit cycles of Z/m."""
    m = (k - 1) // 2
    # order of 2-bar in G = Z_k^*/{+-1}
    x, o = 2 % k, 1
    while min(x, k - x) != 1:
        x = x * 2 % k
        o += 1
    ncyc = gcd(o, m) and m // o          # cycles = m / ord(2bar), each length ord(2bar)
    return ncyc * ((o + 1) // 2), o
for k in (29, 31, 37, 41):
    bnd, o = cycle_bound(k)
    viol = blk = 0
    pool = list(range(1, k))
    minc = None
    for _ in range(200000):
        R = tuple(sorted(random.sample(pool, 13)))
        if blocked(R, k):
            blk += 1
            nc = len(classes_of(R, k))
            minc = nc if minc is None else min(minc, nc)
            if nc < bnd:
                viol += 1
    print(f"k={k}: ord(2bar)={o}, bound >= {bnd}: blocked {blk}, violations {viol}, "
          f"min observed {minc}")

# ---------------------------------------------------------------- general ledger formula
print()
print("=" * 100)
print("GENERAL LEDGER: |A_r \\ 0| = (g-1) + 2g*floor(j/g)  -- verify all k in [29,42], all g | k")
print("=" * 100)
ok = True
for k in range(29, 43):
    j = danger_j(k)
    bad = set(range(0, j + 1)) | set(range(k - j, k))
    for r in range(1, k):
        g = gcd(r, k)
        A = sum(1 for s in range(1, k) if (r * s) % k in bad)
        pred = (g - 1) + 2 * g * (j // g)
        if A != pred:
            ok = False
            print(f"  MISMATCH k={k} r={r} g={g}: |A\\0|={A} pred={pred}")
print(f"ledger formula exact for all k in [29,42], all r: {ok}")

# ---------------------------------------------------------------- THE EXTENDED DODGER [15,42]
print()
print("=" * 100)
print("EXTENDED DODGER: covering sets blocking ALL descents k in [15,42]; capture diagnosis")
print("=" * 100)
def covering(S):  return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1
def rulers_of(S): return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})

def open_descents(S, kmax=42):
    n = 0
    for q in rulers_of(S):
        for k in range(15, min(kmax, q - 1) + 1):
            if q % k == 0:
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    n += 1
    return n

def diagnose(S):
    """What catches this set: C0? C1? descent k<=42? descent k>42? or exact-max only?"""
    out = []
    q0 = min(S) + max(S)
    if max(S) <= 13 * min(S):
        out.append(f"C0 (q={q0})")
    for q in rulers_of(S):
        for k in range(15, q):
            if q % k == 0:
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    out.append(f"C2 k={k} (q={q}){' ABOVE-42' if k > 42 else ''}")
                    return out
    # exact liveness on pair-sum rulers
    for q in rulers_of(S):
        for p in range(1, q):
            if all(14 * (v * p % q) >= q and 14 * (q - v * p % q) >= q for v in S):
                out.append(f"exact-only live q={q} p={p}")
                return out
    out.append("NOT LONELY?!")
    return out

for cap in (150, 250):
    best = None
    for restart in range(6):
        while True:
            S = sorted(random.sample(range(1, cap + 1), 13))
            if covering(S) and primitive(S):
                break
        cur = open_descents(S)
        for step in range(250):
            T = list(S)
            T[random.randrange(13)] = random.randrange(1, cap + 1)
            T = sorted(set(T))
            if len(T) != 13 or not covering(T) or not primitive(T):
                continue
            sc = open_descents(T)
            if sc <= cur:
                S, cur = T, sc
            if cur == 0:
                break
        if best is None or cur < best[0]:
            best = (cur, S)
        if cur == 0:
            break
    print(f"cap {cap}: min OPEN [15,42]-descents = {best[0]}  S = {best[1]}")
    if best[0] == 0:
        print(f"   FULL [15,42] dodge -> capture: {diagnose(best[1])}")
    else:
        print(f"   (no full dodge found at this cap)")

# ---------------------------------------------------------------- spread-vs-concentration stats
print()
print("=" * 100)
print("TENSION STATS: covering sets at cap 250 -- occupied-class counts mod primes 29..41")
print("=" * 100)
from collections import Counter
cnt = {29: Counter(), 31: Counter(), 37: Counter(), 41: Counter()}
nblocked = Counter()
for _ in range(400):
    while True:
        S = sorted(random.sample(range(1, 251), 13))
        if covering(S) and primitive(S):
            break
    for k in cnt:
        nc = len(classes_of(S, k))
        cnt[k][nc] += 1
        if blocked([v % k for v in S if v % k] or [1], k):
            nblocked[k] += 1
for k in cnt:
    bnd, _ = cycle_bound(k)
    dist = sorted(cnt[k].items())
    lo = min(cnt[k]); hi = max(cnt[k])
    print(f"k={k}: occupied classes range [{lo},{hi}] (domination needs >= {bnd} + pattern); "
          f"blocked instances: {nblocked[k]}/400")
print()
print("Done.")
