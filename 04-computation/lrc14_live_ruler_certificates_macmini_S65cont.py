#!/usr/bin/env python3
"""
lrc14_live_ruler_certificates_macmini_S65cont.py -- HYP-5730: attacking the Schur-budget theorem

TARGET (the live-ruler theorem = LRC(14) covering branch, via THM-668):
  every primitive covering 13-set S has a pair-sum modulus q = v_i+v_j and multiplier p with all
  13 residues v_l*p mod q in the band [q/14, 13q/14].

THREE PROVED CERTIFICATES (elementary; proofs in THM-668 addendum):
  C1 LEDGER   : bad multipliers for runner l on ruler q form B_l = {p : v_l*p mod q in D},
                D = {d : |d| <= m}, m = ceil(q/14)-1.  |B_l| = g*(2*floor(m/g)+1), g = gcd(v_l,q).
                B_l = B_k iff v_l = +-v_k (mod q)  (D symmetric).  All B_l contain 0.
                => |bad ∩ [1,q-1]| <= Sum_classes (|B|-1).   LIVE if that sum < q-1.
  C2 DESCENT  : if k | q, k > 14, and exists s in [1,k-1] with all v_l*s mod k in
                [ceil(k/14), k-ceil(k/14)], then p = (q/k)*s is banded mod q.  (For 14 < k <= 28
                the band condition is exactly "avoid {0, +-1} mod k".)  Recursion on divisors.
  C3 SIX-PAIR : q prime in (Vmax, 2Vmax], q = t (mod 14) with t >= 3, and q has SIX pair-sum
                representations => >= t-1 live multipliers.  (7 merged classes, each unit-|D|;
                union bound gives exactly 14(ceil(q/14)-1) = q-t bad in [0,q-1] incl 0.)

This script: verify certificates against EXACT truth; census who certifies what; mine residuals.
All integer arithmetic.
"""
from itertools import combinations
from math import gcd
from functools import reduce
import random, sys

random.seed(66)

def covering(S):
    return all(any(v % k == 0 for v in S) for k in range(2, 15))

def primitive(S):
    return reduce(gcd, S) == 1

# ---------------------------------------------------------------- exact truth on a ruler
def live_exact(S, q):
    """Return a banded multiplier p in [1,q-1] or None.  Exact loop with early exit."""
    for p in range(1, q):
        ok = True
        for v in S:
            r = v * p % q
            if 14 * r < q or 14 * (q - r) < q:
                ok = False
                break
        if ok:
            return p
    return None

# ---------------------------------------------------------------- C1 ledger
def c1_ledger(S, q):
    """Return (fires: bool, ledger_sum, n_classes).  Requires undead (no residue 0)."""
    m = -(-q // 14) - 1                      # ceil(q/14) - 1
    classes = {}
    for v in S:
        r = v % q
        if r == 0:
            return (False, None, None)       # dead
        key = min(r, q - r)
        classes.setdefault(key, gcd(r, q))
    total = 0
    for key, g in classes.items():
        B = g * (2 * (m // g) + 1)
        total += B - 1                       # remove the shared 0
    return (total < q - 1, total, len(classes))

# ---------------------------------------------------------------- C2 descent
def c2_descent(S, q):
    """Smallest divisor k of q with 14 < k < q and band mod k solvable; None if none.
    (k = q allowed too but that's the exact check; we want PROPER descent k < q.)"""
    divs = sorted(k for k in range(15, q) if q % k == 0)
    for k in divs:
        for s in range(1, k):
            ok = True
            for v in S:
                r = v * s % k
                if 14 * r < k or 14 * (k - r) < k:
                    ok = False
                    break
            if ok:
                return k
    return None

# ---------------------------------------------------------------- C3 six-pair prime
def is_prime(n):
    if n < 2: return False
    for d in range(2, int(n**0.5) + 1):
        if n % d == 0: return False
    return True

def c3_sixpair(S, q):
    if q <= max(S) or not is_prime(q): return False
    t = q % 14
    if t < 3: return False
    reps = sum(1 for i in range(13) for j in range(i + 1, 13) if S[i] + S[j] == q)
    return reps == 6

# ---------------------------------------------------------------- per-set audit
def audit(S):
    S = sorted(S)
    rulers = sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})
    rows = []
    for q in rulers:
        if any(v % q == 0 for v in S):
            rows.append((q, 'DEAD', None, False, None, False))
            continue
        p = live_exact(S, q)
        fires1, ledger, ncls = c1_ledger(S, q)
        k2 = c2_descent(S, q)
        fires3 = c3_sixpair(S, q)
        rows.append((q, 'LIVE' if p else 'EMPTY', p, fires1, k2, fires3))
    return rows

def summarize(rows):
    live = [r for r in rows if r[1] == 'LIVE']
    c1 = [r for r in rows if r[3]]
    c2 = [r for r in rows if r[4] is not None]
    c3 = [r for r in rows if r[5]]
    cert = [r for r in rows if r[3] or r[4] is not None or r[5]]
    return live, c1, c2, c3, cert

# ---------------------------------------------------------------- soundness spot check
print("=" * 100)
print("SOUNDNESS: certificates never fire on a non-live ruler (500 random sets, all rulers)")
print("=" * 100)
bad = 0
for _ in range(500):
    S = sorted(random.sample(range(1, 61), 13))
    for (q, st, p, f1, k2, f3) in audit(S):
        if (f1 or k2 or f3) and st != 'LIVE':
            bad += 1
            print(f"  UNSOUND: S={S} q={q} status={st} C1={f1} C2={k2} C3={f3}")
print(f"unsound firings: {bad} / 500 sets x all rulers")

# ---------------------------------------------------------------- census [1,18]
print()
print("=" * 100)
print("CENSUS: 966 primitive covering 13-subsets of [1,18] -- do the certificates cover?")
print("=" * 100)
n = nC1 = nC2 = nC3 = nCert = 0
residuals = []
for S in combinations(range(1, 19), 13):
    if not (covering(S) and primitive(S)):
        continue
    n += 1
    rows = audit(list(S))
    live, c1, c2, c3, cert = summarize(rows)
    assert live, f"LRC violation?! {S}"
    nC1 += bool(c1); nC2 += bool(c2); nC3 += bool(c3); nCert += bool(cert)
    if not cert:
        residuals.append((list(S), rows))
print(f"{n} sets: C1 fires somewhere on {nC1} ({nC1/n:.1%}); C2 on {nC2} ({nC2/n:.1%}); "
      f"C3 on {nC3} ({nC3/n:.1%}); ANY on {nCert} ({nCert/n:.1%})")
print(f"RESIDUAL sets (live but uncertified): {len(residuals)}")
for S, rows in residuals[:6]:
    lv = [(q, p) for (q, st, p, *_ ) in rows if st == 'LIVE']
    print(f"  S={S}")
    for (q, st, p, f1, k2, f3) in rows:
        if st == 'LIVE':
            m = -(-q // 14) - 1
            _, ledger, ncls = c1_ledger(S, q)
            print(f"    live q={q} p={p}  ledger={ledger} vs q-1={q-1}  classes={ncls}  "
                  f"gcds={sorted(set(gcd(v, q) for v in S), reverse=True)[:5]}")

# ---------------------------------------------------------------- census cap 60 random
print()
print("=" * 100)
print("RANDOM census cap 60 (600 covering sets): certificate coverage at larger Vmax")
print("=" * 100)
n = nC1 = nC2 = nC3 = nCert = 0
res60 = []
while n < 600:
    S = sorted(random.sample(range(1, 61), 13))
    if not (covering(S) and primitive(S)):
        continue
    n += 1
    rows = audit(S)
    live, c1, c2, c3, cert = summarize(rows)
    if not live:
        print(f"  !!! M<1/14 candidate {S}")
        continue
    nC1 += bool(c1); nC2 += bool(c2); nC3 += bool(c3); nCert += bool(cert)
    if not cert:
        res60.append(S)
print(f"{n} sets: C1 {nC1/n:.1%}; C2 {nC2/n:.1%}; C3 {nC3/n:.1%}; ANY {nCert/n:.1%}; "
      f"residual {len(res60)}")

# ---------------------------------------------------------------- adversarial: defeat all certificates
print()
print("=" * 100)
print("ADVERSARIAL: hill-climb MINIMIZING #certified rulers (cap 60, 8 restarts x 200 steps)")
print("=" * 100)
def ncert_of(S):
    rows = audit(S)
    return sum(1 for r in rows if r[3] or r[4] is not None or r[5])

best = None
for restart in range(8):
    while True:
        S = sorted(random.sample(range(1, 61), 13))
        if covering(S) and primitive(S):
            break
    cur = ncert_of(S)
    for step in range(200):
        T = list(S)
        T[random.randrange(13)] = random.randrange(1, 61)
        T = sorted(set(T))
        if len(T) != 13 or not covering(T) or not primitive(T):
            continue
        nc = ncert_of(T)
        if nc <= cur:
            S, cur = T, nc
    print(f"  restart {restart}: min #certified rulers = {cur}  S={S}")
    if best is None or cur < best[0]:
        best = (cur, S)
print(f"ADVERSARIAL MIN #certified = {best[0]} at S = {best[1]}")
if best[0] == 0:
    rows = audit(best[1])
    print("  -> DEFEATS all certificates; its live rulers:")
    for (q, st, p, f1, k2, f3) in rows:
        if st == 'LIVE':
            _, ledger, ncls = c1_ledger(best[1], q)
            print(f"    q={q} p={p} ledger={ledger}/{q-1} classes={ncls}")
print()
print("Done.")
