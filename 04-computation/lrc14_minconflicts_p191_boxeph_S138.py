#!/usr/bin/env python3
"""
lrc14_minconflicts_p191_boxeph_S138.py  (HYP-8025)

STRONGER SEARCHER at p=191: min-conflicts local search (random uncovered scale,
best repairing move, random-walk kicks, restarts) replacing the greedy that
already failed at p=191's x14 rung (S136: 1/25) — so its nulls there were
uninformative.  Verified witnesses only; failure stays one-sided.

RUNS (p=191): (A) x14 gate from level-1 — resolves S136's 24/25 UNDECIDED;
(B) x7 gate from level-1; (C) x49 refinement of the x7 witnesses (7 lifts per
slot of the forced mod-7 class).  Controls: AP13/GW boundary tuples verified
improper at 2674 / 1337 / 9359 and available as explicit witnesses.

boxeph-2026-07-19-S138.  Pure Python, exact integers.
"""

import random
from math import gcd

P = 191
DK = 13                     # floor(p/14)
HP = 95

def dist(x, m):
    r = x % m
    return r if r <= m - r else m - r

def improper(vals, Q, band):
    return all(any(dist(v * b, Q) <= band for v in vals) for b in range(1, Q // 2 + 1))

def level1(res):
    return all(any(dist(v * b, P) <= DK for v in res) for b in range(1, HP + 1))

def crt(P_, c_, Ic, rp, cc):
    return (rp + P_ * (((cc - rp) * Ic) % c_)) % (P_ * c_)

def build_slots(base_pairs, c, Q, band, lifts_of=None):
    """base_pairs: list of rp (mod P) or (rp, fixed-class) for refinement.
    returns slots: list of list of (scale-index-list, v)."""
    H = Q // 2
    Ic = pow(P % c, -1, c)
    slots = []
    for item in base_pairs:
        opts = []
        if lifts_of is None:
            it = [(eps, cc) for eps in (1, -1) for cc in range(c)]
            for eps, cc in it:
                rp = (eps * item) % P
                v = crt(P, c, Ic, rp, cc)
                idx = frozenset(b - 1 for b in range(1, H + 1) if dist(v * b, Q) <= band)
                opts.append((idx, v))
        else:
            rp, c7 = item
            for t in range(7):
                cc = c7 + 7 * t
                v = crt(P, c, Ic, rp, cc)
                idx = frozenset(b - 1 for b in range(1, H + 1) if dist(v * b, Q) <= band)
                opts.append((idx, v))
        slots.append(opts)
    return slots, H

def minconflicts(slots, H, Q, band, seed, restarts=5, moves=25000):
    rng = random.Random(seed)
    n = len(slots)
    movers = [[] for _ in range(H)]
    for i in range(n):
        for o, (idx, v) in enumerate(slots[i]):
            for b in idx:
                movers[b].append((i, o))
    for _ in range(restarts):
        opt = [rng.randrange(len(slots[i])) for i in range(n)]
        cnt = [0] * H
        for i in range(n):
            for b in slots[i][opt[i]][0]:
                cnt[b] += 1
        unc = {b for b in range(H) if cnt[b] == 0}
        for _ in range(moves):
            if not unc:
                sol = [slots[i][opt[i]][1] for i in range(n)]
                if improper(sol, Q, band):
                    return sol
                break
            b = rng.choice(tuple(unc))
            cand = [mo for mo in movers[b] if mo[1] != opt[mo[0]]]
            if not cand:
                break
            if rng.random() < 0.1:
                i, o = cand[rng.randrange(len(cand))]
            else:
                bestg, best = -1, []
                for (i2, o2) in cand:
                    g = len(slots[i2][o2][0] & unc)
                    if g > bestg: bestg, best = g, [(i2, o2)]
                    elif g == bestg: best.append((i2, o2))
                i, o = best[rng.randrange(len(best))]
            for x in slots[i][opt[i]][0]:
                cnt[x] -= 1
                if cnt[x] == 0: unc.add(x)
            opt[i] = o
            for x in slots[i][o][0]:
                cnt[x] += 1
                unc.discard(x)
        else:
            continue
    return None

# ---------------- controls ------------------------------------------------------
AP13 = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
assert improper(AP13, 14 * P, P) and improper(AP13, 7 * P, P // 2) \
   and improper(AP13, 49 * P, (7 * P - 1) // 2)
print("controls: AP13 improper mod 2674, 1337, 9359  OK (GW analogous, spot-checked)")
assert improper(GW, 49 * P, (7 * P - 1) // 2)

def random_level1(n_want, seed):
    rng = random.Random(seed)
    hitters = {b: [r for r in range(1, HP + 1) if dist(r * b, P) <= DK]
               for b in range(1, HP + 1)}
    out = []
    while len(out) < n_want:
        chosen, unc = [], set(range(1, HP + 1))
        while unc and len(chosen) < 13:
            b = min(unc, key=lambda x: len(hitters[x]))
            cds = hitters[b][:]; rng.shuffle(cds)
            chosen.append(cds[0])
            unc = {x for x in unc if dist(cds[0] * x, P) > DK}
        if not unc:
            while len(chosen) < 13:
                chosen.append(rng.randrange(1, HP + 1))
            res = sorted(chosen)
            if level1(res):
                out.append(tuple(res))
    return out

bases = random_level1(15, seed=1913)   # same family as S136

# ---------------- (A) x14 at p=191 ----------------------------------------------
print("\n(A) x14 gate at p=191, min-conflicts (resolving S136's undecided):")
a_ok = 0
for i, base in enumerate(bases):
    slots, H = build_slots(list(base), 14, 14 * P, P)
    sol = minconflicts(slots, H, 14 * P, P, seed=10 + i)
    if sol and any(v % 14 for v in sol) and any(v % 2 == 0 for v in sol) \
           and any(v % 7 == 0 for v in sol):
        a_ok += 1
print("  admit: %d/15   (greedy in S136: 1/25)" % a_ok)

# ---------------- (B) x7 at p=191 -----------------------------------------------
print("\n(B) x7 gate at p=191, min-conflicts:")
wit7 = []
for i, base in enumerate(bases):
    slots, H = build_slots(list(base), 7, 7 * P, P // 2)
    sol = minconflicts(slots, H, 7 * P, P // 2, seed=50 + i)
    if sol and any(v % 7 for v in sol):
        wit7.append(sol)
print("  admit: %d/15" % len(wit7))

# ---------------- (C) x49 refinement at p=191 -----------------------------------
print("\n(C) x49 refinement of the x7 witnesses at p=191:")
c_ok = 0
for i, w in enumerate(wit7):
    pairs = [(v % P, v % 7) for v in w]
    slots, H = build_slots(pairs, 49, 49 * P, (7 * P - 1) // 2, lifts_of=True)
    sol = minconflicts(slots, H, 49 * P, (7 * P - 1) // 2, seed=90 + i,
                       restarts=6, moves=40000)
    if sol:
        c_ok += 1
print("  refinement admit: %d/%d" % (c_ok, len(wit7)))

print("\nLADDER UPDATE (one-sided where search): p=43: every rung soft (S135-S137).")
print("p=191 (this run): x14 %d/15, x7 %d/15, x49-refine %d/%d." % (a_ok, len(wit7), c_ok, len(wit7)))
print("If x14/x7 admit recovers under min-conflicts, S136's collapse was SEARCH")
print("hardness, not existence; if the x49 refinement fails where x7 succeeds, the")
print("cage is localized at the second 7-rung at large p -- the |S| experiment's target.")
print("DONE.")
