#!/usr/bin/env python3
"""
lrc14_p191_brackets_conjunction_boxeph_S136.py  (HYP-7990)

Owner-directed: (1) the p=191 brackets, (2) the 43 AND 61 conjunction test.

SCALE-FREENESS CORRECTION (to my own S135 letter, point (i)): the x14 gate's raw
counting does NOT tighten with p — scales = 7p, per-option coverage ~ p, so the
coverage ratio is 13/7 and the per-scale hit rate ~1/7, both p-independent.  If
the p=191 brackets match p=43 (LB=0, small UB), the gate softness is STRUCTURAL
and all selectivity lives in the cross-prime CONJUNCTION — whose coupling runs
through the SHARED mod-14 parts (13 classes serving every prime's slaving law
at once).  Part 2 takes the first datum on that coupling.

PART 1 (p=191, Q=2674): zero-forcing LB + randomized-greedy UB (verified
witnesses) over random level-1-improper bases; controls AP13/GW by explicit
witness (z=0).
PART 2 (conjunction 43 AND 61): take improper-mod-602 witnesses (base + greedy
lift from the S135 machinery), FIX their mod-14 parts, free the mod-61 parts,
and measure (a) MC random-61-assignment survival mod 854, (b) greedy admit rate
(one-sided: greedy failure is NOT a refutation), (c) the level-1-mod-61
requirement on the 61-parts (pure-61 scales) as the visible sub-law.

boxeph-2026-07-19-S136.  Pure Python, exact integers.
"""

import random
from math import gcd

def dist(x, m):
    r = x % m
    return r if r <= m - r else m - r

def improper_mod(vals, Q, band):
    return all(any(dist(v * b, Q) <= band for v in vals) for b in range(1, Q // 2 + 1))

def make_ctx(P):
    Q, DK, HALF, HP = 14 * P, P // 14, 7 * P, (P - 1) // 2
    INV = pow(P % 14, -1, 14)
    return P, Q, DK, HALF, HP, INV

def crt14p(P, INV, rp, c14):
    return (rp + P * (((c14 - rp) * INV) % 14)) % (14 * P)

def slot_options(P, INV, r):
    Q, HALF = 14 * P, 7 * P
    opts = []
    for eps in (1, -1):
        rp = (eps * r) % P
        for c14 in range(14):
            v = crt14p(P, INV, rp, c14)
            m = 0
            for b in range(1, HALF + 1):
                if dist(v * b, Q) <= P:
                    m |= 1 << (b - 1)
            opts.append((m, v, v % 14 == 0))
    return opts

def level1_improper(res, P, DK):
    return all(any(dist(v * b, P) <= DK for v in res) for b in range(1, (P - 1) // 2 + 1))

def random_level1(P, DK, n_want, seed):
    rng = random.Random(seed)
    HP = (P - 1) // 2
    hitters = {b: [r for r in range(1, HP + 1) if dist(r * b, P) <= DK]
               for b in range(1, HP + 1)}
    out = []
    while len(out) < n_want:
        chosen, unc = [], set(range(1, HP + 1))
        while unc and len(chosen) < 13:
            b = min(unc, key=lambda x: len(hitters[x]))
            cands = hitters[b][:]
            rng.shuffle(cands)
            chosen.append(cands[0])
            unc = {x for x in unc if dist(cands[0] * x, P) > DK}
        if not unc:
            while len(chosen) < 13:
                chosen.append(rng.randrange(1, HP + 1))
            res = sorted(chosen)
            if level1_improper(res, P, DK):
                out.append(tuple(res))
    return out

def brackets(P, base, restarts, seed):
    _, Q, DK, HALF, HP, INV = make_ctx(P)
    FULL = (1 << HALF) - 1
    slots = [slot_options(P, INV, r) for r in base]
    nz = 0
    for so in slots:
        for (m, v, isz) in so:
            if not isz: nz |= m
    zero_only = FULL & ~nz
    LB = 1 if zero_only else 0
    rng = random.Random(seed)
    bestz = None
    for _ in range(restarts):
        covered, used, free, zeros, ok = 0, [None] * 13, list(range(13)), 0, True
        while covered != FULL:
            bestsc, cands = -10**9, []
            for i in free:
                for o in range(28):
                    m, v, isz = slots[i][o]
                    g = bin(m & ~covered).count('1')
                    sc = g - (3 if isz else 0)
                    if sc > bestsc: bestsc, cands = sc, [(i, o)]
                    elif sc == bestsc: cands.append((i, o))
            if not free or bestsc <= -10**9: ok = False; break
            i, o = cands[rng.randrange(len(cands))]
            m, v, isz = slots[i][o]
            if not (m & ~covered): ok = False; break
            used[i] = o; free.remove(i); covered |= m; zeros += isz
        if ok:
            sol = [slots[i][used[i]][1] if used[i] is not None
                   else crt14p(P, INV, base[i] % P, 1) for i in range(13)]
            if any(v % 2 == 0 for v in sol) and any(v % 7 == 0 for v in sol) \
               and any(v % 14 for v in sol) and improper_mod(sol, Q, P):
                if bestz is None or zeros < bestz:
                    bestz = zeros
    return LB, bestz, bin(zero_only).count('1')

# =============================== PART 1: p=191 =================================
print("=" * 92)
print("PART 1: p=191 brackets (Q=2674, dk=13; scale-free ratio 13/7 -- see header)")
print("=" * 92)
P1, Q1, DK1 = 191, 2674, 13
AP13 = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
assert level1_improper([v % P1 for v in AP13], P1, DK1)
assert improper_mod(AP13, Q1, P1) and improper_mod(GW, Q1, P1)
print("controls: AP13/GW improper mod 2674, z=0 explicit => min-z = 0 exactly")
bases191 = random_level1(P1, DK1, 25, seed=1913)
lbh, ubh, zoc = {}, {}, []
for i, base in enumerate(bases191):
    LB, UB, zo = brackets(P1, base, restarts=350, seed=100 + i)
    lbh[LB] = lbh.get(LB, 0) + 1
    ubh[UB] = ubh.get(UB, 0) + 1
    zoc.append(zo)
print("25 random level-1 bases at p=191:")
print("  zero-forcing LB histogram: %s   (zero-only scales: max %d)" % (dict(sorted(lbh.items())), max(zoc)))
print("  greedy UB histogram:       %s" % dict(sorted(ubh.items(), key=lambda x: (x[0] is None, x[0]))))

# =============================== PART 2: 43 AND 61 =============================
print("\n" + "=" * 92)
print("PART 2: the 43 AND 61 conjunction (shared mod-14 parts; mod-61 parts free)")
print("=" * 92)
P43, P61 = 43, 61
Q43, Q61 = 602, 854
DK43, DK61 = 3, 4
INV43 = pow(P43 % 14, -1, 14)
# improper-mod-602 witnesses from the S135 machinery
bases43 = random_level1(P43, DK43, 20, seed=1443)
witnesses = []
for i, base in enumerate(bases43):
    _, Q, DK, HALF, HP, INV = make_ctx(P43)
    FULL = (1 << HALF) - 1
    slots = [slot_options(P43, INV, r) for r in base]
    rng = random.Random(500 + i)
    for _ in range(400):
        covered, used, free, ok = 0, [None] * 13, list(range(13)), True
        while covered != FULL:
            bestsc, cands = -10**9, []
            for ii in free:
                for o in range(28):
                    m, v, isz = slots[ii][o]
                    sc = bin(m & ~covered).count('1') - (3 if isz else 0)
                    if sc > bestsc: bestsc, cands = sc, [(ii, o)]
                    elif sc == bestsc: cands.append((ii, o))
            if not free: ok = False; break
            ii, o = cands[rng.randrange(len(cands))]
            m, v, isz = slots[ii][o]
            if not (m & ~covered): ok = False; break
            used[ii] = o; free.remove(ii); covered |= m
        if ok:
            sol = [slots[ii][used[ii]][1] if used[ii] is not None
                   else crt14p(P43, INV, base[ii] % P43, 1) for ii in range(13)]
            if improper_mod(sol, Q43, P43) and any(v % 14 for v in sol):
                witnesses.append(sol); break
print("improper-mod-602 witnesses built: %d/20" % len(witnesses))

INV61 = pow(P61 % 14, -1, 14)
HALF61 = 427
FULL61 = (1 << HALF61) - 1
mc_ok = mc_try = 0
admit = 0
lvl1_note = 0
rng = random.Random(4361)
for w in witnesses:
    c14s = [v % 14 for v in w]
    # option masks: for each slot, c61 in Z/61 -> v61 = CRT(c14 fixed, c61) mod 854
    slots61 = []
    for c14 in c14s:
        opts = []
        for c61 in range(61):
            v = (c61 * 14 * pow(14, -1, 61) + c14 * 61 * pow(61, -1, 14)) % Q61
            m = 0
            for b in range(1, HALF61 + 1):
                if dist(v * b, Q61) <= P61:
                    m |= 1 << (b - 1)
            opts.append((m, v))
        slots61.append(opts)
    # (a) MC random assignment
    for _ in range(400):
        vals = [slots61[i][rng.randrange(61)][1] for i in range(13)]
        mc_try += 1
        if improper_mod(vals, Q61, P61):
            mc_ok += 1
    # (b) greedy admit (one-sided)
    got = False
    for t in range(250):
        covered, used, free, ok = 0, [None] * 13, list(range(13)), True
        while covered != FULL61:
            bestsc, cands = -10**9, []
            for ii in free:
                for o in range(61):
                    m, v = slots61[ii][o]
                    sc = bin(m & ~covered).count('1')
                    if sc > bestsc: bestsc, cands = sc, [(ii, o)]
                    elif sc == bestsc: cands.append((ii, o))
            if not free: ok = False; break
            ii, o = cands[rng.randrange(len(cands))]
            m, v = slots61[ii][o]
            if not (m & ~covered): ok = False; break
            used[ii] = o; free.remove(ii); covered |= m
        if ok:
            vals = [slots61[ii][used[ii]][1] if used[ii] is not None
                    else slots61[ii][1][1] for ii in range(13)]
            if improper_mod(vals, Q61, P61):
                got = True
                if level1_improper([v % P61 for v in vals], P61, DK61):
                    lvl1_note += 1
                break
    admit += got
print("(a) MC random-61-assignment survival mod 854: %d/%d" % (mc_ok, mc_try))
print("(b) greedy conjunction admit (one-sided >=): %d/%d witnesses" % (admit, len(witnesses)))
print("    (each admit's 61-parts are level-1-improper mod 61: %d/%d -- the pure-61 sub-law)" % (lvl1_note, admit))
print("\nREADING: if p=191 brackets match p=43 (LB=0, small UB), gate softness is")
print("structural (scale-free 13/7); the conjunction admit rate is then the first")
print("measured rung of the cage -- the shared mod-14 parts must serve BOTH primes'")
print("slaving laws, and each added prime multiplies the constraint on them.")
print("DONE.")
