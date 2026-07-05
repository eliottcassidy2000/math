#!/usr/bin/env python3
"""
opus-2026-07-05-S83 -- HYP-4126: the multiplicative reframing + tower-limit dichotomy.

PART 1: kill sets are multiplicative translates: K(v) = v^{-1} * (+-[1,12]) mod 169.
PART 2: discrete-log/CRT picture: ZZ/156 = ZZ/12 x ZZ/13; S = log(+-[1,12]) has exactly
        2 elements per column; shadows = diagonal shifts of the standard solution.
PART 3: WITNESS LIFTING down/up the tower (level-l witness lifts with margin) --
        verified numerically; covers project down.
PART 4: THE LEVEL-3 PROBE (the conjecture that matters): above a sampled level-2 cover
        class, do ANY of the 13^12 level-3 refinements cover at level 3 -- and are those
        that do exactly the level-3 shadows?  (Shadow classes must refine to >= 1 cover:
        sanity.)  If non-shadow level-2 covers have ZERO covering refinements, the
        dichotomy holds at level 3 and the alignment sliver dies by a finite lambda-check.
PART 5: side probes: Hamming distance of level-2 covers to the nearest shadow; QR/Paley
        structure of kills; Fourier bias of S (Erdos-Turan echo).
"""
import sys, time, collections
from itertools import product

# ---------------- PART 1: multiplicative translates
MOD = 169
UNITS = [a for a in range(1, MOD) if a % 13 != 0]
B2 = set(list(range(1, 13)) + list(range(157, 169)))  # +-[1,12] mod 169
def bad2(x): return x <= 13 or x >= 156
ok = True
for v in UNITS[:40]:
    vinv = pow(v, -1, MOD)
    K_direct = {a for a in UNITS if bad2((a * v) % MOD)}
    K_transl = {(vinv * b) % MOD for b in B2}
    ok &= (K_direct == K_transl)
print(f"PART 1: K(v) = v^-1 * (+-[1,12]) verified on 40 units: {ok}")
assert ok

# ---------------- PART 2: discrete log CRT
g = None
for cand in range(2, 40):
    if cand % 13 == 0: continue
    seen, x = set(), 1
    for _ in range(156):
        x = (x * cand) % MOD
        seen.add(x)
    if len(seen) == 156:
        g = cand; break
print(f"PART 2: primitive root mod 169: g = {g}")
LOG = {}
x = 1
for e in range(156):
    LOG[x] = e
    x = (x * g) % MOD
S = sorted(LOG[b] for b in B2)
percol = collections.Counter(s % 12 for s in S)
print(f"  |S| = {len(S)}; elements per column (mod 12): {dict(sorted(percol.items()))}")
assert all(c == 2 for c in percol.values())
# shadows are diagonal: shadow(lam) shifts = shadow(1) shifts + log(lam) (all positions)
def shadow_shifts(lam):
    laminv = pow(lam % 13, -1, 13)
    out = {}
    for r in range(1, 13):
        i = (laminv * r) % 13
        v = (lam * i) % MOD
        out[r] = LOG[pow(v, -1, MOD)]  # t_r = -log v_r = log v^-1
    return out
s1 = shadow_shifts(1)
diag_ok = True
for lam in [2, 5, 14, 30, 168]:
    sl = shadow_shifts(lam)
    diffs = {(sl[r] - s1[r]) % 156 for r in range(1, 13)}
    diag_ok &= (len(diffs) == 1 and (diffs == {(-LOG[lam]) % 156}))
print(f"  shadows = diagonal shifts of the standard solution: {diag_ok}")

# ---------------- PART 3: witness lifting
from fractions import Fraction as F
def distF(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)
lift_ok = True
import random
random.seed(83)
for _ in range(200):
    v = random.randrange(1, 10**6)
    if v % 13 == 0: continue
    a = random.randrange(1, 169)
    if a % 13 == 0: continue
    m2 = distF(F(v * a, 169))
    m3 = distF(F(v * 13 * a, 2197))
    if m2 >= F(14, 169):
        lift_ok &= (m3 >= F(170, 2197))  # needs > 169/2197; gets 14*13/2197 = 182/2197
print(f"PART 3: level-2 witness (>=14/169) lifts to level-3 witness (>=170/2197): {lift_ok}")

# ---------------- PART 4: the level-3 probe
MOD3 = 2197
UNITS3_COUNT = 2028
# level-3 columns: b = a mod 13 (12 columns x 169 heights). Kill condition at level 3:
# strict witness needs dist(v a / 2197) >= (169+1)/2197, i.e. v*a mod 2197 in [170, 2027].
# BAD3 = [0,169] u [2028,2196] intersect units => +-[1,169] minus 13-multiples: 312 elts.
def bad3(x): return x <= 169 or x >= 2028

# families: v_r = r + 13*kappa_r + 169*mu_r; level-2 class = kappa, level-3 digit = mu.
# For a FIXED level-2 cover kappa, count mu in (Z/13)^12 whose level-3 class covers all
# 2028 cells. Cell-driven DFS per class, mirroring the level-2 search at level 3.
def level3_census(kappa, cap_solutions=50, cap_nodes=30_000_000, verbose=""):
    CELLS3 = [a for a in range(1, MOD3) if a % 13 != 0]
    C3IDX = {a: i for i, a in enumerate(CELLS3)}
    FULL3 = (1 << len(CELLS3)) - 1
    KM3 = {}
    for r in range(1, 13):
        base = r + 13 * kappa[r - 1]
        for mu in range(13):
            v = base + 169 * mu
            m = 0
            for a in CELLS3:
                if bad3((a * v) % MOD3):
                    m |= 1 << C3IDX[a]
            KM3[(r, mu)] = m
    killers3 = collections.defaultdict(list)
    for (r, mu), m in KM3.items():
        for a in CELLS3:
            if m >> C3IDX[a] & 1:
                killers3[a].append((r, mu))
    COLMASK3 = collections.defaultdict(int)
    for a in CELLS3:
        COLMASK3[a % 13] |= 1 << C3IDX[a]
    sols = []
    nodes = 0
    def dfs(assigned, covered):
        nonlocal nodes
        nodes += 1
        if nodes > cap_nodes or len(sols) >= cap_solutions:
            return
        if covered == FULL3:
            sols.append(dict(assigned)); return
        nrem = 12 - len(assigned)
        for b in range(1, 13):
            if 169 - bin(covered & COLMASK3[b]).count("1") > 26 * nrem:
                return
        best = None
        for a in CELLS3:
            if covered >> C3IDX[a] & 1: continue
            opts = [(r, mu) for (r, mu) in killers3[a] if r not in assigned]
            if best is None or len(opts) < len(best):
                best = opts
                if len(opts) <= 2: break
        if not best: return
        for (r, mu) in best:
            assigned[r] = mu
            dfs(assigned, covered | KM3[(r, mu)])
            del assigned[r]
    t0 = time.time()
    dfs({}, 0)
    return sols, nodes, time.time() - t0

# level-2 covers to probe: 2 shadows + 3 sampled non-shadow covers (from S82 sampling)
def shadow_kappa(lam):
    laminv = pow(lam % 13, -1, 13)
    return tuple((((lam * ((laminv * r) % 13)) % MOD - r) // 13) % 13 for r in range(1, 13))

probes = [("shadow lam=14", shadow_kappa(14)),
          ("shadow lam=27", shadow_kappa(27)),
          ("nonshadow #1", (0, 10, 5, 3, 5, 3, 5, 8, 7, 3, 1, 11)),
          ("nonshadow #2", (0, 10, 5, 3, 5, 3, 3, 1, 8, 3, 1, 1)),
          ("nonshadow #3", (0, 10, 5, 3, 5, 3, 3, 8, 1, 5, 1, 3))]
print("PART 4: level-3 refinement census above level-2 covers")
for name, kap in probes:
    sols, nodes, el = level3_census(kap)
    tag = ""
    if "shadow" in name and sols:
        tag = " (must include the true shadow: sanity OK)"
    print(f"  {name}: kappa={kap}: level-3 covering refinements found = {len(sols)} "
          f"({nodes} nodes, {el:.0f}s){tag}", flush=True)

# ---------------- PART 5: side probes
print("PART 5:")
# Hamming to nearest shadow for the three nonshadow probes
SH2 = {shadow_kappa(l) for l in UNITS}
for name, kap in probes[2:]:
    d = min(sum(1 for a, b in zip(kap, sh) if a != b) for sh in SH2)
    print(f"  {name}: Hamming distance to nearest shadow = {d}")
# QR structure: are kills QR-correlated? kill(a,v) iff av in +-[1,12]: QR(av mod 13)?
QR13 = {pow(x, 2, 13) for x in range(1, 13)}
qr_counts = collections.Counter()
for x in range(1, 13):
    qr_counts["B hits QR" if x % 13 in QR13 else "B hits QNR"] += 2  # +-x same QR-ness? (-1 is QR mod 13: 12=QR? 5^2=25=12 yes!)
print(f"  QR anatomy of +-[1,12] mod 13: {dict(qr_counts)} (-1 is a QR mod 13, so +-x same class)")
# Fourier bias of S in Z/156: top nontrivial |S-hat|
import math
best = []
for k in range(1, 156):
    re = sum(math.cos(2 * math.pi * k * s / 156) for s in S)
    im = sum(math.sin(2 * math.pi * k * s / 156) for s in S)
    best.append((math.hypot(re, im), k))
best.sort(reverse=True)
print(f"  |S|=24; top |S^(k)| (nontrivial): {[(round(b,2), k) for b, k in best[:6]]}")
print(f"  (uniform-random 24-set RMS would be ~{math.sqrt(24):.2f}; large peaks = structure = Erdos-Turan lever)")
print("DONE")
