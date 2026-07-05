#!/usr/bin/env python3
"""
opus-2026-07-05-S82 -- HYP-4119 part 2: THE ORBIT CENSUS of full covers at level 169.

The unit group U(169) (156 elements) acts on families by dilation v -> u*v mod 169;
on patterns: position r moves to u*r mod 13, digit kappa_r moves to
digit(u*(r+13*kappa_r) mod 169).  Cover-hood is invariant (witness a for u*v is
u^{-1}*a for v).  So the full-cover census lives on ORBITS (size <= 156), and one
kernel row per orbit certifies 156 pattern classes (the witness dilates).

This DFS enumerates full covers DIRECTLY in canonical form: branch position-by-position
in a fixed order with the column-deficit prune, collect full 12-digit patterns whose
canonical form is new.  To keep the tree finite we enumerate CORES (irredundant covers,
cell-driven) and canonicalize completions lazily with a global seen-set of canonical
forms; memory = #orbits.
"""
import sys, time, json, collections
from itertools import product

MOD = 169
CELLS = [a for a in range(1, 169) if a % 13 != 0]
CIDX = {a: i for i, a in enumerate(CELLS)}
FULL = (1 << 156) - 1

def bad(x): return x <= 13 or x >= 156

KM = [[0] * 13 for _ in range(13)]
for r in range(1, 13):
    for kap in range(13):
        v = r + 13 * kap
        m = 0
        for a in CELLS:
            if bad((a * v) % MOD):
                m |= 1 << CIDX[a]
        KM[r][kap] = m

def cover_mask(pat):
    m = 0
    for r in range(1, 13):
        m |= KM[r][pat[r - 1]]
    return m

# ---- the unit action on patterns
UNITS = [u for u in range(1, MOD) if u % 13 != 0]
ACTION = {}  # ACTION[u] = list over r=1..12 of (new_position, function kappa -> new_kappa) as table
for u in UNITS:
    tab = []
    for r in range(1, 13):
        s = (u * r) % 13
        row = []
        for kap in range(13):
            val = (u * (r + 13 * kap)) % MOD
            assert val % 13 == s
            row.append(((val - s) // 13) % 13)
        tab.append((s, row))
    ACTION[u] = tab

def transform(pat, u):
    out = [0] * 12
    for r in range(1, 13):
        s, row = ACTION[u][r - 1]
        out[s - 1] = row[pat[r - 1]]
    return tuple(out)

def canon(pat):
    return min(transform(pat, u) for u in UNITS)

# ---- shadows in canonical form
def shadow(lam):
    laminv = pow(lam % 13, -1, 13)
    return tuple((((lam * ((laminv * r) % 13)) % MOD - r) // 13) % 13 for r in range(1, 13))

sh_canon = {canon(shadow(lam)) for lam in UNITS}
print(f"shadow orbits: {len(sh_canon)} (156 shadows collapse to this many orbits)", flush=True)

# ---- cell-driven DFS over cores; canonicalize completions on the fly
killers = {a: [] for a in CELLS}
for r in range(1, 13):
    for kap in range(13):
        m = KM[r][kap]
        for a in CELLS:
            if m >> CIDX[a] & 1:
                killers[a].append((r, kap))

COLMASK = [0] * 13
for a in CELLS:
    COLMASK[a % 13] |= 1 << CIDX[a]

sys.setrecursionlimit(10000)
t0 = time.time()
nodes = 0
cores_seen = set()
orbit_covers = set()

def process_core(assigned):
    """A core covers fully; every completion is a cover. Canonicalize all completions."""
    key = tuple(sorted(assigned.items()))
    if key in cores_seen:
        return
    cores_seen.add(key)
    free = [r for r in range(1, 13) if r not in assigned]
    base = [assigned.get(r, 0) for r in range(1, 13)]
    for combo in product(range(13), repeat=len(free)):
        pat = list(base)
        for r, kv in zip(free, combo):
            pat[r - 1] = kv
        orbit_covers.add(canon(tuple(pat)))

def dfs(assigned, covered):
    global nodes
    nodes += 1
    if nodes % 2000000 == 0:
        print(f"  ...{nodes} nodes, {len(cores_seen)} cores, {len(orbit_covers)} orbits, "
              f"{time.time()-t0:.0f}s", flush=True)
    if covered == FULL:
        process_core(assigned)
        return
    nrem = 12 - len(assigned)
    for b in range(1, 13):
        if 13 - bin(covered & COLMASK[b]).count("1") > 2 * nrem:
            return
    best_opts = None
    for a in CELLS:
        if covered >> CIDX[a] & 1:
            continue
        opts = [(r, kap) for (r, kap) in killers[a] if r not in assigned]
        if best_opts is None or len(opts) < len(best_opts):
            best_opts = opts
            if len(opts) <= 2:
                break
    if not best_opts:
        return
    for (r, kap) in best_opts:
        assigned[r] = kap
        dfs(assigned, covered | KM[r][kap])
        del assigned[r]

dfs({}, 0)
el = time.time() - t0
print(f"DFS complete: {nodes} nodes, {el:.1f}s", flush=True)
print(f"distinct cores: {len(cores_seen)}; distinct cover ORBITS: {len(orbit_covers)}")

non_shadow = orbit_covers - sh_canon
print(f"shadow orbits among covers: {len(orbit_covers & sh_canon)}; NON-SHADOW orbits: {len(non_shadow)}")

hist = collections.Counter(sum(1 for x in p if x) for p in orbit_covers)
print(f"visible-count histogram over cover-orbit REPRESENTATIVES: {dict(sorted(hist.items()))}")
# NOTE: visible count is NOT orbit-invariant (kappa=0 can move); compute orbit-max visible:
def orbit_max_visible(p):
    return max(sum(1 for x in transform(p, u) if x) for u in UNITS)
def orbit_min_visible(p):
    return min(sum(1 for x in transform(p, u) if x) for u in UNITS)
h2 = collections.Counter(orbit_min_visible(p) for p in orbit_covers)
print(f"orbit MIN-visible histogram: {dict(sorted(h2.items()))}")

with open("05-knowledge/results/shadow_orbit_census_opus_S82.json", "w") as f:
    json.dump({"orbits": sorted(map(list, orbit_covers)),
               "shadows": sorted(map(list, sh_canon)),
               "non_shadow": sorted(map(list, non_shadow))}, f)
print("orbit lists dumped to 05-knowledge/results/shadow_orbit_census_opus_S82.json")

# small-rep viability for non-shadow orbits: does the small representative
# (v_r = r + 13 kappa_r <= 168) have a strict witness at SOME denominator <= 336
# (merge-grid attainment: M attained at m/(vi+vj))?
from fractions import Fraction as F
def exact_loose_check(pat):
    v = [r + 13 * pat[r - 1] for r in range(1, 13)]
    best = F(0)
    seen = set()
    for i in range(12):
        for j in range(i, 12):
            s = v[i] + v[j]
            for m in range(1, s):
                t = F(m, s)
                if t in seen or t > F(1, 2):
                    continue
                seen.add(t)
                q = min(min((vi * t) % 1, 1 - (vi * t) % 1) for vi in v)
                if q > best:
                    best = q
    return best, v

checked = 0
for p in sorted(non_shadow)[:12]:
    b, v = exact_loose_check(p)
    tag = "LOOSE (row at merge-grid denom)" if b > F(1, 13) else "*** TIGHT?!"
    print(f"  non-shadow small-rep {v}: M = {b} {tag}", flush=True)
    checked += 1
print(f"small-rep check done on {checked} non-shadow orbits; total {time.time()-t0:.1f}s")
