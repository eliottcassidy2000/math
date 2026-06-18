# lrc14_refute_danger-interval_M4sectionrotatio_kps-S2-wf.py
# ADVERSARIAL REFUTATION of the "danger-interval / M4 section rotational" forbidden-class claim.
#
# CLAIM UNDER TEST:
#   M4 = section rotational map. At grid time tau=a/14, gcd(a,14)=1, runner i sits in
#   SECTION r_i=(v_i*a) mod 14. Arc i->j iff (r_i-r_j) mod 14 in {1,2,3,4,5,6};
#   d=0 or d=7 broken by smaller speed (i<j => i->j when speeds[i]<speeds[j]).
#   CLAIMED: at n=5 this FORBIDS the maximal-H non-rotational class (H=15, c3=4,
#   score=(1,2,2,2,3)) over LRC-constrained inputs, even at large vmax; structural
#   cause = M4 is an induced subtournament of the C_14 rotational tournament R_14.
#
# REFUTATION GOAL: realize that class even ONCE with this exact map over a broad
# exact search => holds=False. If we cannot, report search bound + whether the
# abstract ceiling (over ALL residue multisets incl collisions & all tie-breaks)
# even permits the class.
#
# Everything exact (Fraction). No floats for decisions.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, combinations_with_replacement, product
import sys

# ---------------- EXACT M TOOL (verbatim) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def Mfull(S):
    b = F(0); ats = []
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; ats = [t]
        elif v == b: ats.append(t)
    return b, sorted(ats)

def gcd_list(xs):
    g0 = 0
    for x in xs: g0 = gcd(g0, x)
    return g0

# ---------------- tournament invariants (n=5 fixed) ----------------
PERM5 = list(permutations(range(5)))
def canon5(adj):
    best = None
    for p in PERM5:
        bits = 0
        for a in range(5):
            for b in range(5):
                if a != b:
                    bits = (bits << 1) | (1 if adj[p[a]][p[b]] else 0)
        if best is None or bits < best: best = bits
    return best
def score5(adj):
    return tuple(sorted(sum(1 for j in range(5) if j != i and adj[i][j]) for i in range(5)))
def c3_5(adj):
    c = 0
    for a, b, cc in combinations(range(5), 3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        elif adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def h5(adj):
    return sum(1 for p in PERM5 if all(adj[p[k]][p[k+1]] for k in range(4)))

def invariants5(adj):
    return (h5(adj), c3_5(adj), score5(adj))

TARGET = (15, 4, (1, 2, 2, 2, 3))  # the claimed forbidden maximal-H non-rotational class

# ---------------- section-rotational M4 (EXACT, matches method4_adj) ----------------
def m4_section(speeds, a):
    """speeds: list of ints; a: grid time (gcd(a,14)=1). Returns adj on the given order."""
    m = len(speeds)
    r = [(speeds[i] * a) % 14 for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            d = (r[i] - r[j]) % 14
            if d == 0: adj[i][j] = (speeds[i] < speeds[j])
            elif 1 <= d <= 6: adj[i][j] = True
            elif d == 7: adj[i][j] = (speeds[i] < speeds[j])
            else: adj[i][j] = False
    return adj

# For OFF-GRID: there is no canonical "section". The map M4 is defined ON the grid
# (tau = a/14). The honest off-grid analogue: round tau to the nearest grid line, OR
# use the section the runner is CLOSEST to. We define an off-grid section as
# r_i = round(14 * frac(v_i*tau)) mod 14 (nearest of the 14 grid sections). We test
# this too, in case the optimal lonely tau is off-grid.
def frac(x):
    r = x - int(x); return r + 1 if r < 0 else r
def m4_section_offgrid(speeds, tau):
    m = len(speeds)
    r = []
    for v in speeds:
        # nearest integer to 14*frac(v*tau), exact via Fraction
        x = 14 * frac(v * tau)  # in [0,14)
        fl = int(x)
        # nearest int
        rr = fl if (x - fl) <= F(1, 2) else fl + 1
        r.append(rr % 14)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            d = (r[i] - r[j]) % 14
            if d == 0: adj[i][j] = (speeds[i] < speeds[j])
            elif 1 <= d <= 6: adj[i][j] = True
            elif d == 7: adj[i][j] = (speeds[i] < speeds[j])
            else: adj[i][j] = False
    return adj

# =====================================================================
# PHASE 0: ABSTRACT CEILING. Can the section-rotational map EVER produce TARGET,
# over ALL residue 5-multisets mod 14 (incl 0, incl collisions) and ALL tie-break
# speed orders? If NO, the class is forbidden for STRUCTURAL reasons independent of
# any LRC constraint -> the claim's mechanism is sound and refutation is impossible
# in principle (for this map). If YES, we hunt for an LRC realization.
# =====================================================================
def phase0():
    print("="*72)
    print("PHASE 0: abstract ceiling -- ALL residue 5-multisets (incl 0 & collisions),")
    print("         ALL tie-break orders. Can section-M4 produce TARGET at all?")
    print("="*72)
    reached = {}
    target_witnesses = []
    # tie-break only matters within equal-residue groups and across d=7 pairs.
    # We realize ALL tie-break orders by trying ALL distinct speed orderings (perm of 1..5).
    # SPEEDUP: dedupe by adjacency BIT PATTERN before the expensive canonicalization.
    tieperms = list(set(permutations(range(5))))
    cnt = 0
    seen_bits = set()
    def bits_of(adj):
        b = 0
        for a in range(5):
            for c in range(5):
                if a != c:
                    b = (b << 1) | (1 if adj[a][c] else 0)
        return b
    for res in combinations_with_replacement(range(14), 5):
        for perm in tieperms:
            speeds = [0]*5
            for rank, pos in enumerate(perm): speeds[pos] = rank + 1
            adj = _adj_from_res(list(res), speeds)
            bb = bits_of(adj)
            cnt += 1
            if bb in seen_bits:
                continue
            seen_bits.add(bb)
            inv = invariants5(adj)
            k = canon5(adj)
            if k not in reached: reached[k] = (inv, tuple(res), tuple(speeds))
            if inv == TARGET and len(target_witnesses) < 3:
                target_witnesses.append((tuple(res), tuple(speeds)))
    print(f"  scanned {cnt} (residue-multiset, tiebreak) pairs; {len(seen_bits)} distinct labeled tournaments")
    print(f"  distinct iso classes reachable: {len(reached)}")
    target_reachable = any(v[0] == TARGET for v in reached.values())
    print(f"  TARGET (H=15,c3=4,score=(1,2,2,2,3)) reachable abstractly? {target_reachable}")
    if target_witnesses:
        for res, sp in target_witnesses:
            print(f"     witness residues={res}, speed-order={sp}")
    return target_reachable

def _adj_from_res(res, speeds):
    m = len(res)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            d = (res[i] - res[j]) % 14
            if d == 0: adj[i][j] = (speeds[i] < speeds[j])
            elif 1 <= d <= 6: adj[i][j] = True
            elif d == 7: adj[i][j] = (speeds[i] < speeds[j])
            else: adj[i][j] = False
    return adj

# =====================================================================
# PHASE 1: LRC-constrained on-grid realization hunt. Broad primitive speed sets.
# For each primitive 5-set, for each grid time a in (Z/14)*, build section-M4 and
# check if it == TARGET. We do NOT restrict a to be a lonely time (the claim's own
# evidence scans ALL a in (Z/14)*). We report whether TARGET ever appears.
# =====================================================================
def phase1(vmaxes):
    print("="*72)
    print("PHASE 1: LRC on-grid hunt. Primitive 5-sets, all a in (Z/14)*.")
    print("="*72)
    units = [a for a in range(1, 14) if gcd(a, 14) == 1]
    for vmax in vmaxes:
        sets = [c for c in combinations(range(1, vmax+1), 5) if gcd_list(c) == 1]
        realized = {}
        target_hit = None
        for S in sets:
            S = list(S)
            for a in units:
                adj = m4_section(S, a)
                inv = invariants5(adj)
                k = canon5(adj)
                if k not in realized: realized[k] = inv
                if inv == TARGET and target_hit is None:
                    target_hit = (tuple(S), a)
        forb = 12 - len(realized)
        msg = "*** TARGET REALIZED ***" if target_hit else "target NOT realized"
        print(f"  vmax={vmax}: {len(sets)} primitive 5-sets, realized {len(realized)}/12, "
              f"forbidden={forb}.  {msg}")
        if target_hit:
            print(f"     WITNESS: S={target_hit[0]}, a={target_hit[1]}")

# =====================================================================
# PHASE 2: OFF-GRID hunt at the TRUE optimal lonely tau (and all optimal taus).
# THM-524: optimum is at a binding-pair crossing. The optimal tau is generally
# OFF the 14-grid. Use the nearest-section off-grid M4. Broad primitive sets,
# including covering/sporadic ones.
# =====================================================================
def phase2(vmaxes):
    print("="*72)
    print("PHASE 2: OFF-GRID at optimal lonely tau (nearest-section M4).")
    print("="*72)
    for vmax in vmaxes:
        sets = [c for c in combinations(range(1, vmax+1), 5) if gcd_list(c) == 1]
        realized = {}
        target_hit = None
        for S in sets:
            S = list(S)
            gap, taus = Mfull(S)
            for tau in taus:
                adj = m4_section_offgrid(S, tau)
                inv = invariants5(adj)
                k = canon5(adj)
                if k not in realized: realized[k] = inv
                if inv == TARGET and target_hit is None:
                    target_hit = (tuple(S), tau, gap)
        msg = "*** TARGET REALIZED ***" if target_hit else "target NOT realized"
        print(f"  vmax={vmax}: {len(sets)} sets, realized {len(realized)}/12.  {msg}")
        if target_hit:
            print(f"     WITNESS: S={target_hit[0]}, tau={target_hit[1]}, gap={target_hit[2]}")

# =====================================================================
# PHASE 3: OFF-GRID at ARBITRARY tau (not just the optimum). The map M4 is only
# *defined* on-grid, but the off-grid nearest-section extension is well-defined for
# ANY tau. Scan many tau values per set to maximize chances. This is the most
# generous off-grid search.
# =====================================================================
def phase3(vmax, tau_den_max):
    print("="*72)
    print(f"PHASE 3: OFF-GRID arbitrary tau (den<= {tau_den_max}), vmax={vmax}.")
    print("="*72)
    sets = [c for c in combinations(range(1, vmax+1), 5) if gcd_list(c) == 1]
    realized = {}
    target_hit = None
    taus = []
    for den in range(2, tau_den_max+1):
        for num in range(1, den):
            t = F(num, den)
            if t < F(1, 2) or t == F(1, 2):
                taus.append(t)
    taus = sorted(set(taus))
    for S in sets:
        S = list(S)
        for tau in taus:
            adj = m4_section_offgrid(S, tau)
            inv = invariants5(adj)
            k = canon5(adj)
            if k not in realized: realized[k] = inv
            if inv == TARGET and target_hit is None:
                target_hit = (tuple(S), tau)
        if target_hit: break
    msg = "*** TARGET REALIZED ***" if target_hit else "target NOT realized"
    print(f"  {len(sets)} sets x {len(taus)} taus, realized {len(realized)}/12.  {msg}")
    if target_hit:
        print(f"     WITNESS: S={target_hit[0]}, tau={target_hit[1]}")

# =====================================================================
# PHASE 4: ALL UNITS / SPORADIC / COVERING / TIGHT-style inputs at n=5 analogues.
# Hand-picked LRC-hard 5-sets: AP {1,2,3,4,5}, covering-flavored, with a "parked"
# multiple of 14, sporadic. On-grid all a + off-grid optimum.
# =====================================================================
def phase4():
    print("="*72)
    print("PHASE 4: sporadic / covering / parked / tight 5-set probes.")
    print("="*72)
    probes = [
        (1,2,3,4,5),        # AP
        (1,2,3,4,6),
        (1,2,4,5,7),
        (1,3,4,5,7),
        (1,2,3,5,8),
        (1,2,3,4,14),       # parked multiple of 14
        (1,2,3,5,14),
        (1,2,4,7,14),
        (1,5,7,11,13),
        (2,3,5,7,11),       # primes
        (1,2,3,4,28),       # parked 2*14
        (1,3,5,9,11),       # QR-ish mod 14 units
        (1,2,6,9,13),
        (3,5,6,7,8),
        (1,4,6,9,11),
    ]
    units = [a for a in range(1, 14) if gcd(a, 14) == 1]
    any_hit = False
    for S in probes:
        if gcd_list(S) != 1: continue
        Sl = list(S)
        hits = []
        # on-grid all a
        for a in units:
            adj = m4_section(Sl, a)
            if invariants5(adj) == TARGET: hits.append(("grid", a))
        # off-grid optimum
        gap, taus = Mfull(Sl)
        for tau in taus:
            adj = m4_section_offgrid(Sl, tau)
            if invariants5(adj) == TARGET: hits.append(("opt", tau))
        flag = "  <-- TARGET HIT" if hits else ""
        if hits: any_hit = True
        print(f"  S={S}: gap={gap}  hits={hits}{flag}")
    print(f"  any TARGET hit among probes? {any_hit}")

if __name__ == "__main__":
    ceiling = phase0()
    print()
    if ceiling:
        print(">> Abstract ceiling ALLOWS the target. Hunting for LRC realization...\n")
    else:
        print(">> Abstract ceiling FORBIDS the target structurally (no residue multiset, no")
        print(">> tie-break can make it). LRC realization is impossible for THIS map.\n")
    phase1([10, 13, 16, 20])
    print()
    phase2([10, 13, 16])
    print()
    phase3(11, 30)
    print()
    phase4()
