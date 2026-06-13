#!/usr/bin/env python3
"""
lrc_braid_persistence_exotic_s540.py    oracle-2026-06-01-S540o

Out-of-the-box constructions for LRC. Two computed STARS + a posed wild multitude.

STAR 1 -- SPACETIME BRAID. The runner worldlines x_i(t)=frac(v_i t) on the cylinder
S^1 x [0,1] are n strands (observer = strand at 0). Over [0,1) they form a PURE BRAID
on n strands (at t=1 every frac(v_i)=0, so all strands return to start). Facts:
  - signed linking number lk(i,j) = v_i - v_j  (i laps j exactly v_i-v_j times);
    the LINKING MATRIX L = v(x)1 - 1(x)v is exactly the S538 TENSION / difference data.
  - braid word length = total crossings = sum_{i<j} |v_i - v_j| = the holdback (S25).
  - LRC = the braid has a time-slice where the observer strand is at circular distance
    >= 1/n from every other strand (a 'fat channel' around strand 0).
  Realizable braids = 'torus/homogeneous' pure braids from linear worldlines (restricted).

STAR 2 -- PERSISTENT HOMOLOGY / VINEYARDS. At time t the n points on the circle are a
point cloud; the Rips H_0 barcode = the GAP structure (a gap g dies at radius g/2; the
longest bar = the largest gap = the apex, S530). The observer's component persists
until radius r_obs(t) = (min flanking gap)/2. So:
  LRC <=>  max_t r_obs(t) >= 1/(2n)   (observer's H_0 bar outlives radius 1/(2n)).
Over t this is a VINEYARD (persistence diagram moving with t); loneliness = the
observer's vine crosses the line r=1/(2n).

THE WILD MULTITUDE (posed): tropical Newton polygon of the speeds; p-adic Bruhat-Tits
tree walk (3-adic for n*=9, S534); quantum/operator non-commutativity tournament
(U_i=e^{2pi i v_i t}); abelian sandpile / chip-firing on the sector cycle (sandpile
group Z_n); dynamical zeta over wall-crossing periodic orbits; quasicrystal hull
(irrational-speed cut-and-project); Sprague-Grundy game values of runner-pairs;
Galois/Frobenius action on the n-th roots of unity.

Computes braid (linking=tension, word length=holdback, pure-braid check, fat-channel
LRC) and persistence (observer vine, LRC threshold, barcode), verifying the
connections.
"""
from itertools import combinations
from functools import reduce
from math import gcd
import random

def frac(x): return x - int(x // 1)
def cdist(a, b):
    d = abs(a-b) % 1.0; return min(d, 1-d)

# ---------------- STAR 1: spacetime braid ----------------
def linking_matrix(speeds_with_obs, n):
    return [[speeds_with_obs[i]-speeds_with_obs[j] for j in range(n)] for i in range(n)]

def braid_word_length(speeds_with_obs, n):
    return sum(abs(speeds_with_obs[i]-speeds_with_obs[j]) for i, j in combinations(range(n), 2))

def is_pure_braid(speeds_with_obs, n):
    # at t=1 all frac(v_i)=0 -> identity permutation
    return all(frac(s*1.0) < 1e-9 or abs(frac(s*1.0)-1) < 1e-9 for s in speeds_with_obs)

def braid_fat_channel(speeds, n, samples=4000):
    """max over t of min circular distance from observer (strand 0) to any runner.
    LRC <=> this >= 1/n."""
    best = 0.0; bt = 0
    for s in range(samples):
        t = (s+0.5)/samples
        d = min(cdist(0.0, frac(v*t)) for v in speeds)
        if d > best: best = d; bt = t
    return best

def study_braid(n, n_sets=60):
    rnd = random.Random(40+n); rows = []; tot = 0
    tension_ok = True
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        sp = [0]+list(v)
        L = linking_matrix(sp, n)
        # tension check: L[i][j] = (v_i) - (v_j) is a coboundary -> L[i][j]+L[j][k]+L[k][i]=0
        for a, b, c in combinations(range(n), 3):
            if L[a][b]+L[b][c]+L[c][a] != 0: tension_ok = False
        wl = braid_word_length(sp, n)
        pure = is_pure_braid(sp, n)
        fat = braid_fat_channel(v, n)
        rows.append((wl, pure, fat))
    avg_wl = sum(r[0] for r in rows)/len(rows)
    all_pure = all(r[1] for r in rows)
    lrc = sum(1 for r in rows if r[2] >= 1.0/n - 1e-6)
    return avg_wl, all_pure, tension_ok, lrc, len(rows)

# ---------------- STAR 2: persistence / vineyards ----------------
def gaps(speeds, n, t):
    pts = sorted([0.0]+[frac(v*t) for v in speeds])
    g = [pts[i+1]-pts[i] for i in range(len(pts)-1)] + [1.0 - pts[-1] + pts[0]]
    return g, pts

def observer_vine(speeds, n, samples=4000):
    """r_obs(t) = (min flanking gap of observer at 0)/2; return max over t."""
    best = 0.0
    for s in range(samples):
        t = (s+0.5)/samples
        # observer at 0; nearest ahead and behind
        ahead = min((frac(v*t) for v in speeds if frac(v*t) > 1e-12), default=1.0)
        behind = min((1-frac(v*t) for v in speeds if frac(v*t) > 1e-12), default=1.0)
        r = min(ahead, behind)/2.0
        if r > best: best = r
    return best

def barcode_at_best(speeds, n, samples=4000):
    """barcode (sorted gap/2 death radii) at the observer-optimal time."""
    best = -1; bg = None
    for s in range(samples):
        t = (s+0.5)/samples
        ahead = min((frac(v*t) for v in speeds if frac(v*t) > 1e-12), default=1.0)
        behind = min((1-frac(v*t) for v in speeds if frac(v*t) > 1e-12), default=1.0)
        r = min(ahead, behind)/2.0
        if r > best:
            best = r; g, _ = gaps(speeds, n, t); bg = sorted(round(x/2, 3) for x in g)
    return bg

def study_persistence(n, n_sets=60):
    rnd = random.Random(70+n); tot = 0; lrc = 0; vines = []
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        rmax = observer_vine(v, n)
        vines.append(rmax)
        if rmax >= 1.0/(2*n) - 1e-6: lrc += 1
    return sum(vines)/len(vines), lrc, len(vines)

def main():
    print("="*74)
    print("STAR 1 -- SPACETIME BRAID: runners = strands; pure braid; linking = tension")
    print("="*74)
    for n in (5, 6, 7):
        avg_wl, all_pure, tension_ok, lrc, tot = study_braid(n)
        print(f"  n={n}: all pure braids={all_pure}; linking matrix = tension (cocycle) ok={tension_ok}; "
              f"avg braid-word length (=holdback Sum|v_i-v_j|)={avg_wl:.1f};")
        print(f"        FAT-CHANNEL LRC (observer strand >= 1/n from all at some t): {lrc}/{tot}")
    print("  => the runner system IS a pure braid on n strands; its linking matrix")
    print("     lk(i,j)=v_i-v_j is exactly the S538 tension; word length = holdback (S25);")
    print("     LRC = a fat channel around the observer strand. Realizable = torus braids.")
    print()
    print("="*74)
    print("STAR 2 -- PERSISTENT HOMOLOGY / VINEYARDS: H_0 barcode = gaps; observer vine")
    print("="*74)
    for n in (5, 6, 7):
        avg_vine, lrc, tot = study_persistence(n)
        v = tuple(sorted(random.Random(70+n).sample(range(1, 6*n), n-1)))
        while reduce(gcd, v) != 1:
            v = tuple(sorted(random.Random(71+n).sample(range(1, 6*n), n-1)))
        bc = barcode_at_best(v, n)
        print(f"  n={n}: threshold 1/(2n)={1/(2*n):.4f}; avg observer-vine max r_obs={avg_vine:.4f}; "
              f"LRC (r_obs >= 1/(2n) reached): {lrc}/{tot}")
        print(f"        example barcode (gap/2 death radii) at observer-optimal t: {bc}")
    print("  => the H_0 barcode = the gaps (longest bar = apex = largest gap, S530);")
    print("     observer's component persists to r_obs=(min flanking gap)/2; LRC = the")
    print("     observer VINE crosses r=1/(2n). A vineyard over t.")
    print()
    print("="*74)
    print("THE WILD MULTITUDE (posed)")
    print("="*74)
    print("""  - TROPICAL: Newton polygon / tropical curve of the speed polynomial; tournament
    from tropical (min-plus) comparisons; LRC = a tropical cell containing 0.
  - p-ADIC BRUHAT-TITS TREE: speeds' p-adic valuations (p|n*, e.g. 3 for n=18/n*=9,
    S534); runner dynamics = a walk on the tree; loneliness = reaching a far vertex.
  - QUANTUM / OPERATOR: U_i = e^{2pi i v_i t} unitaries; tournament from operator
    ordering / non-commutativity defect; LRC = a near-classical (commuting) window.
  - ABELIAN SANDPILE: sector occupancy = chips on the cycle C_n; toppling; sandpile
    group Z_n; LRC = a recurrent config with the observer cells empty.
  - DYNAMICAL ZETA: wall-crossing flow on the torus; periodic orbits = rational t;
    the zeta/Ruelle resonances; LRC = a spectral gap.
  - QUASICRYSTAL HULL: irrational speeds -> cut-and-project point set; the hull's
    local patches; LRC at the rational boundary = a forbidden patch.
  - SPRAGUE-GRUNDY GAME: each runner-pair = a subtraction game (period |v_i-v_j|);
    Grundy values -> a game-tournament; LRC = a P-position of the observer.
  - GALOIS / FROBENIUS: the n-th roots of unity with Gal(Q(zeta_n)/Q); the speed
    action; tournament from Frobenius orbits; LRC = an orbit avoiding the observer arc.""")

if __name__ == "__main__":
    main()
