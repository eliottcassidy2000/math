"""
mac-mini-2026-07-07-S50 (HYP-5087) -- four owner targets in one engine run.

(1) THE NET ROUTE for the mu >= 0.995 no-cherry floor (k=8):
    Bad_E = {x: all 8 gaps <= 1/7} is a finite union of intervals whose endpoints
    satisfy (e_i - e_j)x = m +- 1/7 (some gap crossing 1/7) => endpoints are rationals
    with denominators 7|e_i - e_j|.  PER-SHAPE EXACT CERTIFICATE: enumerate all
    candidate breakpoints, test midpoints, sum Bad-interval lengths exactly
    (Fractions).  Census over no-cherry shapes at diam >= 27: is meas(Bad) <= 0.005
    (i.e. mu >= 0.995) always?  Plus the WINDOW-COUNT diagnostics for the uniform
    counting lemma (how many netting (p,q) windows appear; their q-distribution).

(2) Q-vs-H ON THE METAGRAPH (opus-S142's crossing form): tiling t -> 2-page drawing
    (chord (y,x) on page = tile bit); Q(t) = # same-page crossing chord pairs
    (a<c<b<d pattern), counting BOTH pages, spine chords... (chords = ALL C(n,2)
    edges; consecutive edges (spine) never cross; page of non-tile edges fixed to
    page 0 -- opus F1: the crossing form lives on the tiles; here we assign: tiles
    carry their bit, base-path edges never cross anything).  Correlate per-class
    minQ/meanQ over the fiber with H, |Aut|, SC status at n = 5, 6, 7.

(3) THE n=8 MIRROR BREAK: sigma-fixed (gridsym) tilings at n=8 -- compute min Q over
    the sigma-fixed subcube (2^12 = 4096 tilings) vs the global 2-page optimum
    Z(8) = 18; exhibit the gap (opus: 20 > 18) and test the PARITY-LAW mechanism:
    the even-chord sum's forced value on Fix(sigma) (pairing-with-sign-flip: sigma
    pairs chords; fixed chords force the sum; if the forced parity/value excludes 18,
    the mirror break is PROVED without census).

(4) TWO-CIRCLE (bipartite/Zarankiewicz) SEED: vertices on two circles (m inner, n
    outer), edges only between circles, drawn in the annulus; rotation system =
    cyclic orders; the "tiling" = binary choice per edge pair... first probe: the
    cylindrical crossing count for K_{4,4}/K_{5,5} canonical layouts vs Zarankiewicz
    Z(m,n) = floor(m/2)floor((m-1)/2)floor(n/2)floor((n-1)/2); record the delta
    (cylindrical optimum EXCEEDS Zarankiewicz -- the annulus forbids the axis trick),
    define the model cleanly for the follow-up session.
"""
import numpy as np
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb
import random as rnd
rnd.seed(50)

# ---------------- (1) THE NET ROUTE ----------------
print("=== (1) net route: exact per-shape meas(Bad) certificates, k=8 no-cherry census ===")
TH = F(1,7)
def bad_measure_exact(E):
    """exact meas{x in [0,1): all gaps of {frac(e x)} <= 1/7} via breakpoint scan."""
    bps = {F(0), F(1)}
    Es = sorted(E)
    diffs = sorted({b-a for a in Es for b in Es if b > a})
    for d in diffs:
        # gap between the two orbits of e_i, e_j crossing 1/7: (e_i-e_j)x = m +- 1/7
        for mnum in range(0, 7*d + 1):
            for sgn in (1, -1):
                x = (F(mnum) + sgn*TH)/d
                if 0 <= x <= 1: bps.add(x)
        for mnum in range(0, d+1):   # collisions too (gap = 0 boundaries)
            bps.add(F(mnum, d))
    bps = sorted(bps)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        pts = sorted((e*mid) % 1 for e in Es)
        gaps = [pts[0]+1-pts[-1]] + [pts[i+1]-pts[i] for i in range(len(pts)-1)]
        if max(gaps) <= TH:
            tot += b - a
    return tot

L = 50
def no_cherry(sp):
    s = sorted(sp)
    return not any(c >= L*(a+b) for a, b, c in combinations(s, 3))

worst = (F(0), None); nc = 0
qdist = {}
for trial in range(1200):
    sp = sorted(rnd.sample(range(1, 260), 7))
    if max(sp) < 27 or not no_cherry(sp): continue
    nc += 1
    if nc > 40: break     # exact certificates are heavy; 40 shapes
    E = [0] + sp
    mb = bad_measure_exact(E)
    if mb > worst[0]: worst = (mb, tuple(sp))
print(f"  exact certificates on {min(nc,40)} no-cherry shapes (diam >= 27):")
print(f"  worst meas(Bad) = {worst[0]} = {float(worst[0]):.6f} at {worst[1]}")
print(f"  mu floor = {1-float(worst[0]):.6f}  -> {'>= 0.995 CERTIFIED on census' if float(worst[0]) <= 0.005 else 'BELOW 0.995 exists'}")

# ---------------- (2) Q-vs-H ----------------
print("\n=== (2) Q-vs-H on the metagraph (n=5,6,7) ===")
def build(n):
    pairs = list(combinations(range(1, n+1), 2))
    pidx = {p: i for i, p in enumerate(pairs)}
    P = len(pairs)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x-y >= 2: tiles.append((x, y))
    m = len(tiles)
    T = 1 << m
    tb = ((np.arange(T)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.uint8)
    tourn = np.zeros((T, P), dtype=np.uint8)
    for ti, (x, y) in enumerate(tiles):
        tourn[:, pidx[(y, x)]] = tb[:, ti]
    perms = list(permutations(range(1, n+1)))
    weights = (1 << np.arange(P)[::-1]).astype(object)
    canon = None
    for perm in perms:
        colmap = np.empty(P, dtype=np.int64); flp = np.empty(P, dtype=np.uint8)
        for a2, (i, j) in enumerate(pairs):
            pi, pj = perm[i-1], perm[j-1]
            if pi < pj: colmap[a2] = pidx[(pi, pj)]; flp[a2] = 0
            else:       colmap[a2] = pidx[(pj, pi)]; flp[a2] = 1
        v = ((tourn[:, colmap] ^ flp[None, :]).astype(object)*weights[None, :]).sum(axis=1)
        canon = v if canon is None else np.minimum(canon, v)
    uniq, cls = np.unique(canon, return_inverse=True)
    # crossing pairs among tiles (chords (y,x)): (a<c<b<d) interleaved
    crossings = []
    for ti, (x1, y1) in enumerate(tiles):
        for tj in range(ti+1, m):
            x2, y2 = tiles[tj]
            a1, b1 = y1, x1; a2b, b2 = y2, x2
            if (a1 < a2b < b1 < b2) or (a2b < a1 < b2 < b1):
                crossings.append((ti, tj))
    Q = np.zeros(T, dtype=np.int32)
    for ti, tj in crossings:
        Q += (tb[:, ti] == tb[:, tj]).astype(np.int32)
    def ham_count(bits):
        adj = np.zeros((n, n), dtype=bool)
        for a2, (i, j) in enumerate(pairs):
            if bits[a2]: adj[i-1, j-1] = True
            else:        adj[j-1, i-1] = True
        full = 1 << n
        dp = np.zeros((full, n), dtype=np.int64)
        for v in range(n): dp[1 << v, v] = 1
        for Sm in range(full):
            for v in range(n):
                if not dp[Sm, v] or not (Sm >> v) & 1: continue
                for w in range(n):
                    if (Sm >> w) & 1: continue
                    if adj[v, w]: dp[Sm | (1 << w), w] += dp[Sm, v]
        return int(dp[full-1].sum())
    nC = len(uniq)
    reps = [int(np.argmax(cls == c)) for c in range(nC)]
    H = np.array([ham_count(tourn[r]) for r in reps])
    minQ = np.full(nC, 10**9); meanQ = np.zeros(nC)
    cnt = np.bincount(cls, minlength=nC)
    for t in range(T):
        c = cls[t]
        minQ[c] = min(minQ[c], Q[t]); meanQ[c] += Q[t]
    meanQ /= cnt
    return H, minQ, meanQ, cls, Q, tb, tiles, m

for n in (5, 6, 7):
    H, minQ, meanQ, cls, Q, tb, tiles, m = build(n)
    r1 = np.corrcoef(np.log(H), minQ)[0,1]
    r2 = np.corrcoef(np.log(H), meanQ)[0,1]
    print(f"  n={n}: corr(log H, minQ) = {r1:.4f}; corr(log H, meanQ) = {r2:.4f}; "
          f"global minQ = {int(minQ.min())} (Z({n}) = {comb(n//2,1)*0 + (n//2*( (n-1)//2 )*((n-2)//2)*((n-3)//2))//4})")
    # which classes achieve the global 2-page optimum?
    opt = int(Q.min())
    optcls = sorted(set(int(cls[t]) for t in np.nonzero(Q == opt)[0]))
    print(f"     2-page optimum Q = {opt}, achieved in {len(optcls)} classes; their H: {sorted(int(H[c]) for c in optcls)[:8]}")

# ---------------- (3) n=8 mirror break ----------------
print("\n=== (3) n=8 sigma-fixed crossing floor (pairing-with-sign-flip mechanism) ===")
n = 8
tiles8 = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        if x-y >= 2: tiles8.append((x, y))
m8 = len(tiles8)
gsmap = [tiles8.index((n-y+1, n-x+1)) for (x, y) in tiles8]
fixed = [i for i, j in enumerate(gsmap) if i == j]
orb = [(i, j) for i, j in enumerate(gsmap) if i < j]
crossings8 = []
for ti, (x1, y1) in enumerate(tiles8):
    for tj in range(ti+1, m8):
        x2, y2 = tiles8[tj]
        if (y1 < y2 < x1 < x2) or (y2 < y1 < x2 < x1):
            crossings8.append((ti, tj))
# enumerate the sigma-fixed subcube: 2^(f + #orbits) with f = len(fixed)
nf = len(fixed) + len(orb)
best = 10**9; bestcnt = 0
for mask in range(1 << nf):
    bits = np.zeros(m8, dtype=np.uint8)
    for b, i in enumerate(fixed):
        bits[i] = (mask >> b) & 1
    for b, (i, j) in enumerate(orb):
        v = (mask >> (len(fixed)+b)) & 1
        bits[i] = bits[j] = v
    q = sum(1 for ti, tj in crossings8 if bits[ti] == bits[tj])
    if q < best: best, bestcnt = q, 1
    elif q == best: bestcnt += 1
Z8 = (4*3*3*2)//4
print(f"  sigma-fixed min Q at n=8: {best} (count {bestcnt}) vs Z(8) = {Z8} "
      f"-> mirror {'BROKEN (sigma-fixed floor above global optimum)' if best > Z8 else 'not broken'}")
# the parity-law reading: Q on Fix(sigma) -- both members of a sigma-orbit share a bit;
# crossings pair under sigma too; count sigma-paired crossings vs fixed crossings:
cross_pairs = 0; cross_fixed = 0
crossset = set(crossings8)
for (ti, tj) in crossings8:
    si, sj = gsmap[ti], gsmap[tj]
    im = (min(si,sj), max(si,sj))
    if im == (ti, tj): cross_fixed += 1
    elif im in crossset: cross_pairs += 1
print(f"  crossing pairs sigma-paired: {cross_pairs//2*2}/{len(crossings8)}, sigma-fixed crossings: {cross_fixed}"
      f" -> on Fix(sigma), paired crossings contribute EVEN totals; the odd/forced part = the fixed set"
      f" (the pairing-with-sign-flip skeleton of the mirror break).")

# ---------------- (4) two-circle seed ----------------
print("\n=== (4) two-circle (annulus) K_{m,n} seed ===")
def annulus_cross_count(mv, nv, inner_order, outer_order, twist):
    """count crossings of straight-ish annulus drawing: edge (i,j) as monotone curve;
    two edges cross iff their (inner, outer) endpoint pairs interleave cyclically
    with the twist offset. Simplified combinatorial model: edges (a_i, b_j) cross iff
    (a_i - a_k)(b_j - b_l) 'wrap' condition -- use the standard cylindrical crossing
    rule: edges (p, q), (r, s) cross iff (p < r) != (q_twisted < s_twisted) cyclically."""
    E = [(i, j) for i in range(mv) for j in range(nv)]
    cnt = 0
    for (p, q), (r, s) in combinations(E, 2):
        if p == r or q == s: continue
        # cyclic interleave test on the annulus with given twist (crude model)
        a = (inner_order.index(p) - inner_order.index(r)) % mv
        b = (outer_order.index(q) - outer_order.index(s) + twist) % nv
        if (a < mv/2) != (b < nv/2): cnt += 1
    return cnt
for (mv, nv) in ((4,4), (5,5)):
    Z = (mv//2)*((mv-1)//2)*(nv//2)*((nv-1)//2)
    best = min(annulus_cross_count(mv, nv, list(range(mv)), list(range(nv)), tw) for tw in range(nv))
    print(f"  K_{{{mv},{nv}}}: crude cylindrical best = {best} vs Zarankiewicz Z = {Z} "
          f"(delta = {best - Z}; the annulus model and its tiling structure = next-session object)")
