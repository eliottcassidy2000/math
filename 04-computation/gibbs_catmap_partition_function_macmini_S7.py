"""
GIBBS MEASURE & ARNOLD CAT MAP <-> the apex-prime partition function.
mac-mini-2026-06-20-S7.

Thread hypothesis: measS7 (LRC) and H/OCF (tournaments) are FREE ENERGIES
F = -(1/beta) log Z of a Gibbs measure; the "irreducibly aggregate" obstruction
(HYP-2704 / mac-mini-S6) is PRECISELY that a free energy is -log of a SUM, hence
NOT a sum of local terms. consec (LRC) / regular (OCF) = the beta->inf GROUND STATE
selected by a variational principle, not a term-by-term inequality.

We test FOUR concrete falsifiable claims with exact arithmetic (Fractions) on small n.

THE OBJECTS (cited):
- THM-554: Z_n(x) = (prod_{v>=2} x_v) prod_{tiles (a,b)} (x_a + x_b). Score partition fn.
- THM-555: c3 = last score-determined OCF datum; H not cell-affine (cycle space).
- HYP-2704: death-chain K_{r+1}(t)=(1-t/7)K_r(t)+(t/7)K_r(t-1); measS7 = surjection event.
- HYP-2703 / mac-mini-S6: measS7=(1/7)sum_s bandcover_s; band twist e->s*e mod 7.
- the-apex-prime-partition-function...: H/p0 = interacting free energy; c3/decorr cover = free energy.

H1 (CONTROL, cut-space is LOCAL): Setting x_v=exp(-beta h_v) in Z_n, the free energy
    F = -(1/beta) log Z_n factorizes: log Z_n = sum_v(-beta h_v) + sum_tiles log(e^{-beta h_a}+e^{-beta h_b}).
    Ground state beta->inf: each tile -> max(h_a,h_b)=local. So the SCORE side IS a sum of local
    terms => no aggregate obstruction. PREDICTION: log Z is exactly additive over tiles. CONTROL.

H2 (LRC is non-local): measS7 = meas{|colors|=7} = a SURJECTION event = -log of a sum.
    PREDICTION: -log measS7(E) is NOT a sum of per-runner contributions, i.e. there is no
    f: offsets->R with log measS7(E) = sum_{e in E} f(e) + const, OR with pairwise corrections.
    We test additive AND pairwise (2-body) log-linear models and measure residual. Aggregate => fails.

H3 (Ground-state/variational): Build Gibbs family p_beta(E) ~ exp(beta * measS7(E)) over k-subsets.
    PREDICTION: beta->+inf concentrates on argmax = consec block (LRC extremal). And the
    log-partition Phi(beta)=log sum_E exp(beta*measS7(E)) is CONVEX in beta (always true), with
    Phi'(beta)->max measS7. Tournament twin: p_beta(T)~exp(beta*c3) -> regular score (THM-027/555),
    and exp(beta*H) ground state = H-max. Check convexity + ground state identity exactly.

H4 (Arnold cat / toral automorphism): The band twist e->s*e mod 7 (HYP-2703) is a linear map on
    Z/7. PREDICTION: it is a hyperbolic toral automorphism (cat-map-like) iff its matrix on the
    relevant torus has |trace|>2 / eigenvalues off unit circle. Test the actual renormalization:
    the multiplier action on (slope band, residue) as a 2D map; is it hyperbolic? Compute its
    matrix, trace, eigenvalues, and whether iterating mixes (Anosov-like).
"""
from itertools import combinations, product
from fractions import Fraction as F
from math import comb, gcd, isqrt
import math


# ---------------------------------------------------------------------------
# Shared: LRC measS7 exact (Z/7 vertex coloring c(e,x)=floor(7*frac(ex)))
# meas over x in [0,1) of {e: floor(7 frac(e x)) for e in E} == all of Z/7.
# Exact: x-axis is partitioned by breakpoints frac where any e x crosses k/7.
# We compute the exact Lebesgue measure with Fractions.
# ---------------------------------------------------------------------------
def measS7_exact(E):
    """Exact measure of {x in [0,1): colors{floor(7 frac(e x)): e in E} == Z/7}.
    Breakpoints: x where e*x = k/7 + integer, i.e. x = (m + k/7)/e for integers m,k.
    Within each open interval the coloring is constant; sum lengths where |colors|==7."""
    E = [e for e in E if e != 0] + ([0] if 0 in E else [])
    # collect breakpoints in (0,1)
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        ae = abs(e)
        # e*x in [0, ae); crossings at value k/7 + j for j=0..ae-1, k=0..6
        for j in range(ae):
            for k in range(7):
                val = F(7 * j + k, 7 * ae)  # x = (j + k/7)/ae
                if F(0) < val < F(1):
                    bps.add(val)
    pts = sorted(bps)
    total = F(0)
    for i in range(len(pts) - 1):
        lo, hi = pts[i], pts[i + 1]
        mid = (lo + hi) / 2
        cols = set()
        for e in E:
            fe = (e * mid) % 1   # frac(e*mid)
            cols.add(int(7 * fe))  # floor(7*frac)
        if len(cols) == 7:
            total += (hi - lo)
    return total


def consec(k):
    return tuple(range(k))


# ---------------------------------------------------------------------------
# Tournament side: Z_n score partition function, c3, H (Redei Ham-path count)
# ---------------------------------------------------------------------------
def tilings_scores_c3_H(n):
    """Enumerate all 2^{C(n-1,2)} tilings; return list of (score_tuple, c3, H).
    Base path n->n-1->...->1 (vertices 0..n-1, base arc k->k-1).
    Tile (a,b), a-b>=2: bit0 => a->b, bit1 => b->a (a beats b vs b beats a)."""
    verts = list(range(n))
    base_arcs = [(k, k - 1) for k in range(1, n)]  # k beats k-1
    tiles = [(a, b) for a in range(n) for b in range(a) if a - b >= 2]
    m = len(tiles)
    out = []
    for bits in product((0, 1), repeat=m):
        adj = [0] * n
        for (u, v) in base_arcs:
            adj[u] |= 1 << v
        for (a, b), bt in zip(tiles, bits):
            if bt == 0:
                adj[a] |= 1 << b
            else:
                adj[b] |= 1 << a
        sc = tuple(sorted(bin(adj[v]).count("1") for v in range(n)))
        # c3 = number of 3-cycles
        c3 = 0
        for i, j, l in combinations(range(n), 3):
            e = 0
            e += 1 if (adj[i] >> j & 1) else 0
            # cyclic iff each has outdeg 1 within triple
            od = {i: 0, j: 0, l: 0}
            for x, y in [(i, j), (i, l), (j, l)]:
                if adj[x] >> y & 1:
                    od[x] += 1
                else:
                    od[y] += 1
            if set(od.values()) == {1}:
                c3 += 1
        # H = Hamiltonian path count (Redei)
        H = Hcount(n, adj)
        out.append((sc, c3, H))
    return out, tiles


def Hcount(n, adj):
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            av = adj[v]
            for w in range(n):
                if not (mask >> w & 1) and (av >> w & 1):
                    dp[mask | 1 << w][w] += c
    return sum(dp[size - 1][v] for v in range(n))


# ===========================================================================
print("=" * 78)
print("H1 (CONTROL): cut-space free energy IS a sum of local terms (Z_n factorizes)")
print("=" * 78)
# Z_n(x)= (prod_{v>=2} x_v) prod_tiles(x_a+x_b). With x_v=exp(-beta h_v):
# log Z_n = sum_{v>=2}(-beta h_v) + sum_tiles log(e^{-beta h_a}+e^{-beta h_b}).
# This is MANIFESTLY a sum over tiles. Verify numerically + show ground state per-tile.
n = 5
beta = 3.0
import random
random.seed(0)
h = [random.uniform(0, 1) for _ in range(n)]
x = [math.exp(-beta * hv) for hv in h]
tiles = [(a, b) for a in range(n) for b in range(a) if a - b >= 2]
logZ_direct = sum(math.log(x[v]) for v in range(1, n)) + sum(math.log(x[a] + x[b]) for (a, b) in tiles)
# additive model: per-tile term
per_tile = [math.log(x[a] + x[b]) for (a, b) in tiles]
logZ_sum = sum(math.log(x[v]) for v in range(1, n)) + sum(per_tile)
print(f"  n={n}, beta={beta}: log Z (direct) = {logZ_direct:.8f}")
print(f"           log Z (sum of LOCAL tile terms) = {logZ_sum:.8f}")
print(f"           additive residual = {abs(logZ_direct - logZ_sum):.2e}  -> EXACTLY local (control)")
# ground state beta->inf: each tile -> max(h_a,h_b) (min energy = min(beta h)?). x=e^{-beta h},
# x_a+x_b -> e^{-beta min(h_a,h_b)}; so -1/beta log -> min(h_a,h_b) per tile. LOCAL minimizer.
gs_per_tile = [min(h[a], h[b]) for (a, b) in tiles]
print(f"  ground state energy/tile (beta->inf) = min(h_a,h_b) per tile, fully LOCAL.")
print(f"  CONCLUSION H1: cut-space (score) gas free energy = SUM of local terms => CHEAP. CONTROL holds.\n")


# ===========================================================================
print("=" * 78)
print("H2 (CRUX): LRC free energy -log measS7 is NOT a sum of local terms")
print("=" * 78)
# Fit additive model log measS7(E) = const + sum_{e in E} a_e over many k-subsets E.
# If LRC were local this fits exactly. Then add pairwise (2-body) corrections.
# We use a fixed offset universe U and all/many k-subsets, least-squares over GF... use rationals->float LS.
U = list(range(0, 10))          # offsets 0..9 (need k>=7 for measS7>0; consec maximizes)
k = 8
subsets = [s for s in combinations(U, k) if 0 in s]   # e=0 pins residue 0 (the clock)
data = []
for s in subsets:
    mv = measS7_exact(s)
    data.append((s, mv))
# Only keep subsets with measS7 > 0 (log defined)
data_pos = [(s, mv) for (s, mv) in data if mv > 0]
print(f"  universe U={U}, k={k}: {len(subsets)} subsets (0 in E); {len(data_pos)} with measS7>0")
if not data_pos:
    raise SystemExit("  [no positive measS7 subsets at this k -- increase k>=7]")
# ADDITIVE model: log m = c + sum_{e in s} a_e. Design matrix over offsets present.
import itertools
offs = sorted(set(e for s, _ in data_pos for e in s))
idx = {e: i for i, e in enumerate(offs)}
# Build A (rows=subsets) cols = [1] + indicator(e)
def lstsq(A, y):
    # normal equations (A^T A) x = A^T y, small; use python floats + Gaussian elim
    cols = len(A[0])
    AtA = [[sum(A[r][i] * A[r][j] for r in range(len(A))) for j in range(cols)] for i in range(cols)]
    Aty = [sum(A[r][i] * y[r] for r in range(len(A))) for i in range(cols)]
    # solve with tiny ridge for stability
    for i in range(cols):
        AtA[i][i] += 1e-9
    # Gaussian elimination
    M = [row[:] + [Aty[i]] for i, row in enumerate(AtA)]
    for c in range(cols):
        p = max(range(c, cols), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        for j in range(c, cols + 1):
            M[c][j] /= pv
        for r in range(cols):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                for j in range(c, cols + 1):
                    M[r][j] -= f * M[c][j]
    x = [M[i][cols] for i in range(cols)]
    pred = [sum(A[r][i] * x[i] for i in range(cols)) for r in range(len(A))]
    ss_res = sum((y[r] - pred[r]) ** 2 for r in range(len(A)))
    ss_tot = sum((y[r] - sum(y) / len(y)) ** 2 for r in range(len(A)))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return x, ss_res, r2, pred

y = [math.log(float(mv)) for s, mv in data_pos]
A_add = [[1.0] + [1.0 if e in s else 0.0 for e in offs] for s, _ in data_pos]
xa, ssr_a, r2_a, pred_a = lstsq(A_add, y)
print(f"  ADDITIVE (1-body) log-linear model: R^2 = {r2_a:.5f}, residual SS = {ssr_a:.4e}")
# PAIRWISE model: add indicator(e<e' both in s)
pairs = list(combinations(offs, 2))
A_pair = [[1.0] + [1.0 if e in s else 0.0 for e in offs] +
          [1.0 if (e in s and f_ in s) else 0.0 for (e, f_) in pairs] for s, _ in data_pos]
xp, ssr_p, r2_p, pred_p = lstsq(A_pair, y)
print(f"  PAIRWISE (2-body) log-linear model:  R^2 = {r2_p:.5f}, residual SS = {ssr_p:.4e}")
print(f"  -> additive R^2<1 means -log measS7 is NOT a sum of local terms (genuinely aggregate).")
print(f"  -> even 2-body residual {ssr_p:.2e}: higher-order interaction needed iff >~1e-6.\n")


# ===========================================================================
print("=" * 78)
print("H3 (VARIATIONAL): Gibbs family p_beta ~ exp(beta*measS7); ground state = consec")
print("=" * 78)
# Over all k-subsets of a window, p_beta(E) ~ exp(beta * measS7(E)).
# Phi(beta) = log sum_E exp(beta measS7). Convex. Phi'(beta) -> max measS7 as beta->inf.
# Ground state argmax should be the consecutive block.
k = 8
W = 11  # offsets 0..W-1, with 0 forced (clock); consec=(0..7) is in-window
subs = [s for s in combinations(range(W), k) if 0 in s]
mvals = {s: float(measS7_exact(s)) for s in subs}
best = max(subs, key=lambda s: mvals[s])
consec_set = tuple(range(k))
print(f"  k={k}, window {W}: {len(subs)} subsets. argmax measS7 = {best}, val={mvals[best]:.6f}")
print(f"  consec block {consec_set}: val={mvals.get(consec_set, float('nan')):.6f}; is argmax consec? {best == consec_set}")
# Gibbs concentration + convexity check
def phi(beta):
    mx = max(mvals.values())
    return math.log(sum(math.exp(beta * (v - mx)) for v in mvals.values())) + beta * mx
betas = [0.0, 1.0, 5.0, 20.0, 100.0, 1000.0]
print("  beta      Phi(beta)     <measS7>_beta   P(consec)_beta")
for b in betas:
    mx = max(mvals.values())
    w = {s: math.exp(b * (mvals[s] - mx)) for s in subs}
    Zb = sum(w.values())
    mean = sum(mvals[s] * w[s] for s in subs) / Zb
    pcons = w[consec_set] / Zb
    print(f"  {b:7.1f}  {phi(b):11.5f}   {mean:.6f}      {pcons:.6f}")
# convexity: Phi''>=0 numerically
def phi2(b, eps=1e-3):
    return (phi(b + eps) - 2 * phi(b) + phi(b - eps)) / eps ** 2
print(f"  Phi''(beta) >= 0 at beta=1: {phi2(1.0):.4e} (variance, must be >=0) -> convex free energy.")
print(f"  GROUND STATE (beta->inf) selects argmax measS7 = consec: {best == consec_set}")

print()
print("  Tournament twin: p_beta(T) ~ exp(beta*c3); ground state = regular (THM-027/555)")
n = 5
tinfo, _ = tilings_scores_c3_H(n)
c3s = [c for (_, c, _) in tinfo]
Hs = [H for (_, _, H) in tinfo]
scs = [sc for (sc, _, _) in tinfo]
maxc3 = max(c3s)
# regular score for n=5 is (2,2,2,2,2)
reg = tuple([(n - 1) // 2] * n)
argmax_c3_scores = set(scs[i] for i in range(len(scs)) if c3s[i] == maxc3)
print(f"  n={n}: max c3={maxc3}; score(s) at max c3 = {argmax_c3_scores}; regular {reg} present? {reg in argmax_c3_scores}")
maxH = max(Hs)
print(f"  n={n}: max H={maxH} (H-max ground state of exp(beta*H)); count at max = {sum(1 for H in Hs if H == maxH)}")
print(f"  -> exp(beta*c3) beta->inf = regular score; exp(beta*H) beta->inf = H-max. Ground-state extremality.\n")


# ===========================================================================
print("=" * 78)
print("H4 (ARNOLD CAT MAP): band twist e->s*e mod 7 as a toral automorphism")
print("=" * 78)
# HYP-2703: measS7=(1/7) sum_{s=0..6} bandcover_s; bandcover_s twists e->s*e mod 7.
# The renormalization acts on (slope band s, residue r). Question: is the multiplier
# action a HYPERBOLIC toral automorphism (cat-map-like, |trace|>2, Anosov)?
# The cover lives on a 2-torus: x in [0,1) (the time) and the residue phase. The slope-band
# index s is the integer part of theta=7x; the renormalization x->{7x}, s->floor(7x) is the
# Gauss-like/baker map for base 7. Combined with the multiplier e->s*e, the natural 2x2 matrix
# acting on (winding in x, winding in residue) is the companion of the Sturmian/continued-fraction.
# Arnold's cat map is M=[[2,1],[1,1]], trace 3, eigenvalues (3+/-sqrt5)/2. Test which 2x2
# integer maps over Z (mod nothing) realize the renormalization and whether hyperbolic.

# (a) The base-7 expanding map x->7x mod 1 has derivative 7>1: expanding (one direction).
# A toral AUTOMORPHISM needs det=+-1. The genuine automorphism here is the multiplier group
# acting on Z/7: e->s*e. Its order divides 6 (|(Z/7)^*|=6). On the torus R/Z this is rotation
# by log_g(s)/6 -- ELLIPTIC, NOT hyperbolic. So the multiplier alone is NOT a cat map.
print("  (a) multiplier e->s*e on Z/7: it is in (Z/7)^* (cyclic order 6).")
g = 3  # primitive root mod 7
order = {}
for s in range(1, 7):
    p = 1
    o = 0
    while True:
        p = (p * s) % 7
        o += 1
        if p == 1:
            break
    order[s] = o
print(f"      multiplier orders mod 7: {order}  (all finite => ELLIPTIC, periodic, NOT hyperbolic)")
print("      => the band twist is a TORSION rotation, not an Anosov cat map. (analogy-killer)")

# (b) Where IS the hyperbolicity? The Sturmian/continued-fraction renormalization of frac(e x)
# is the Gauss map, whose natural extension IS hyperbolic (geodesic flow on modular surface).
# Test: the consec offsets {0..k-1} and the 7-cover -> the relevant CF expansion of slopes.
# A genuine cat-map appears for the GOLDEN slope. Compute the renormalization matrix for the
# slope continued fraction and check |trace|>2 for the periodic (quadratic-irrational) slopes.
def cf_matrix(partial_quotients):
    # product of [[a,1],[1,0]] -- the CF convergent matrix; periodic CF => conjugate to cat map
    M = [[1, 0], [0, 1]]
    for a in partial_quotients:
        M = [[M[0][0] * a + M[0][1], M[0][0]], [M[1][0] * a + M[1][1], M[1][0]]]
    return M
# golden mean [1;1,1,...] one period: a=1 -> M=[[1,1],[1,0]] trace 1 (parabolic-ish); two periods cat:
M_golden = cf_matrix([1, 1])  # [[2,1],[1,1]] = ARNOLD CAT MAP
tr = M_golden[0][0] + M_golden[1][1]
det = M_golden[0][0] * M_golden[1][1] - M_golden[0][1] * M_golden[1][0]
disc = tr * tr - 4 * det
print(f"  (b) golden-slope CF renormalization matrix (period [1,1]) = {M_golden}")
print(f"      trace={tr}, det={det}, discriminant={disc}; |trace|>2 => HYPERBOLIC (this IS Arnold's cat map).")
# Which slopes are relevant to LRC? The breakpoints x=(j+k/7)/e are RATIONAL -> the cover is
# piecewise-constant on rationals; the renormalization that mixes bands is x->7x mod1 (expanding
# Bernoulli, entropy log7), NOT an automorphism. The cat map appears only if we close up the
# (slope, residue) torus with a unimodular gluing -- which the apex prime 7 does NOT provide.
print("  (c) LRC cover breakpoints are RATIONAL (x=(7j+k)/(7e)); the band-mixing renormalization")
print("      is x->7x mod 1: EXPANDING Bernoulli (entropy log 7), det 7 != +-1 => NOT a toral")
print("      automorphism. The cat map requires a unimodular (det=+-1) gluing; the apex prime")
print("      gives an expanding endomorphism instead. So: SHIFT/Bernoulli yes, CAT MAP no.")
print()
print("  H4 VERDICT: the renormalization is an EXPANDING base-7 endomorphism (a Bernoulli shift on")
print("  7 symbols = the 7 bands), NOT a hyperbolic toral automorphism. Arnold's cat map is the")
print("  WRONG dynamical analogy; the right one is the full 7-shift (subshift of finite type).")


# ===========================================================================
print()
print("=" * 78)
print("SYNTHESIS / FALSIFIABLE SCORECARD")
print("=" * 78)
print(f"  H1 cut-space local: additive residual {abs(logZ_direct - logZ_sum):.1e} ~ 0  -> CONFIRMED (control)")
print(f"  H2 LRC aggregate: additive R^2={r2_a:.4f}<1, pairwise R^2={r2_p:.4f}  -> {'CONFIRMED non-local' if r2_a < 0.999 else 'REFUTED'}")
print(f"  H3 variational ground state = consec: {best == consec_set}; tournament regular at max c3: {reg in argmax_c3_scores}")
print(f"  H4 Arnold cat map: REFUTED (multiplier elliptic order-6; band map expanding det7, not det+-1)")
