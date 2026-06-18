#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_tourmap_character-spectral_kps-S2-wf.py        kind-pasteur S2-wf

CREATIVE TOURNAMENT-GENERATION from LRC(14), theme = "CHARACTER-SPECTRAL".

We map convoluted spectral/character aspects of the Lonely Runner setup to a
tournament on small vertex sets, GATE each map for NON-triviality (does it ever
make a 3-cycle / H>1?), and then EXHAUSTIVELY enumerate which tournament ISO
CLASSES are realized as (S, tau) range over LRC-constrained inputs.  A class in
the FREE set (all A000568 iso classes on n vertices) that is NEVER realized =
a FORBIDDEN CLASS.

All decisions use EXACT rationals (fractions.Fraction) on Z[i] / cyclotomic
sums; nothing depends on a float comparison.

THEME = character sums e(v tau) = exp(2 pi i v tau).  At a rational tau = p/q the
value e(v tau) is a q-th root of unity, EXACTLY representable.  We build several
SIGNED spectral quantities and turn their signs / orderings / parities into arcs.

------------------------------------------------------------------------------
THE SIX METHODS (each: vertices, arc rule, source of asymmetry)
------------------------------------------------------------------------------
M1  CORRELATION-INTEGRAL tournament (pairs of runners).
    Vertices = runners.  arc i->j by the SIGN of the exact cross-correlation
    integral  C_{ij} = Re integral_0^1 e(v_i t) conj(e(v_j t)) w(t) dt  with a
    NON-symmetric danger weight w(t) (indicator of the danger band of a FIXED
    reference runner).  This is a genuine signed pairwise spectral statistic.
    [predicted: likely degenerate because integral of e((v_i-v_j)t) over a band
    is governed by |v_i-v_j| -> we TEST.]

M2  WEYL-PHASE / ARGUMENT tournament (runners at the lonely tau*).
    Vertices = runners.  At tau*, runner i has phase theta_i = frac(v_i tau*) in
    [0,1).  arc i->j by sign of Im( e(theta_i) conj(e(theta_j)) ) = sign of
    sin(2pi(theta_i - theta_j)) = "is j within half a turn CCW ahead of i".
    This is the OVERTAKING snapshot in disguise -> predicted TRANSITIVE.  GATE
    and discard if trivial.

M3  RADEMACHER / WALSH-SIGN tournament (runners, parity of danger-crossings).
    Vertices = runners.  For runner i define the Rademacher-style sign vector
    over the grid times a in (Z/14)*: s_i(a) = +1 if runner i is in a "near"
    section at a/14 (section in {1,13}, i.e. within one section of 0 -> the
    danger band on the grid), else -1.  arc i->j by sign of the Walsh
    correlation  <s_i, s_j> = sum_a s_i(a) s_j(a); tie-break by speed.  This is a
    Walsh/Rademacher signed-correlation tournament built from the SECTION data.

M4  SECTION-CHARACTER (Gauss-sum) tournament on SECTIONS (vertices = Z/14
    sections, 14 vertices, but we look at the induced sub-tournament on the
    OCCUPIED sections).  For two sections r,s use the additive character
    chi(r,s) = sum_{a in units} e( a (r - s) / 14 ) which is a real Gauss-type
    sum (Ramanujan sum c_14(r-s)).  arc r->s by sign of an ANTISYMMETRIZED
    version: sign of sum_{a} a-weighted e(a(r-s)/14)   (the weight a breaks the
    r<->s symmetry).  Genuine character-sum source.

M5  BINDING-PARITY tournament (runners; uses THM-524 binding structure).
    Vertices = runners.  At tau* each runner i has a "crossing index"
    k_i defined by frac(v_i tau*) and a winding number w_i = floor(v_i tau*).
    arc i->j by the parity-twisted comparison
        (w_i + w_j) even ?  (theta_i > theta_j) : (theta_i < theta_j),
    i.e. the order is REVERSED whenever the total winding parity is odd.  Parity
    twists CAN break transitivity -> TEST.

M6  RESONANCE / DIVISOR tournament (runners; pure arithmetic of speeds, no tau).
    Vertices = runners.  arc i->j by sign of the signed resonance functional
        R_{ij} = c_14(v_i)*v_j - c_14(v_j)*v_i        (Ramanujan sum c_14)
    antisymmetric by construction; ties -> speed order.  This injects the
    arithmetic CHARACTER c_14(v) (a Ramanujan/Gauss sum) of each speed; not a
    snapshot order, so it can be non-transitive.

For every NON-trivial map we range over LRC-constrained input families at small
vertex count n in {3,4,5} (sometimes 6), canonicalize each tournament to an iso
class (brute over all n! relabelings), and bucket.  We compare the realized set
to the FULL free set A000568(n).

LRC-CONSTRAINED INPUT FAMILIES (so the enumeration is honest):
  - "all n-subsets of small speeds" with each evaluated at its OWN exact lonely
    tau* (so tau is genuinely the LRC optimum, not arbitrary);
  - "primitive n-subsets of {1..K}" for K up to a cap;
  - for tau-free maps (M6, M3, M4) we range over residue/speed patterns directly.

stdlib only.  Exact rationals + Z[i]/cyclotomic via Fraction pairs.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
from collections import Counter, defaultdict

# ======================================================================
# exact M tool (verbatim, validated)
# ======================================================================
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

# ======================================================================
# tournament canonicalization (brute, n<=6) + iso-class bookkeeping
# ======================================================================
def is_tournament(adj, V):
    """adj is dict (a,b)->bool meaning a->b. check antisymmetry+completeness."""
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            a, b = V[i], V[j]
            if adj.get((a, b), None) == adj.get((b, a), None):
                return False
    return True

def adj_matrix(adj, V):
    """0/1 matrix mat[i][j]=1 iff V[i]->V[j]."""
    n = len(V)
    mat = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and adj[(V[i], V[j])]:
                mat[i][j] = 1
    return mat

def canon(mat):
    """canonical form under vertex relabeling: min over all permutations of the
    flattened upper-coding. n<=6 so 6!=720 fine."""
    n = len(mat)
    best = None
    for p in permutations(range(n)):
        # relabeled matrix entry (p[i],p[j]) -> i,j position; build code
        code = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    code.append(mat[p[i]][p[j]])
        code = tuple(code)
        if best is None or code < best:
            best = code
    return best

def ham_count(mat):
    n = len(mat)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for k in range(n): dp[1 << k][k] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if mat[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full][k] for k in range(n))

def num_3cycles(mat):
    n = len(mat); c = 0
    for a in range(n):
        for b in range(n):
            for cc in range(n):
                if a < b and b < cc:
                    # count directed 3-cycles on {a,b,cc} (0 or 1)
                    if mat[a][b] and mat[b][cc] and mat[cc][a]:
                        c += 1
                    if mat[a][cc] and mat[cc][b] and mat[b][a]:
                        c += 1
    return c

def score_seq(mat):
    return tuple(sorted(sum(row) for row in mat))

# enumerate the FULL free set: all iso classes of tournaments on n vertices
def free_classes(n):
    """return dict canon-code -> representative invariants (score,#3cyc,H)."""
    classes = {}
    # iterate over all tournaments: choose orientation of each of C(n,2) edges
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(pairs)):
        mat = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                mat[i][j] = 1
            else:
                mat[j][i] = 1
        c = canon(mat)
        if c not in classes:
            classes[c] = (score_seq(mat), num_3cycles(mat), ham_count(mat))
    return classes

# precompute free sets
print("Precomputing free iso-class sets (A000568)...")
FREE = {}
for n in (3, 4, 5):
    FREE[n] = free_classes(n)
    print(f"  n={n}: {len(FREE[n])} free iso classes")
# n=6 has 56 classes; building all 2^15 = 32768 tournaments * 720 perms is heavy
# but feasible; do it lazily only if a method needs it.
print()

def class_label(mat):
    c = canon(mat)
    sc = score_seq(mat); t3 = num_3cycles(mat); h = ham_count(mat)
    return c, (sc, t3, h)

# ======================================================================
# exact cyclotomic helper: e(p/q) as (cos,sin) is irrational, but we only ever
# need SIGNS of real linear combos of cos(2pi a/q) and sin(2pi a/q).  To keep
# everything EXACT we represent the needed character sums symbolically when they
# reduce to rationals (Ramanujan/Gauss sums are INTEGERS), and otherwise use a
# high-precision RATIONAL approximation of cos/sin via the Fraction of a Taylor /
# Chebyshev evaluation -- BUT to stay honest we only base ARC DECISIONS on
# quantities that are provably exact integers (Ramanujan sums) or on rational
# comparisons of fractional parts.  Where a genuine transcendental sign is
# needed (M1 correlation integral, M2 phase), the integral evaluates to an
# EXACT rational because integral of cos(2pi k t) over rational endpoints folds
# into sin at rational multiples -> we handle those cases exactly below.
# ======================================================================

def ramanujan_c(q, n):
    """Ramanujan sum c_q(n) = sum_{a in (Z/q)*} e(a n / q).  Integer-valued.
    Compute via the standard multiplicative formula using gcd:
       c_q(n) = mu(q/d) * phi(q) / phi(q/d),  d = gcd(n,q)."""
    d = gcd(n % q, q) if (n % q) != 0 else q
    qd = q // d
    return _mobius(qd) * _phi(q) // _phi(qd)

def _phi(m):
    if m == 1: return 1
    result = m; p = 2; mm = m
    while p*p <= mm:
        if mm % p == 0:
            while mm % p == 0: mm //= p
            result -= result // p
        p += 1
    if mm > 1: result -= result // mm
    return result

def _mobius(m):
    if m == 1: return 1
    res = 1; mm = m; p = 2
    while p*p <= mm:
        if mm % p == 0:
            mm //= p
            if mm % p == 0: return 0
            res = -res
        p += 1
    if mm > 1: res = -res
    return res

UNITS14 = [a for a in range(1, 14) if gcd(a, 14) == 1]   # {1,3,5,9,11,13}

# ======================================================================
# REALIZATION ENUMERATOR
# ======================================================================
def enumerate_realized(arc_fn, families, n, name, note=""):
    """arc_fn(S, tau) -> (adj dict, V list) or None if config invalid.
    families: iterable of S (speed tuples) of length n.  For each S we evaluate
    at its OWN lonely tau* (arc_fn may ignore tau).  Bucket iso classes.
    Returns (realized_dict canon->invariants, n_inputs, n_nontrivial, examples)."""
    realized = {}
    n_inputs = 0; n_nontrivial = 0
    examples = {}  # canon -> sample S
    for S in families:
        Mv, tau = M(S)
        out = arc_fn(S, tau)
        if out is None:
            continue
        adj, V = out
        if not is_tournament(adj, V):
            continue
        n_inputs += 1
        mat = adj_matrix(adj, V)
        c = canon(mat)
        h = ham_count(mat)
        if h > 1:
            n_nontrivial += 1
        if c not in realized:
            realized[c] = (score_seq(mat), num_3cycles(mat), h)
            examples[c] = S
    return realized, n_inputs, n_nontrivial, examples

def report_method(name, arc_fn, families_by_n, note=""):
    print("#" * 78)
    print(f"METHOD {name}")
    if note: print("  " + note)
    print("#" * 78)
    any_nontrivial = False
    forbidden_total = {}
    for n, fams in families_by_n.items():
        fams = list(fams)
        realized, ninp, nnt, ex = enumerate_realized(arc_fn, fams, n, name)
        free_n = FREE.get(n)
        if free_n is None:
            free_n = free_classes(n); FREE[n] = free_n
        realized_set = set(realized.keys())
        free_set = set(free_n.keys())
        forbidden = free_set - realized_set
        if nnt > 0: any_nontrivial = True
        print(f"\n  n={n}: {ninp} valid tournaments from {len(fams)} input configs; "
              f"{nnt} NON-transitive (H>1).")
        print(f"    realized {len(realized_set)}/{len(free_set)} free iso classes.")
        # show realized invariants
        inv_list = sorted((free_n[c][2], free_n[c][1], free_n[c][0]) for c in realized_set)
        print(f"    realized (H,#3cyc,score): {inv_list}")
        if forbidden:
            forb_inv = sorted((free_n[c][2], free_n[c][1], free_n[c][0]) for c in forbidden)
            print(f"    FORBIDDEN {len(forbidden)} classes (H,#3cyc,score): {forb_inv}")
            forbidden_total[n] = forb_inv
        else:
            print(f"    FORBIDDEN: none (realizes everything).")
    print(f"\n  --> METHOD {name}: nontrivial? {any_nontrivial}\n")
    return any_nontrivial, forbidden_total

# ======================================================================
# INPUT FAMILIES (LRC-constrained)
# ======================================================================
def primitive_subsets(K, n, cap=None):
    """all n-subsets of {1..K} that are PRIMITIVE (gcd=1), as sorted tuples."""
    out = []
    for combo in combinations(range(1, K+1), n):
        if gcd_list(combo) == 1:
            out.append(combo)
            if cap and len(out) >= cap:
                break
    return out

def gcd_list(xs):
    g = 0
    for x in xs: g = gcd(g, x)
    return g

# families per vertex count (kept exhaustive but small)
def make_families():
    fam = {}
    # n=3: all primitive 3-subsets of {1..18}
    fam[3] = primitive_subsets(18, 3)
    # n=4: all primitive 4-subsets of {1..15}
    fam[4] = primitive_subsets(15, 4)
    # n=5: all primitive 5-subsets of {1..13}  (manageable)
    fam[5] = primitive_subsets(13, 5)
    return fam

FAM = make_families()
print("Input family sizes (primitive subsets):")
for n in FAM: print(f"  n={n}: {len(FAM[n])} configs")
print()

# ======================================================================
# M1  CORRELATION-INTEGRAL tournament  (EXACT)
# integral_0^1 e(v_i t) conj(e(v_j t)) w(t) dt with w = indicator of danger band
# of a FIXED reference (the slowest runner v0): band B = { t : ||v0 t|| < 1/14 }.
# Re part = integral_B cos(2pi (v_i - v_j) t) dt.  Over a union of rational
# intervals this is an EXACT rational combination of sin at rational multiples of
# 2pi... which is NOT rational in general.  To stay EXACT we instead use the
# DISCRETE correlation over the grid sample t in {a/14 : a in units}, weighted by
# the danger weight, with the REAL part taken as the exact integer
# cos-sum surrogate sum_a w(a) * Re[ e((v_i-v_j) a /14) ].  Re[e(m a/14)] sums
# over a to a Ramanujan-type INTEGER.  We antisymmetrize by an explicit weight.
# ======================================================================
def m1_arc(S, tau):
    V = list(S)
    v0 = min(S)
    # danger weight on the grid: w(a) = 1 if v0*a mod 14 in {1,2,12,13} (near 0), else 0
    # (near band ~ within 2 sections of 0)
    def w(a):
        r = (v0 * a) % 14
        return 1 if r in (1, 2, 12, 13) else 0
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            m = (x - y)
            # exact integer surrogate: S_xy = sum_a w(a)*cos(2pi m a/14) ... cos
            # is not integer per-term, but SUM over all a in units of cos(2pi m a/14)
            # = ramanujan_c(14, m).  With weight w(a) we instead use a RATIONAL
            # exact value: we know cos(2pi m a /14) for a in units only takes a few
            # algebraic values; to keep EXACT decisions we use the SIGNED count
            # surrogate: D_xy = sum_a w(a)*( (v0*a%14 distance pattern) ) ...
            # SIMPLEST EXACT ASYMMETRIC functional:
            #   D_xy = sum_a w(a) * sign( nrm(x*F(a,14)) - nrm(y*F(a,14)) )
            # this is the danger-weighted count of grid times where x is FARTHER
            # from 0 than y.  Antisymmetric: D_xy = -D_yx.  Exact (rational nrm).
            D = 0
            for a in UNITS14:
                if w(a):
                    dx = nrm(x * F(a, 14)); dy = nrm(y * F(a, 14))
                    if dx > dy: D += 1
                    elif dx < dy: D -= 1
            if D > 0: adj[(x, y)] = True
            elif D < 0: adj[(x, y)] = False
            else:
                adj[(x, y)] = x > y  # tie-break by speed
    return adj, V

# ======================================================================
# M2  WEYL-PHASE tournament: i->j iff sin(2pi(theta_i-theta_j))>0,
# theta=frac(v tau*).  sign of sin(2pi d) for d in (0,1): >0 iff d<1/2.
# EXACT: arc x->y iff frac((x-y) tau*) in (0,1/2).  This is the overtaking
# snapshot -> predicted transitive.
# ======================================================================
def m2_arc(S, tau):
    V = list(S)
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            d = ((x - y) * tau) % 1
            if d == 0:
                adj[(x, y)] = x > y
            elif d < F(1, 2):
                adj[(x, y)] = True
            elif d > F(1, 2):
                adj[(x, y)] = False
            else:  # d == 1/2 exactly
                adj[(x, y)] = x > y
    return adj, V

# ======================================================================
# NOTE (design correction):  a SYMMETRIC Walsh correlation <s_x,s_y> = <s_y,s_x>
# CANNOT orient arcs (it gives the same sign for (x,y) and (y,x)) -> never a
# tournament (verified: 0 valid configs).  An arc statistic MUST be ANTISYMMETRIC
# under x<->y.  We therefore build the M3 family from ANTISYMMETRIC Walsh/Rademacher
# functionals: a danger-WEIGHTED majority vote of "who is farther from 0" across
# the grid, where the weight is a character/Rademacher sign.  D_xy = -D_yx by
# construction.
# ======================================================================
# M3  WALSH-WEIGHTED LEAD tournament.  Reference = slowest runner v0.  Weight each
# grid time a in units by the Rademacher sign w(a) of whether v0 is in the NEAR
# band {1,13} (so danger times count +, safe times count -).  D_xy = sum_a w(a)*
# sign( ||x a/14|| - ||y a/14|| ).  arc x->y iff D_xy>0.  Antisymmetric.
def m3_arc(S, tau):
    V = list(S); v0 = min(S)
    def w(a):
        r = (v0 * a) % 14
        return 1 if r in (1, 13) else -1
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            D = 0
            for a in UNITS14:
                dx = nrm(x * F(a, 14)); dy = nrm(y * F(a, 14))
                if dx > dy: D += w(a)
                elif dx < dy: D -= w(a)
            if D > 0: adj[(x, y)] = True
            elif D < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# M3b: ANTISYMMETRIC section-parity character vote.  weight w(a)=(-1)^{section of
# v0 at a} (Rademacher of the parity character), compare section-NUMBERS of x,y.
def m3b_arc(S, tau):
    V = list(S); v0 = min(S)
    def w(a):
        r = (v0 * a) % 14
        return 1 if r % 2 == 0 else -1
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            D = 0
            for a in UNITS14:
                sx = (x * a) % 14; sy = (y * a) % 14
                # signed distance of section from 0 on the cycle of 14
                cx = min(sx, 14 - sx); cy = min(sy, 14 - sy)
                if cx > cy: D += w(a)
                elif cx < cy: D -= w(a)
            if D > 0: adj[(x, y)] = True
            elif D < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# M3c: ANTISYMMETRIC Rademacher half-circle vote.  weight w(a)=sign(frac(v0 a/14)
# - 1/2); compare half-circle Rademacher of x,y across grid.  D_xy = sum_a w(a)*
# (r_x(a) - r_y(a)) where r(a)=sign(frac(v a/14)-1/2).
def m3c_arc(S, tau):
    V = list(S); v0 = min(S)
    def rad(v, a):
        fr = (v * F(a, 14)) % 1
        return 1 if fr > F(1, 2) else -1
    def w(a): return rad(v0, a)
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            D = 0
            for a in UNITS14:
                D += w(a) * (rad(x, a) - rad(y, a))
            if D > 0: adj[(x, y)] = True
            elif D < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# ======================================================================
# M5  BINDING-PARITY tournament.  At tau*: theta_i=frac(v_i tau*),
# w_i=floor(v_i tau*) winding.  arc i->j: if (w_i+w_j) even -> (theta_i>theta_j)
# else reversed.
# ======================================================================
def m5_arc(S, tau):
    V = list(S)
    th = {v: (v * tau) % 1 for v in V}
    wd = {v: int(v * tau) for v in V}   # floor, tau in [0,1) so this is 0 mostly;
    # to make winding meaningful use the crossing count to tau* of the runner vs
    # observer = floor(v*tau*)+ (#half-laps).  We use a richer winding: number of
    # times runner i has passed 0 strictly before tau* = floor(v_i * tau*).
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            base = th[x] > th[y]
            if th[x] == th[y]:
                adj[(x, y)] = x > y; continue
            par = (wd[x] + wd[y]) % 2
            if par == 1:
                base = not base
            adj[(x, y)] = base
    return adj, V

# M5b: parity twist by the CROSSING count c_xy = floor(|x-y| tau*) (overtake count)
def m5b_arc(S, tau):
    V = list(S)
    th = {v: (v * tau) % 1 for v in V}
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            base = th[x] > th[y]
            if th[x] == th[y]:
                adj[(x, y)] = x > y; continue
            cc = int(abs(x - y) * tau)   # overtaking count on (0,tau*]
            if cc % 2 == 1:
                base = not base
            adj[(x, y)] = base
    return adj, V

# ======================================================================
# M6  RESONANCE/DIVISOR tournament (tau-free, pure arithmetic character).
# arc i->j by sign of R_{ij} = c_14(v_i)*v_j - c_14(v_j)*v_i, antisymmetric.
# c_14 = Ramanujan sum (integer character of the speed mod 14).
# ======================================================================
def m6_arc(S, tau):
    V = list(S)
    c = {v: ramanujan_c(14, v) for v in V}
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            R = c[x]*y - c[y]*x
            if R > 0: adj[(x, y)] = True
            elif R < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# M6b: resonance with a DIFFERENT character pairing -- use c_14(v_i - v_j) and
# c_14(v_i + v_j) (the BINDING moduli from THM-524!) to decide:
#   arc i->j iff c_14(v_i+v_j) + c_14(... )  -- need antisymmetry, so use
#   R = c_14(v_i+v_j)*(v_i - v_j) ... but c_14(v_i+v_j) is symmetric in i,j, so
#   multiply by the antisymmetric (v_i - v_j).  arc i->j iff
#   c_14(v_i+v_j)*(v_i-v_j) > 0  (sign flips with i<->j -> valid tournament).
def m6b_arc(S, tau):
    V = list(S)
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            cs = ramanujan_c(14, x + y)   # binding-sum character (symmetric)
            R = cs * (x - y)
            if R > 0: adj[(x, y)] = True
            elif R < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# M6c: binding-DIFFERENCE character: arc i->j iff c_14(v_i - v_j) twist.
# c_14(v_i - v_j) = c_14(v_j - v_i) (even), so again multiply by sign(x-y):
#   R = c_14(|x-y|) * (x - y).  But this depends only on |x-y| sign -> basically
#   a speed order unless c_14 changes sign with the gap.  TEST.
def m6c_arc(S, tau):
    V = list(S)
    adj = {}
    for x in V:
        for y in V:
            if x == y: continue
            cd = ramanujan_c(14, abs(x - y))
            R = cd * (1 if x > y else -1)
            if R > 0: adj[(x, y)] = True
            elif R < 0: adj[(x, y)] = False
            else: adj[(x, y)] = x > y
    return adj, V

# ======================================================================
# RUN ALL METHODS
# ======================================================================
results = {}

def run(name, fn, note=""):
    nontrivial, forb = report_method(name, fn, FAM, note)
    results[name] = (nontrivial, forb)

run("M1 correlation-integral (danger-weighted)", m1_arc,
    "arc x->y by sign of danger-band-weighted grid count of {x farther from 0 than y}.")
run("M2 Weyl-phase snapshot", m2_arc,
    "arc x->y iff frac((x-y)tau*) in (0,1/2) -- the overtaking snapshot (predicted transitive).")
run("M3 Rademacher/Walsh sign (near-section {1,13})", m3_arc,
    "arc by sign of Walsh correlation of near-band indicators across grid times.")
run("M3b Walsh sign (section parity)", m3b_arc,
    "arc by Walsh correlation of (-1)^section across grid times.")
run("M3c Rademacher half-circle sign", m3c_arc,
    "arc by Walsh correlation of sign(frac - 1/2) across grid times.")
run("M5 binding-parity (winding twist)", m5_arc,
    "arc x->y = (theta_x>theta_y) XOR (winding parity).")
run("M5b overtaking-parity (crossing twist)", m5b_arc,
    "arc x->y = (theta_x>theta_y) XOR (overtake-count parity).")
run("M6 resonance c_14(v) functional", m6_arc,
    "arc by sign of c_14(v_i)v_j - c_14(v_j)v_i (Ramanujan character).")
run("M6b binding-sum character", m6b_arc,
    "arc by sign of c_14(v_i+v_j)*(v_i-v_j).")
run("M6c binding-diff character", m6c_arc,
    "arc by sign of c_14(|v_i-v_j|)*sign(v_i-v_j).")

# ======================================================================
# SUMMARY
# ======================================================================
print("=" * 78)
print("SUMMARY: nontriviality + forbidden classes")
print("=" * 78)
for name, (nt, forb) in results.items():
    flag = "NON-TRIVIAL" if nt else "trivial (transitive only) -- DEAD"
    print(f"\n{name}: {flag}")
    if nt:
        for n, fl in forb.items():
            if fl:
                print(f"    n={n}: forbidden {len(fl)} classes -> {fl}")
            else:
                print(f"    n={n}: realizes ALL free classes (no forbidden)")
print("\nDONE.")
