"""
lrc14_tourmap_diophantine_kps-S2-wf.py

CREATIVE TOURNAMENT-GENERATION from LRC data, THEME = DIOPHANTINE APPROXIMATION.

We build several DISTINCT maps  (LRC speed set S, lonely time tau) -> tournament T
on the runners (or on pairs / sections), where each arc i->j is decided by a
Diophantine comparison (continued fractions, three-distance gaps, Stern-Brocot
ancestry, approximation quality ||v_i tau|| etc.) that need NOT be a total order.

For each method:
 (a) define the arc rule exactly,
 (b) NON-TRIVIALITY GATE: does it ever make a 3-cycle (H>1)?  If always transitive -> dead.
 (c) for non-trivial ones, enumerate realized iso classes over LRC-constrained inputs
     at small n (3,4,5) and compare to the FREE set (A000568 = 2,4,12).
     A FREE class never realized = FORBIDDEN.

All decisions use exact Fractions. Floats are NEVER used for a decision.
"""

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

# ----------------------------------------------------------------------
# EXACT M TOOL (validated, verbatim from task spec)
# ----------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def signed_frac(x):
    # signed distance to nearest integer in (-1/2, 1/2], used by some methods
    r = x - int(x); r = r + 1 if r < 0 else r          # frac in [0,1)
    return r if r <= F(1, 2) else r - 1                 # in (-1/2,1/2]

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2))
    return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

# ----------------------------------------------------------------------
# TOURNAMENT CANONICALIZATION
# A tournament on vertices 0..n-1 is a dict adj[(i,j)] in {0,1}: 1 means i->j.
# We store as a frozenset of arcs (i,j) meaning i beats j.
# Canonical form = lexicographically minimal adjacency over all n! relabelings.
# ----------------------------------------------------------------------
def arcs_to_matrix(arcs, n):
    Mx = [[0]*n for _ in range(n)]
    for (i, j) in arcs:
        Mx[i][j] = 1
    return Mx

def canon(arcs, n):
    """Canonical adjacency string over all relabelings. Returns a tuple key."""
    Mx = arcs_to_matrix(arcs, n)
    best = None
    for perm in permutations(range(n)):
        # relabel: new vertex k = perm[k]; build adjacency under this labeling
        flat = []
        for a in range(n):
            for b in range(n):
                flat.append(Mx[perm[a]][perm[b]])
        t = tuple(flat)
        if best is None or t < best:
            best = t
    return best

def is_tournament(arcs, n):
    """Check exactly one arc between each pair (no missing, no double)."""
    seen = {}
    for (i, j) in arcs:
        if i == j:
            return False
        key = frozenset((i, j))
        if key in seen:
            return False
        seen[key] = 1
    return len(seen) == n*(n-1)//2

def num_3cycles(arcs, n):
    A = arcs_to_matrix(arcs, n)
    c = 0
    for a, b, cc in combinations(range(n), 3):
        # count cyclic triangles among a,b,cc
        # a tournament triangle is cyclic iff not transitive
        edges = []
        for (x, y) in ((a, b), (b, cc), (cc, a),
                       (b, a), (cc, b), (a, cc)):
            pass
        # determine orientation
        def beats(x, y): return A[x][y] == 1
        # count out-degree within triple
        outs = {a: 0, b: 0, cc: 0}
        for (x, y) in combinations([a, b, cc], 2):
            if beats(x, y): outs[x] += 1
            else: outs[y] += 1
        # cyclic iff all out-degrees == 1
        if sorted(outs.values()) == [1, 1, 1]:
            c += 1
    return c

def is_transitive(arcs, n):
    return num_3cycles(arcs, n) == 0

# ----------------------------------------------------------------------
# Reference: full free set sizes (A000568) for n=3,4,5  -> 2,4,12
# We'll generate all iso classes for the free set by brute force for n<=5.
# ----------------------------------------------------------------------
def all_tournament_classes(n):
    """All iso classes of tournaments on n vertices (brute)."""
    pairs = list(combinations(range(n), 2))
    classes = {}
    for bits in product([0, 1], repeat=len(pairs)):
        arcs = []
        for (p, (i, j)) in zip(bits, pairs):
            if p == 0:
                arcs.append((i, j))
            else:
                arcs.append((j, i))
        key = canon(arcs, n)
        if key not in classes:
            classes[key] = num_3cycles(arcs, n)
    return classes  # key -> #3cycles

# ----------------------------------------------------------------------
# Helpers for continued fractions / Stern-Brocot
# ----------------------------------------------------------------------
def cont_frac(p, q):
    """Continued fraction expansion of p/q (p,q positive ints)."""
    a = []
    while q:
        a.append(p // q)
        p, q = q, p % q
    return a

def stern_brocot_path(p, q):
    """Path from root to p/q in Stern-Brocot tree as a string of 'L'/'R'.
    Uses continued fraction. p/q in lowest terms, positive rational."""
    gg = gcd(p, q); p //= gg; q //= gg
    cf = cont_frac(p, q)
    # standard: cf = [a0; a1, a2, ...]; path is R^{a0} L^{a1} R^{a2} ...
    # but last term reduced by 1 to terminate
    if len(cf) == 0:
        return ""
    cf = cf[:]
    cf[-1] -= 1
    path = []
    letter = 'R'
    for k, a in enumerate(cf):
        path.append(letter * a)
        letter = 'L' if letter == 'R' else 'R'
    return "".join(path)

# ======================================================================
# LRC INPUT GENERATION
# We need families of (S, tau) where S is a set of n distinct positive speeds,
# tau is the lonely-optimum time (from M), and the observer is "lonely enough".
# For exploration at small n we enumerate many primitive speed sets and use
# tau = argmax (the lonely time).  We also offer the "all crossing times" variant.
# ======================================================================
def speed_sets(n, maxspeed):
    """All n-subsets of {1..maxspeed} that are primitive (gcd 1) and distinct."""
    out = []
    for S in combinations(range(1, maxspeed+1), n):
        if gcd_list(S) == 1:
            out.append(S)
    return out

def gcd_list(xs):
    g0 = 0
    for x in xs:
        g0 = gcd(g0, x)
    return g0

# ======================================================================
# METHOD DEFINITIONS
# Each method: given (S, tau) produce arcs on n vertices (or fewer/more).
# Return None if the rule yields a tie / undefined -> skip that input.
# ======================================================================

# ---- METHOD 1: APPROXIMATION-QUALITY tournament ----------------------
# Vertices = runners. Arc i->j iff ||v_i tau|| < ||v_j tau||  (i is "lonelier",
#   closer to integer => better Diophantine approx of an integer by v_i tau).
# This is the OVERTAKING-style snapshot on the *unsigned* distance.
# Ties broken? If ||v_i||==||v_j|| skip (undefined). Expected: TOTAL ORDER => transitive.
def method1(S, tau):
    n = len(S)
    vals = [nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            if vals[i] == vals[j]:
                return None
            if vals[i] < vals[j]:
                arcs.append((i, j))
            else:
                arcs.append((j, i))
    return arcs

# ---- METHOD 2: SIGNED-SIDE + QUALITY (the "which side, then ratio") ---
# Arc i->j by a 2-key comparison that is NOT a global total order:
#   primary key: sign of signed_frac(v tau)  (which side of the integer:
#       ahead = positive, behind = negative). Group A = ahead, group B = behind.
#   Within a group: closer-to-integer beats farther (smaller ||.||).
#   ACROSS groups: the runner whose APPROX QUALITY (||v tau|| * v, the Diophantine
#       "badness" v*||v tau||) is SMALLER beats the other -- this cross term can
#       cycle because it mixes side with a different scalar.
# Concretely define a score key per runner: q_i = v_i * ||v_i tau||  (irrationality
#   measure proxy). Arc i->j iff:
#       same side: smaller ||.|| wins
#       diff side: smaller q wins
# Non-total because the comparison key changes by pairing.
def method2(S, tau):
    n = len(S)
    dist = [nrm(v*tau) for v in S]
    side = [1 if signed_frac(v*tau) > 0 else (-1 if signed_frac(v*tau) < 0 else 0) for v in S]
    q = [v * nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            if side[i] == side[j]:
                if dist[i] == dist[j]:
                    return None
                if dist[i] < dist[j]:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
            else:
                if q[i] == q[j]:
                    return None
                if q[i] < q[j]:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
    return arcs

# ---- METHOD 3: CONTINUED-FRACTION DEPTH of the speed RATIO -----------
# Vertices = runners. Arc i->j by comparing v_i/v_j to its "complexity".
# Rule: for ordered pair, look at ratio r = v_max/v_min for {i,j}. Compute the
#   continued fraction length L(i,j) = len(cont_frac(v_i, v_j)) (symmetric).
#   Direction: i->j iff (sum of CF partial quotients of v_i/v_j) has the SAME
#   parity as ... no. Make it depend on tau too:
#   Arc i->j iff the LAST convergent denominator q* of v_i/v_j that is <= 1/(2 tau-ish)
#   ... Simpler robust rule:
#   Define for the pair {i,j}: a = max(v_i,v_j), b = min. Let cf = cont_frac(a,b).
#   Let P = sum(cf) (sum of partial quotients), an integer >= 1.
#   Arc points from the SMALLER-INDEX-IN-SORTED-SPEED to larger if P is even,
#   reversed if P is odd. This is the "Stern-Brocot parity" orientation, tau-free.
#   (We also offer a tau-coupled variant 3b.)
def method3(S, tau):
    n = len(S)
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            a, b = max(S[i], S[j]), min(S[i], S[j])
            cf = cont_frac(a, b)
            P = sum(cf)
            # who is "bigger speed"?
            big = i if S[i] > S[j] else j
            sml = j if big == i else i
            if P % 2 == 0:
                arcs.append((sml, big))   # even: small speed beats big speed
            else:
                arcs.append((big, sml))
    return arcs

# ---- METHOD 3b: CF-depth coupled with tau (three-distance flavored) ----
# Arc i->j iff  floor(1/||v_i tau||) > floor(1/||v_j tau||)  with tie -> CF parity
#   of v_i/v_j. floor(1/||.||) is the "Ostrowski depth" = how many integer
#   multiples fit before the gap; mixing this discrete depth with the CF parity
#   tiebreak can cycle.
def method3b(S, tau):
    n = len(S)
    depth = []
    for v in S:
        d = nrm(v*tau)
        if d == 0:
            depth.append(None)  # parked runner: infinite depth
        else:
            depth.append(int(F(1,1)//d))  # floor(1/d)
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            di, dj = depth[i], depth[j]
            if di is None and dj is None:
                return None
            if di is None:  # i parked = infinite depth = i wins
                arcs.append((i, j)); continue
            if dj is None:
                arcs.append((j, i)); continue
            if di > dj:
                arcs.append((i, j))
            elif dj > di:
                arcs.append((j, i))
            else:
                # tiebreak by CF parity of speeds
                a, b = max(S[i], S[j]), min(S[i], S[j])
                P = sum(cont_frac(a, b))
                big = i if S[i] > S[j] else j
                sml = j if big == i else i
                if P % 2 == 0:
                    arcs.append((sml, big))
                else:
                    arcs.append((big, sml))
    return arcs

# ---- METHOD 4: STERN-BROCOT ANCESTRY tournament ----------------------
# Vertices = runners. Place each "phase point" p_i = frac(v_i tau) in (0,1) as a
#   rational. Arc i->j iff p_i is an ANCESTOR of p_j in the Stern-Brocot tree,
#   OR (if neither is ancestor of the other) by which subtree is "more left".
#   Ancestry is a PARTIAL order (forest), so to make a tournament we need a rule
#   for incomparable pairs. Rule:
#     - if SB-path(p_i) is a prefix of SB-path(p_j): i->j (ancestor beats descendant)
#     - elif SB-path(p_j) is a prefix of SB-path(p_i): j->i
#     - else (diverge): at the first differing letter, the one taking 'L' beats
#       the one taking 'R' (left beats right). This is a total order on leaves
#       (the in-order = numeric order) so divergent pairs are ordered by VALUE,
#       but ancestor/descendant pairs are ordered by DEPTH -> the mix can cycle.
def method4(S, tau):
    n = len(S)
    paths = []
    for v in S:
        p = nrm_frac01(v*tau)
        if p == 0:
            return None   # parked at integer; SB undefined for 0
        paths.append(stern_brocot_path(p.numerator, p.denominator))
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = paths[i], paths[j]
            if pi == pj:
                return None
            if pj.startswith(pi):     # i ancestor of j
                arcs.append((i, j))
            elif pi.startswith(pj):   # j ancestor of i
                arcs.append((j, i))
            else:
                # first divergence
                k = 0
                while pi[k] == pj[k]:
                    k += 1
                # left beats right
                if pi[k] == 'L':
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
    return arcs

def nrm_frac01(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r  # frac in [0,1)

# ---- METHOD 5: BINDING / CROSSING-PARITY tournament ------------------
# This is the one most aligned with THM-524 (binding pairs). For each pair {i,j},
#   consider the crossing time of the pair: the rule looks at HOW the pair binds
#   at tau. Arc i->j iff at the optimum tau the runner i is "tighter against the
#   wall" measured by the crossing structure:
#   Define for pair {i,j}: the nearest binding crossing value among
#     k/(v_i+v_j) and k/(v_j-v_i). Let c = the crossing time closest to tau.
#   Arc i->j iff v_i*tau is on the FAR side (frac > 1/2) of its nearest half-integer
#     while v_j is on the near side -- i.e., compare parity floor(2 v_i tau) vs
#     floor(2 v_j tau). This counts how many half-integers each runner has passed.
#   Arc i->j iff floor(2 v_i tau) ... no, that's monotone in v. Use PARITY:
#     b_i = floor(2 v_i tau) mod 2  (which "lane", even/odd half-integer band).
#     If b_i != b_j: the one with b=0 (in an even band, closer side) beats b=1.
#     If b_i == b_j: tiebreak by larger v beats smaller? -> would be transitive.
#       Instead tiebreak by the CROSSING denominator: smaller (v_i+v_j vs ...)...
#   To get genuine cycles we use a 3-way-rock-paper rule:
#     c_i = (v_i * tau) mod 1 mapped to one of 3 thirds: which third of [0,1).
#     Arc i->j iff (third_i - third_j) mod 3 == 1  (cyclic rock-paper-scissors),
#       tie (same third) -> closer-to-integer wins.
#   Three thirds + cyclic rule is a classic cycle generator -> definitely H>1.
def method5(S, tau):
    n = len(S)
    third = []
    for v in S:
        f = nrm_frac01(v*tau)  # [0,1)
        if f < F(1,3): third.append(0)
        elif f < F(2,3): third.append(1)
        else: third.append(2)
    dist = [nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            ti, tj = third[i], third[j]
            if ti == tj:
                if dist[i] == dist[j]:
                    return None
                if dist[i] < dist[j]:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
            else:
                # cyclic: 0 beats 2, 2 beats 1, 1 beats 0  (i.e. (a-b) mod 3 ==2 -> a beats b)
                if (ti - tj) % 3 == 2:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
    return arcs

# ---- METHOD 6: THREE-DISTANCE (Steinhaus) GAP tournament -------------
# Vertices = runners. At the lonely tau, the points {frac(v_i tau)} partition the
#   circle. By the three-distance theorem the GAPS take at most 3 lengths. Orient
#   by gap structure: sort points around circle; each runner has a LEFT gap and a
#   RIGHT gap (to neighbors). Arc i->j iff runner i's right-gap is SHORTER than
#   runner j's right-gap (more crowded on its right), tie -> by point value.
#   But "right gap" mixes circular neighbor relations -> can be non-total.
#   Actually we define: g_i = length of the arc from point_i to the NEXT point
#   clockwise. Arc i->j iff g_i < g_j (i more crowded), tie-> smaller ||v tau||.
#   Because g_i depends on circular neighbors not on i alone, but the comparison
#   g_i<g_j is still a total order on the multiset... so likely transitive UNLESS
#   ties chain. We add a twist: i->j iff g_i < g_j AND point_i < point_j is NOT
#   required -> just g comparison (total) => transitive. Mark and test.
def method6(S, tau):
    n = len(S)
    pts = [nrm_frac01(v*tau) for v in S]
    # if any coincide, skip
    if len(set(pts)) != n:
        return None
    order = sorted(range(n), key=lambda k: pts[k])
    # clockwise gap from pts[k] to next
    gap = {}
    for idx in range(n):
        cur = order[idx]
        nxt = order[(idx+1) % n]
        if idx < n-1:
            gap[cur] = pts[nxt] - pts[cur]
        else:
            gap[cur] = (1 - pts[cur]) + pts[nxt]
    dist = [nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            if gap[i] == gap[j]:
                if dist[i] == dist[j]:
                    return None
                if dist[i] < dist[j]:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
            else:
                if gap[i] < gap[j]:
                    arcs.append((i, j))
                else:
                    arcs.append((j, i))
    return arcs

# ======================================================================
# DRIVER
# ======================================================================
METHODS = {
    "M1_approx_quality": method1,
    "M2_side_then_quality": method2,
    "M3_CF_parity_speed": method3,
    "M3b_Ostrowski_depth": method3b,
    "M4_stern_brocot": method4,
    "M5_thirds_rockpaper": method5,
    "M6_three_distance_gap": method6,
}

def realized_classes(method, n, maxspeed, tau_mode="argmax"):
    """
    Run `method` over all primitive n-speed sets in 1..maxspeed.
    tau_mode='argmax': use the single lonely-optimum tau (M(S)[1]).
    tau_mode='allcand': use EVERY candidate crossing time (richer family).
    Returns: dict class_key -> #3cycles realized, plus counts and a few stats.
    """
    realized = {}            # key -> #3cycles
    n_inputs = 0
    n_nontrivial_tournaments = 0
    n_skipped = 0
    sample_nonpaths = []
    for S in speed_sets(n, maxspeed):
        if tau_mode == "argmax":
            _, tau = M(S)
            taus = [tau]
        else:
            taus = sorted(cand(S))
        for tau in taus:
            n_inputs += 1
            arcs = method(S, tau)
            if arcs is None:
                n_skipped += 1
                continue
            if not is_tournament(arcs, n):
                n_skipped += 1
                continue
            c3 = num_3cycles(arcs, n)
            if c3 > 0:
                n_nontrivial_tournaments += 1
            key = canon(arcs, n)
            if key not in realized:
                realized[key] = c3
                if c3 > 0 and len(sample_nonpaths) < 3:
                    sample_nonpaths.append((S, tau, c3))
    return realized, n_inputs, n_nontrivial_tournaments, n_skipped, sample_nonpaths

def main():
    print("="*70)
    print("DIOPHANTINE TOURNAMENT MAPS for LRC(14)")
    print("="*70)

    # Free sets
    free = {}
    for n in (3, 4, 5):
        free[n] = all_tournament_classes(n)
        print(f"FREE set n={n}: {len(free[n])} iso classes "
              f"(A000568 expects {[None,1,1,2,4,12][n]}); "
              f"3cycle-counts present: {sorted(set(free[n].values()))}")
    print()

    # quick non-triviality gate at n=3, both tau modes, modest maxspeed
    print("-"*70)
    print("NON-TRIVIALITY GATE (n=3): does the map ever make a 3-cycle (H>1)?")
    print("-"*70)
    gate = {}
    for name, fn in METHODS.items():
        rc, ni, nt, sk, samp = realized_classes(fn, 3, 14, tau_mode="allcand")
        has_cycle = any(c > 0 for c in rc.values())
        gate[name] = has_cycle
        print(f"{name:26s}: realized {len(rc)}/2 classes, "
              f"nontrivial-tournaments={nt}, ever-3cycle={has_cycle}")
    print()

    # Full realized-class census for nontrivial maps at n=3,4,5
    print("="*70)
    print("REALIZED-CLASS CENSUS (nontrivial maps only)")
    print("  search family: all primitive n-speed sets in 1..MAXSPEED, "
          "BOTH tau modes")
    print("="*70)
    MAXSPEED = {3: 16, 4: 14, 5: 12}   # keep n=5 tractable
    for name, fn in METHODS.items():
        if not gate[name]:
            print(f"\n### {name}: TRANSITIVE-ONLY at n=3 (dead). Skipping census.")
            continue
        print(f"\n### {name}")
        for n in (3, 4, 5):
            ms = MAXSPEED[n]
            # argmax mode
            rc_a, ni_a, nt_a, sk_a, samp_a = realized_classes(fn, n, ms, "argmax")
            # allcand mode
            rc_c, ni_c, nt_c, sk_c, samp_c = realized_classes(fn, n, ms, "allcand")
            combined = dict(rc_a); combined.update(rc_c)
            nfree = len(free[n])
            forbidden = [k for k in free[n] if k not in combined]
            print(f"  n={n} (maxspeed {ms}): "
                  f"realized {len(combined)}/{nfree} free classes  "
                  f"[argmax {len(rc_a)}, allcand {len(rc_c)}]")
            print(f"        nontrivial tournaments: argmax {nt_a}, allcand {nt_c}; "
                  f"inputs(allcand) {ni_c}, skipped {sk_c}")
            if forbidden:
                fb_c3 = sorted(free[n][k] for k in forbidden)
                print(f"        FORBIDDEN classes: {len(forbidden)} "
                      f"(their 3cycle-counts: {fb_c3})")
            else:
                print(f"        FORBIDDEN classes: none (realizes everything)")
            if samp_c:
                print(f"        sample nontrivial (S,tau,#3cyc): {samp_c}")

if __name__ == "__main__":
    main()
