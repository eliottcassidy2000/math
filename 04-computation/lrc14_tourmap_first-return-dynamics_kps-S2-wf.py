# lrc14_tourmap_first-return-dynamics_kps-S2-wf.py
# kind-pasteur-S2-wf : Creative tournament-generation under theme "first-return dynamics".
#
# Goal: map convoluted aspects of LRC to tournament structures, test NON-TRIVIALITY
# (does the map ever produce a non-transitive tournament, H>1?), and for the non-trivial
# ones, enumerate which iso classes are realized as (S, tau) range over LRC-constrained
# inputs -- compared to the full "free" set (A000568). A class never realized = FORBIDDEN.
#
# EXACT rationals everywhere (fractions.Fraction). No floats for decisions.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

# ----------------------------------------------------------------------------
# LRC primitives (validated tool, copied verbatim)
# ----------------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def frac(x):
    # fractional part in [0,1)
    r = x - int(x); return r + 1 if r < 0 else r

def g(S, t):
    return min(nrm(v * t) for v in S)

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

def Mgap(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

# ----------------------------------------------------------------------------
# Tournament canonicalization tools
# ----------------------------------------------------------------------------
# A tournament on k labeled vertices is a function arc(i,j) in {True,False}:
# arc(i,j)=True means i->j. We require it total + antisymmetric.
# Represent as frozenset of directed edges (i,j) with i->j.

def is_tournament(n, edges):
    for i in range(n):
        for j in range(i+1, n):
            if ((i,j) in edges) == ((j,i) in edges):
                return False
    return True

def num_3cycles(n, edges):
    # count directed 3-cycles
    cnt = 0
    for a,b,c in combinations(range(n), 3):
        # check the 6 orderings collapse: a 3-cycle is a->b->c->a or a->c->b->a
        def arc(x,y): return (x,y) in edges
        if arc(a,b) and arc(b,c) and arc(c,a): cnt += 1
        if arc(a,c) and arc(c,b) and arc(b,a): cnt += 1
    return cnt

def ham_paths(n, edges):
    # count Hamiltonian paths (Redei). H is always odd.
    arc = lambda x,y: (x,y) in edges
    cnt = 0
    for perm in permutations(range(n)):
        ok = all(arc(perm[k], perm[k+1]) for k in range(n-1))
        if ok: cnt += 1
    return cnt

def score_seq(n, edges):
    sc = [0]*n
    for (i,j) in edges:
        sc[i] += 1
    return tuple(sorted(sc))

def canon(n, edges):
    # brute canonical form over all n! relabelings: pick lexicographically minimal
    # adjacency tuple of out-edges.
    best = None
    verts = list(range(n))
    for perm in permutations(verts):
        # perm[k] = new label of old vertex k ; build relabeled edge set
        # we want a canonical string: for new vertices 0..n-1, the upper triangle orientation
        inv = [0]*n
        for k,p in enumerate(perm): inv[p] = k
        # new edge (a,b) exists iff old edge (inv? ) ... use: old i->j becomes perm[i]->perm[j]
        mat = [[0]*n for _ in range(n)]
        for (i,j) in edges:
            mat[perm[i]][perm[j]] = 1
        key = tuple(mat[a][b] for a in range(n) for b in range(n) if a<b)
        if best is None or key < best:
            best = key
    return best

# Precompute the full free set of iso classes for small n by enumerating all tournaments.
def all_tournament_classes(n):
    pairs = list(combinations(range(n), 2))
    classes = {}
    for bits in range(2**len(pairs)):
        edges = set()
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: edges.add((i,j))
            else: edges.add((j,i))
        c = canon(n, edges)
        if c not in classes:
            classes[c] = (num_3cycles(n,edges), ham_paths(n,edges), score_seq(n,edges))
    return classes

FREE = {}
for n in (3,4,5):
    FREE[n] = all_tournament_classes(n)
A000568 = {3:2, 4:4, 5:12}
for n in (3,4,5):
    assert len(FREE[n])==A000568[n], (n, len(FREE[n]), A000568[n])

# Helper: from an arc oracle (function arc(i,j)->bool, total/antisym on subset verts)
# build edge set and canonicalize. Returns (canon, n3, H, score) or None if not a tournament.
def classify(verts, arcfn):
    n = len(verts)
    idx = {v:k for k,v in enumerate(verts)}
    edges = set()
    for a in range(n):
        for b in range(a+1, n):
            va, vb = verts[a], verts[b]
            r = arcfn(va, vb)  # True => va->vb, False => vb->va, None => undefined
            if r is None: return None
            if r: edges.add((a,b))
            else: edges.add((b,a))
    if not is_tournament(n, edges): return None
    return canon(n, edges), num_3cycles(n,edges), ham_paths(n,edges), score_seq(n,edges)

# ----------------------------------------------------------------------------
# LRC input generators (constrained families)
# ----------------------------------------------------------------------------
def primitive(S):
    from functools import reduce
    return reduce(gcd, S) == 1

def gen_speed_sets(n, maxspeed):
    # all primitive n-subsets of {1..maxspeed}
    out = []
    for S in combinations(range(1, maxspeed+1), n):
        if primitive(S):
            out.append(S)
    return out

# ============================================================================
# METHOD 1: FIRST-HITTING-OF-SECTION-0 ORDER (continuous, off-grid).
#   Vertices = runners. As tau increases from 0+, runner v first hits "near 0"
#   (||v tau||=0, i.e. v tau in Z) at tau = 1/v. So first-hit time of integer is 1/v:
#   that's a TOTAL ORDER by speed -> transitive. DEAD (like overtaking).
#   BUT: first-hit of a HALF-INTEGER wall (||v tau||=1/2) is at tau=1/(2v): also total order.
#   So "first hitting of a fixed point" is transitive. We instead use FIRST-RETURN to a
#   MOVING reference: pairwise first MEETING time. Arc i->j iff i laps/meets j-relative...
#   Method 1 proper: PAIRWISE FIRST-MEETING PARITY (see below). The naive single-point
#   first-hit is recorded as M1a (transitive, dead).
# ============================================================================

def M1a_arc(vi, vj):
    # first time ||v tau||=0 is tau=1/v; smaller speed hits later. i->j iff i hits 0 FIRST
    # i.e. 1/vi < 1/vj  <=> vi>vj. Total order. (Transitive.)
    return vi > vj

# ============================================================================
# METHOD 2: PAIRWISE FIRST-MEETING AT A FIXED LONELY TIME tau (hitting-time dice style).
#   At the lonely time tau*, each runner v has a phase phi_v = frac(v*tau*).
#   Define a relative angular "first to return to phase 0 going forward" between pairs.
#   Runner v returns to integer position at multiples of 1/v. The RELATIVE first-meeting
#   of i and j (when frac(vi t)=frac(vj t) for t>tau*) is governed by |vi-vj|. This is
#   symmetric in i,j -> need orientation. We orient by WHICH runner is ahead in phase at
#   the meeting. Concretely: arc i->j iff at tau*, runner i is "chasing" j and catches up,
#   determined by sign of (vi-vj) combined with phase order. Tested below.
# ============================================================================

def M2_arc_factory(S, tau):
    phi = {v: frac(v*tau) for v in S}
    def arc(vi, vj):
        # relative speed and relative phase; "first meeting going forward in time"
        # relative position r(t) = frac((vi-vj) t + (phi_i - phi_j))... orient by who closes the gap.
        d = vi - vj  # >0 means i faster
        gap = frac(phi[vi] - phi[vj])
        # time for i to overtake j (close gap to 0 from above) if d>0: t = gap/d ; if d<0: t=(1-gap)/(-d)...
        # define "first to meet" winner = the one that reaches the meeting as the leader.
        # Simplest non-snapshot orientation: i->j iff (vi - vj) and (phi_i - phi_j over 1/2) agree.
        if d > 0:
            tmeet = gap / d if gap != 0 else F(1, abs(d))  # gap=0 means already meeting; use full lap
        else:
            tmeet = (1 - gap) / (-d) if gap != 0 else F(1, abs(d))
        # winner: the faster runner arrives "as leader" -> i->j iff vi>vj XOR (phase says j leads)
        # Use tmeet parity-free orientation: i->j iff vi>vj. (still need a real twist) -- see M3.
        return vi > vj
    return arc

# ============================================================================
# METHOD 3: LAP-COUNT PARITY TOURNAMENT (overtaking COUNTS within a period).
#   Over one full period of the joint motion (period = 1 in tau since speeds integer,
#   everything periodic with period 1), runner i laps runner j a net |vi-vj| times.
#   That's symmetric-ish. Instead count, on the interval (0, tau*], how many times
#   runner i crosses runner j (i.e. frac(vi t)=frac(vj t)): that's the number of
#   integer multiples of 1/|vi-vj| in (0,tau*], = floor(|vi-vj| * tau*).
#   Arc i->j iff this crossing-count parity is ODD (i is "currently ahead" in the
#   meeting cycle), with a tie-break. THIS is a genuine hitting/parity structure.
#   Orientation: i->j iff ( (vi>vj) XOR (floor(|vi-vj|*tau*) is odd) ).
# ============================================================================

def M3_arc_factory(S, tau):
    def arc(vi, vj):
        d = abs(vi - vj)
        crossings = (d * tau)  # Fraction
        k = int(crossings)     # floor since crossings>=0
        odd = (k % 2 == 1)
        base = (vi > vj)
        return base ^ odd
    return arc

# ============================================================================
# METHOD 4: SECTION-0 FIRST-RETURN ORDER ON THE GRID (dihedral clock / wall switches).
#   Use the on-grid lens. For a grid time tau=a/14 (gcd(a,14)=1) runner v sits in
#   section r_v = v*a mod 14. Define a DYNAMIC: as the grid index a runs through the
#   units (Z/14)* in their natural cyclic/clock order 1,3,5,9,11,13, track for each
#   pair (i,j) WHO ENTERS A FIXED TARGET SECTION (say section 1, the "wall") FIRST.
#   For runner v, the set of grid times a with v*a ≡ 1 (mod 14) is a single residue
#   class (when gcd(v,14)=1): a ≡ v^{-1} mod 14. So "first a (in clock order) at which
#   v reaches the wall" = position of v^{-1} in the clock order. Arc i->j iff i reaches
#   the wall before j. Total order by clock-position of inverse -> likely transitive.
#   We TEST. We also test the variant: first a at which v reaches ANY of two walls
#   {1,13} -> min over two residue classes (can break transitivity).
# ============================================================================

CLOCK = [1,3,5,9,11,13]  # natural listing of (Z/14)* ; we test several orders
def inv14(v):
    v %= 14
    for x in range(14):
        if (v*x)%14==1: return x
    return None

def M4_arc_factory(S, wall_set, order):
    pos = {a:i for i,a in enumerate(order)}
    def firsttime(v):
        # first a in `order` with v*a ≡ w (mod14) for some w in wall_set
        best = None
        for a in order:
            if (v*a)%14 in wall_set:
                best = pos[a]; break
        return best
    def arc(vi, vj):
        fi, fj = firsttime(vi), firsttime(vj)
        if fi is None or fj is None: return None
        if fi == fj: return None  # tie -> undefined (skip this set)
        return fi < fj
    return arc

# ============================================================================
# METHOD 5: PAIRWISE BINDING-CROSSING ORDER (THM-524 dynamics).
#   The gap is attained at a binding-pair crossing tau=k/(vi±vj). For each PAIR of
#   runners, the *earliest* crossing wall they share in (0,1/2] is tau_pair = 1/(vi+vj)
#   (first sum-crossing) or 1/(vj-vi). Build tournament on RUNNERS: i->j iff the first
#   "sum wall" 1/(vi+vj)... no, that's pairwise (one value per pair, not an orientation).
#   Real construction: vertex set = runners; orient i->j by comparing each runner's
#   OWN first half-integer wall time 1/(2v) against the SHARED sum-wall -- i.e. who is
#   closer to its private wall when the pair-wall fires. This is the FIRST-RETURN race
#   between a runner's solo clock (period 1/(2v)) and the pair clock. Tested.
#   arc i->j iff frac( (vi+vj) * (1/(2 vi)) ) < frac( (vi+vj)*(1/(2 vj)) )  -- a genuine
#   cross ratio. Equivalent to comparing (vi+vj)/(2vi) vs (vi+vj)/(2vj) fractional parts.
# ============================================================================

def M5_arc(vi, vj):
    a = frac(F(vi+vj, 2*vi))
    b = frac(F(vi+vj, 2*vj))
    if a == b: return None
    return a < b

# ============================================================================
# METHOD 6: RETURN-TIME-TO-NEIGHBORHOOD-OF-LONELY-PHASE (hitting-time dice, the real one).
#   At lonely time tau*, observer phase is 0 and each runner v has phase phi_v=frac(v tau*).
#   Consider the FORWARD return dynamics: increase t past tau*. Runner v next reaches
#   the SAME phase value phi_v again only after period 1/v... but we instead ask: starting
#   from tau*, who FIRST RE-ENTERS section 0 (||v t||<1/14, i.e. comes within 1/14 of an
#   integer). Runner v is within 1/14 of integer when frac(v t) in [0,1/14] U [13/14,1).
#   First t>tau* with v entering that band: determined by phi_v and v. The runner whose
#   band-entry time is earliest "returns first". Orient i->j iff i returns to section 0
#   before j (off-grid, continuous, genuine hitting time). This is the classic
#   non-transitive hitting-time tournament. Tested with band radius 1/14.
# ============================================================================

def M6_arc_factory(S, tau, radius=F(1,14)):
    # next time runner v enters band |frac(v t) - 0| <= radius (mod 1) for t>tau.
    # frac(v t) decreasing? It's increasing in t at rate v. band is [0,radius] U [1-radius,1).
    # entry into [0,radius] happens when v t crosses an integer from below... we want the
    # next t>tau where frac(v t) JUST enters [1-radius,1)U[0,radius]; entry point of the
    # union (going upward in frac) is at frac=1-radius. So we want next t with frac(v t)=1-radius.
    def nexttime(v):
        # solve frac(v t)=1-radius for smallest t>tau: v t = m + (1-radius) for integer m
        target = 1 - radius
        # t = (m + target)/v ; smallest m with t>tau => m > v*tau - target
        m = (v*tau - target)
        mm = int(m)
        while F(mm) <= m: mm += 1  # smallest integer strictly greater
        # ensure t>tau strictly
        t = F(mm) + target
        t = t / v
        while t <= tau:
            mm += 1; t = (F(mm)+target)/v
        return t
    def arc(vi, vj):
        ti, tj = nexttime(vi), nexttime(vj)
        if ti == tj: return None
        return ti < tj  # i returns first => i->j
    return arc

# ----------------------------------------------------------------------------
# DRIVER: for each method, run non-triviality gate + class census on small n.
# ----------------------------------------------------------------------------
def census_method_on_speedsets(name, arc_builder, n, speed_sets, use_tau=False,
                                tau_mode='M'):
    """arc_builder: if use_tau False -> function arc(vi,vj). If use_tau True ->
    factory(S,tau)->arc. tau_mode 'M' uses the gap-optimal tau; 'all' uses all candidate taus."""
    realized = {}  # canon -> (n3, H, score, example)
    nontrivial = False
    examples_3cyc = 0
    total_evaluated = 0
    for S in speed_sets:
        Ssort = sorted(S)
        if use_tau:
            if tau_mode == 'M':
                _, tau = Mgap(Ssort)
                taus = [tau] if tau is not None else []
            else:
                taus = [t for t in cand(Ssort) if t>0]
            for tau in taus:
                arc = arc_builder(Ssort, tau)
                res = classify(Ssort, arc)
                total_evaluated += 1
                if res is None: continue
                c,n3,H,score = res
                if H>1: nontrivial = True; examples_3cyc += 1
                if c not in realized:
                    realized[c] = (n3,H,score,(tuple(S),tau))
        else:
            arc = arc_builder
            res = classify(Ssort, arc)
            total_evaluated += 1
            if res is None: continue
            c,n3,H,score = res
            if H>1: nontrivial = True; examples_3cyc += 1
            if c not in realized:
                realized[c] = (n3,H,score,(tuple(S),None))
    free = FREE[n]
    forbidden = [c for c in free if c not in realized]
    return {
        'name': name, 'n': n, 'realized_count': len(realized),
        'free_count': len(free), 'forbidden': forbidden, 'nontrivial': nontrivial,
        'examples_3cyc': examples_3cyc, 'total_evaluated': total_evaluated,
        'realized': realized,
    }

def summarize(res):
    print(f"\n=== {res['name']}  (n={res['n']}) ===")
    print(f"  evaluated tournaments: {res['total_evaluated']}")
    print(f"  realized classes: {res['realized_count']} / free {res['free_count']}")
    print(f"  NON-TRIVIAL (some H>1): {res['nontrivial']}  (#3-cycle tournaments seen: {res['examples_3cyc']})")
    # H distribution of realized classes
    Hs = sorted(set(v[1] for v in res['realized'].values()))
    print(f"  realized H values: {Hs}")
    if res['forbidden']:
        # describe forbidden by (n3,H,score)
        fdesc = []
        for c in res['forbidden']:
            n3,H,score = FREE[res['n']][c]
            fdesc.append((H,n3,score))
        fdesc.sort()
        print(f"  FORBIDDEN classes ({len(res['forbidden'])}): " +
              "; ".join(f"H={H},c3={n3},score={score}" for H,n3,score in fdesc))
    else:
        print("  FORBIDDEN: none (realizes the full free set)")

if __name__ == '__main__':
    print("LRC-14 first-return-dynamics tournament maps  (kind-pasteur-S2-wf)")
    print("="*70)
    print("Free set sizes (A000568):", {n:len(FREE[n]) for n in (3,4,5)})

    RESULTS = []

    for n in (3,4,5):
        maxspeed = {3:14, 4:14, 5:12}[n]   # keep enumeration tractable but rich
        SS = gen_speed_sets(n, maxspeed)
        print(f"\n--- n={n}: {len(SS)} primitive speed-sets from 1..{maxspeed} ---")

        # METHOD 1a (control, expected transitive/dead)
        r = census_method_on_speedsets("M1a first-hit-section0 (speed order)",
                                        M1a_arc, n, SS, use_tau=False)
        summarize(r); RESULTS.append(r)

        # METHOD 3 lap-count parity (needs tau)
        r = census_method_on_speedsets("M3 lap-count-parity @gap-tau",
                                        M3_arc_factory, n, SS, use_tau=True, tau_mode='M')
        summarize(r); RESULTS.append(r)
        r = census_method_on_speedsets("M3 lap-count-parity @ALL-cand-tau",
                                        M3_arc_factory, n, SS, use_tau=True, tau_mode='all')
        summarize(r); RESULTS.append(r)

        # METHOD 5 cross-ratio first-return (no tau)
        r = census_method_on_speedsets("M5 cross-ratio solo-vs-pair wall",
                                        M5_arc, n, SS, use_tau=False)
        summarize(r); RESULTS.append(r)

        # METHOD 6 return-to-section0 hitting-time dice (needs tau)
        r = census_method_on_speedsets("M6 return-to-section0 hitting @gap-tau",
                                        M6_arc_factory, n, SS, use_tau=True, tau_mode='M')
        summarize(r); RESULTS.append(r)
        r = census_method_on_speedsets("M6 return-to-section0 hitting @ALL-cand-tau",
                                        M6_arc_factory, n, SS, use_tau=True, tau_mode='all')
        summarize(r); RESULTS.append(r)

    # METHOD 4 grid dihedral-clock first-to-wall (n up to 5, speeds coprime to 14 for inverses)
    print("\n" + "="*70)
    print("METHOD 4 (grid dihedral-clock, units mod 14)")
    for n in (3,4,5):
        # speeds must be units mod 14 for inverse to exist; use {1,3,5,9,11,13} and beyond coprime to14
        coprime14 = [v for v in range(1,28) if gcd(v,14)==1]
        SS = [S for S in combinations(coprime14, n) if primitive(S)]
        for wall_set in ({1}, {1,13}, {1,3}):
            for ord_name, order in (("nat", CLOCK), ("sorted", sorted(CLOCK))):
                realized={}; nontrivial=False; ex=0; tot=0
                for S in SS:
                    Ssort=sorted(S)
                    arc=M4_arc_factory(Ssort, wall_set, order)
                    res=classify(Ssort, arc)
                    tot+=1
                    if res is None: continue
                    c,n3,H,score=res
                    if H>1: nontrivial=True; ex+=1
                    if c not in realized: realized[c]=(n3,H,score)
                free=FREE[n]; forb=[c for c in free if c not in realized]
                print(f"  n={n} wall={sorted(wall_set)} order={ord_name}: "
                      f"realized {len(realized)}/{len(free)} nontrivial={nontrivial} "
                      f"(#3cyc={ex}, evaluated={tot}) forbidden={len(forb)}")

    # ---- Focused n=4,5 deep dive on the most promising non-trivial method ----
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print(f"{'method':52s} {'n':>2s} {'real':>5s} {'free':>5s} {'nontriv':>8s} {'forb':>5s}")
    for r in RESULTS:
        print(f"{r['name']:52s} {r['n']:>2d} {r['realized_count']:>5d} "
              f"{r['free_count']:>5d} {str(r['nontrivial']):>8s} {len(r['forbidden']):>5d}")
