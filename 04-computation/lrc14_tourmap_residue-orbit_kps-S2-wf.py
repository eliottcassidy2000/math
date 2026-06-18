"""
lrc14_tourmap_residue-orbit_kps-S2-wf.py

THEME: residue-orbit tournament maps for the Lonely Runner Conjecture (LRC, N=14).

We map LRC data (speed sets S, lonely time tau, sections mod 14) to a TOURNAMENT
(complete directed graph on the runners / pairs / sections), then ask:
  (1) NON-TRIVIALITY: does the map ever produce a non-transitive tournament (H>1)?
  (2) FORBIDDEN CLASSES: which iso classes are realized as (S,tau) range over
      LRC-constrained inputs, vs the full "free" set of all A000568 iso classes?

All arithmetic exact (Fraction). No floats in decisions.

A000568 (number of tournament iso classes by n): n=1..7 -> 1,1,2,4,12,56,456
"""

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

# ---------------- EXACT M TOOL (verbatim from validated prompt) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def gmin(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

def all_opt_taus(S):
    """All tau in candidate set achieving the max gap (the lonely optima)."""
    b, _ = M(S)
    return sorted([t for t in cand(S) if gmin(S, t) == b]), b

# ---------------- TOURNAMENT CANONICALIZATION ----------------
# A tournament on m vertices is a set of arcs; we store as adjacency matrix A[i][j]=1 if i->j.

def canon_key(adj, m):
    """Canonical form: min over all relabelings of the flattened upper structure.
    adj is dict/2D: adj[i][j] in {0,1}, with adj[i][j]+adj[j][i]=1 for i!=j."""
    best = None
    verts = list(range(m))
    for perm in permutations(verts):
        # build tuple of arcs under relabeling: new vertex order = perm
        bits = []
        for i in range(m):
            for j in range(m):
                if i != j:
                    bits.append(adj[perm[i]][perm[j]])
        t = tuple(bits)
        if best is None or t < best:
            best = t
    return best

def h_count(adj, m):
    """Number of Hamiltonian paths (directed). Redei: always odd."""
    cnt = 0
    for perm in permutations(range(m)):
        ok = True
        for k in range(m - 1):
            if adj[perm[k]][perm[k + 1]] != 1:
                ok = False; break
        if ok: cnt += 1
    return cnt

def score_seq(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))

def is_valid_tournament(adj, m):
    for i in range(m):
        for j in range(i + 1, m):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True

def num_3cycles(adj, m):
    c = 0
    for a, b, cc in combinations(range(m), 3):
        # a->b->c->a or a->c->b->a
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c

# Reference: enumerate ALL iso classes on m vertices (the free set).
def all_iso_classes(m):
    seen = {}
    # iterate over all 2^C(m,2) tournaments
    pairs = list(combinations(range(m), 2))
    for bits in product([0, 1], repeat=len(pairs)):
        adj = [[0] * m for _ in range(m)]
        for (i, j), b in zip(pairs, bits):
            if b: adj[i][j] = 1
            else: adj[j][i] = 1
        k = canon_key(adj, m)
        if k not in seen:
            seen[k] = (h_count(adj, m), num_3cycles(adj, m), score_seq(adj, m))
    return seen

# ---------------- LRC INPUT FAMILIES (small) ----------------
# We generate small speed sets S of size m (m = number of runners = vertices for
# runner-vertex maps). LRC uses positive integer speeds, primitive (gcd=1).
# For tractable iso-class enumeration we keep m in {3,4,5}.

def primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

def gen_speed_sets(m, maxspeed):
    """All m-subsets of {1..maxspeed}, primitive, distinct positive speeds."""
    out = []
    for S in combinations(range(1, maxspeed + 1), m):
        if primitive(S):
            out.append(S)
    return out

# ============================================================
# METHOD DEFINITIONS
# Each method: given (S, modulus) produces a tournament adjacency on some vertex set.
# We aggregate over the unit-orbit (Z/mod)* — NOT a single time.
# ============================================================

def units_mod(mod):
    return [a for a in range(1, mod) if gcd(a, mod) == 1]

# ---- METHOD 1: Orbit-majority section comparison (runners as vertices) ----
# Arc i->j iff  #{a in (Z/mod)* : (v_i*a mod mod) < (v_j*a mod mod)}  >
#               #{a in (Z/mod)* : (v_j*a mod mod) < (v_i*a mod mod)}.
# Ties broken (rare for distinct residues) by raw speed.
def method1_orbit_majority(S, mod):
    m = len(S)
    U = units_mod(mod)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            wi = wj = 0
            for a in U:
                ri = (S[i] * a) % mod
                rj = (S[j] * a) % mod
                if ri < rj: wi += 1
                elif rj < ri: wj += 1
            if wi > wj:
                adj[i][j] = 1
            elif wj > wi:
                adj[j][i] = 1
            else:
                # tie: break by speed (deterministic, total order -> no help to nontriv)
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

# ---- METHOD 2: Quadratic-residue character of speed difference ----
# Requires mod = odd prime p. chi(x) = Legendre symbol. Arc i->j iff chi(v_i - v_j)=+1.
# (chi is antisymmetric: chi(-x) = chi(-1)chi(x). For p=3 mod 4, chi(-1)=-1 so this is a
#  well-defined tournament: exactly one of chi(v_i-v_j), chi(v_j-v_i) is +1 when nonzero.)
def legendre(a, p):
    a %= p
    if a == 0: return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1

def method2_qr_char(S, p):
    # p must be prime, p ≡ 3 mod 4 for a clean tournament (chi(-1)=-1).
    m = len(S)
    adj = [[None] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            d = (S[i] - S[j]) % p
            if d == 0:
                return None  # collision mod p -> undefined; skip this input
            l = legendre(d, p)
            if l == 1:
                adj[i][j] = 1; adj[j][i] = 0
            else:
                adj[i][j] = 0; adj[j][i] = 1
    for i in range(m): adj[i][i] = 0
    return adj

# ---- METHOD 3: Lonely-time section-position aggregated over the optimal-tau set ----
# At a lonely optimum tau (gap=M), each runner i has a fractional position frac(v_i*tau)
# in [0,1). We DON'T use a single snapshot order (that's the dead transitive map).
# Instead: over ALL optimal taus T*, arc i->j iff
#   #{tau in T* : frac(v_i tau) and frac(v_j tau) on OPPOSITE sides of 1/2,
#                 with i in lower half} ... aggregated by majority of a SWITCH parity.
# Concretely: define for each optimal tau the "wrap count" w_i(tau)=floor(2*frac(v_i tau))
# (0 if in [0,1/2), 1 if in [1/2,1)). Arc i->j iff sum over T* of (w_i - w_j) ... majority.
# This is a parity/half-plane aggregation, NOT a total order snapshot.
def method3_halfplane_majority(S, mod=None):
    taus, gap = all_opt_taus(S)
    m = len(S)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            score = 0
            for t in taus:
                fi = (S[i] * t) % 1
                fj = (S[j] * t) % 1
                wi = 0 if fi < F(1, 2) else 1
                wj = 0 if fj < F(1, 2) else 1
                # signed comparison of half-plane membership; if equal use frac compare
                if wi != wj:
                    score += (1 if wi < wj else -1)
                else:
                    score += (1 if fi < fj else -1) if fi != fj else 0
            if score > 0: adj[i][j] = 1
            elif score < 0: adj[j][i] = 1
            else:
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

# ---- METHOD 4: Orbit "lower-section more often" with strict section count ----
# Arc i->j iff runner i lands in the LOWER section (residue closer to 0, i.e. its
# normalized distance ||v_i a / mod|| ... ) more often across the orbit than j.
# We use the SECTION value r=v*a mod mod and count how often i's section is in the
# "lower half" {1..(mod-1)/2}-ish, comparing i vs j by who's deeper. Different
# aggregation than method 1 (uses absolute section band, not pairwise compare).
def method4_band_majority(S, mod):
    m = len(S)
    U = units_mod(mod)
    half = mod / 2
    # "depth" of a section r = min(r, mod-r) = how lonely-far this runner is from observer
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            wi = wj = 0
            for a in U:
                ri = (S[i] * a) % mod
                rj = (S[j] * a) % mod
                di = min(ri, mod - ri)
                dj = min(rj, mod - rj)
                if di > dj: wi += 1   # i is farther from observer (more "lonely-supporting")
                elif dj > di: wj += 1
            if wi > wj: adj[i][j] = 1
            elif wj > wi: adj[j][i] = 1
            else:
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

# ---- METHOD 5: Pairwise binding-crossing parity tournament ----
# THM-524: M attained at binding-pair crossing tau=k/(v_a +/- v_b). For each pair
# (i,j), define the FIRST crossing time tau_ij where ||v_i tau||=||v_j tau|| with both
# binding. Build arc i->j by: at the lonely optimum tau*, who is "ahead" in the
# direction of motion of the binding configuration. To avoid transitivity, aggregate
# the SIGN of (v_i - v_j) against the parity of floor(v_i tau*)+floor(v_j tau*).
# Vertex set = runners. This mixes a global tau* with per-pair parity (non-snapshot).
def method5_crossing_parity(S, mod=None):
    taus, gap = all_opt_taus(S)
    m = len(S)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            score = 0
            for t in taus:
                # number of half-integer crossings each runner has made by time t:
                pi = int(2 * S[i] * t)   # floor(2 v_i t) = # of half-turns
                pj = int(2 * S[j] * t)
                # parity-weighted comparison
                par = (pi + pj) % 2
                base = 1 if S[i] > S[j] else -1
                score += base if par == 0 else -base
            if score > 0: adj[i][j] = 1
            elif score < 0: adj[j][i] = 1
            else:
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

# ---- METHOD 6: Section-vertex tournament (vertices = SECTIONS, not runners) ----
# Vertices = the occupied nonzero sections at grid times. For a fixed speed set S and a
# unit a, runner i sits in section v_i*a mod mod. Build a tournament on the SECTIONS
# s in {1..mod-1}: arc s->s' iff, aggregated over the orbit and over S's runners, a
# runner more often "passes from s to s'" as a increments through the units (a dynamic
# transition majority). This is a genuinely different object (mod-1 vertices).
def method6_section_transition(SfamSets, mod):
    """Aggregate transitions over a FAMILY of speed sets to populate section graph.
    Returns adjacency on vertices = nonzero sections 1..mod-1."""
    V = mod - 1
    # map section value -> index 0..V-1
    idx = {s: s - 1 for s in range(1, mod)}
    U = units_mod(mod)
    wins = [[0] * V for _ in range(V)]
    for S in SfamSets:
        for v in S:
            secs = []
            for a in U:
                r = (v * a) % mod
                if r != 0:
                    secs.append(r)
            # consecutive transitions in the orbit-ordering of a
            for k in range(len(secs) - 1):
                s, sp = secs[k], secs[k + 1]
                if s != sp:
                    wins[idx[s]][idx[sp]] += 1
    adj = [[0] * V for _ in range(V)]
    for i in range(V):
        for j in range(i + 1, V):
            if wins[i][j] > wins[j][i]: adj[i][j] = 1
            elif wins[j][i] > wins[i][j]: adj[j][i] = 1
            else: adj[i][j] = 1  # tie -> arbitrary
    return adj

# ============================================================
# DRIVER
# ============================================================
def study_runner_method(name, fn, needs_prime, m_values, maxspeed, mod_or_prime,
                        constrained=True):
    """For a runner-vertex method, enumerate realized iso classes over LRC inputs."""
    print("=" * 70)
    print(f"METHOD: {name}   (constrained-LRC={constrained})")
    results = {}
    for m in m_values:
        free = all_iso_classes(m)
        realized = {}
        nontriv = False
        sets = gen_speed_sets(m, maxspeed)
        n_used = 0
        n_skipped = 0
        for S in sets:
            if constrained:
                # constrain to LRC-relevant: keep all primitive sets (already filtered);
                # this is the family we range over.
                pass
            mod = mod_or_prime
            adj = fn(S, mod)
            if adj is None:
                n_skipped += 1
                continue
            if not is_valid_tournament(adj, m):
                n_skipped += 1
                continue
            n_used += 1
            k = canon_key(adj, m)
            h = h_count(adj, m)
            if h > 1: nontriv = True
            realized.setdefault(k, [0, free.get(k, ("?",))[0] if k in free else None])
            realized[k][0] += 1
        forbidden = set(free.keys()) - set(realized.keys())
        results[m] = {
            "free_count": len(free),
            "realized_count": len(realized),
            "forbidden_count": len(forbidden),
            "nontriv": nontriv,
            "n_used": n_used,
            "n_skipped": n_skipped,
            "forbidden_details": [(free[k][0], free[k][1], free[k][2]) for k in forbidden],
            "realized_details": [(free[k][0] if k in free else "?",
                                  free[k][1] if k in free else "?",
                                  free[k][2] if k in free else "?", cnt)
                                 for k, (cnt, _) in realized.items()],
        }
        print(f"  n={m}: free={len(free)}  realized={len(realized)}  "
              f"forbidden={len(forbidden)}  nontrivial={nontriv}  "
              f"(used {n_used} sets, skipped {n_skipped})")
        if forbidden:
            print(f"     FORBIDDEN classes (H, #3cyc, score): "
                  f"{sorted([(free[k][0], free[k][1], free[k][2]) for k in forbidden])}")
        # show H distribution of realized
        hdist = {}
        for k in realized:
            h = free[k][0] if k in free else h_count_from_key(k, m)
            hdist[h] = hdist.get(h, 0) + 1
        print(f"     realized H-class distribution: {dict(sorted(hdist.items()))}")
    return results

def h_count_from_key(k, m):
    # reconstruct adj from canon key to count H (key is flattened off-diagonal bits)
    adj = [[0] * m for _ in range(m)]
    idx = 0
    for i in range(m):
        for j in range(m):
            if i != j:
                adj[i][j] = k[idx]; idx += 1
    return h_count(adj, m)


if __name__ == "__main__":
    print("A000568 free iso-class counts: n=3->2, n=4->4, n=5->12\n")
    for m in [3, 4, 5]:
        print(f"  free({m}) = {len(all_iso_classes(m))}")
    print()

    # METHOD 1: orbit majority, mod 14
    study_runner_method("M1 orbit-majority-section (mod 14)", method1_orbit_majority,
                        False, [3, 4, 5], 26, 14)

    # METHOD 1 prime variant: mod 13 (prime), to see if group structure matters
    study_runner_method("M1b orbit-majority-section (mod 13 prime)", method1_orbit_majority,
                        False, [3, 4, 5], 26, 13)

    # METHOD 2: QR character, prime p=7 (3 mod 4) and p=11 (3 mod 4)
    study_runner_method("M2 QR-character (p=7)", method2_qr_char, True, [3, 4, 5], 20, 7)
    study_runner_method("M2b QR-character (p=11)", method2_qr_char, True, [3, 4, 5], 24, 11)
    study_runner_method("M2c QR-character (p=23)", method2_qr_char, True, [3, 4, 5], 30, 23)

    # METHOD 3: half-plane majority over optimal taus
    study_runner_method("M3 halfplane-majority (opt taus)", method3_halfplane_majority,
                        False, [3, 4, 5], 16, None)

    # METHOD 4: band majority (depth), mod 14
    study_runner_method("M4 band-depth-majority (mod 14)", method4_band_majority,
                        False, [3, 4, 5], 26, 14)

    # METHOD 5: crossing-parity over optimal taus
    study_runner_method("M5 crossing-parity (opt taus)", method5_crossing_parity,
                        False, [3, 4, 5], 16, None)

    # METHOD 6: section-transition tournament (separate analysis, vertices=sections)
    print("=" * 70)
    print("METHOD 6: section-transition tournament (vertices = sections mod 14)")
    fam = gen_speed_sets(5, 20)
    adj6 = method6_section_transition(fam, 14)
    V = 13
    print(f"  V={V} sections. H={h_count(adj6, V) if V<=8 else 'too big'}  "
          f"#3cyc={num_3cycles(adj6, V)}  score={score_seq(adj6, V)}")
    # smaller mod for iso enumeration
    print("  (mod 6 variant, V=5 for iso enumeration)")
    fam6 = gen_speed_sets(4, 12)
    adj6b = method6_section_transition(fam6, 6)
    V6 = 5
    print(f"  V={V6}: H={h_count(adj6b,V6)}  #3cyc={num_3cycles(adj6b,V6)}  score={score_seq(adj6b,V6)}")
