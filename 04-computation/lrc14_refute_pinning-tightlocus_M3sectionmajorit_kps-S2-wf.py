"""
Adversarial refutation of the 'pinning-tightlocus' forbidden-class claim.

CLAIM (to refute): The M3 section-majority-vote map, applied to LRC-mod-14
loneliness on the grid tau=a/14, a in (Z/14)* = {1,3,5,9,11,13}, NEVER
realizes the regular tournament SCORE SEQUENCE (2,2,2,2,2) at k=5 (the unique
forbidden score), over the genuine primitive-speed family. Also claimed
(weaker) forbidden at k=5: score (1,1,2,3,3) [c3=3,H=9] and two (1,2,2,2,3)
classes.

M3 MAP (exact, copied from the claim):
  Vertices = the n runners (speeds v_i).
  On the grid tau = a/14 for a in U14 = {1,3,5,9,11,13}, runner i sits in
  SECTION  r_i(a) = (v_i * a) mod 14, with DEPTH d_i(a) = min(r, 14-r)
  (distance to section 0; small depth = closer to lonely-violation).
  Arc i->j  iff  margin(i,j) = #{a: d_i(a) > d_j(a)} - #{a: d_i(a) < d_j(a)}
  is POSITIVE (i 'deeper/safer' across more columns).
  Tie (margin == 0) broken by SMALLER speed -> arc points FROM smaller speed
  (i.e. if v_i < v_j and tie, arc i->j). [we make tie-break deterministic on speed]

We build the tournament adjacency, compute the SCORE SEQUENCE (out-degrees),
the number of directed 3-cycles c3, and H (#Hamiltonian paths), and check
whether score (2,2,2,2,2) is ever realized.

NOTE on depth: d_i(a) depends ONLY on (v_i mod 14). So the M3 tournament
depends only on the multiset of residues r_i = v_i mod 14. We enumerate over
residue assignments directly for the residue-level question, AND over genuine
primitive speed sets for the honest empirical question.

All exact integer arithmetic (no floats).
"""

from itertools import combinations, permutations, product
from math import gcd
from functools import reduce

U14 = [1, 3, 5, 9, 11, 13]  # (Z/14)*

def depth(r):
    """depth of residue r mod 14 = min(r, 14-r), r in 0..13."""
    r = r % 14
    return min(r, 14 - r)

def depth_profile(v):
    """6-vector of depths d(v*a mod 14) for a in U14. Depends only on v mod 14."""
    return tuple(depth((v * a) % 14) for a in U14)

def m3_tournament(speeds):
    """
    Build M3 adjacency. speeds: list of n distinct positive ints (the runners).
    Returns adjacency dict adj[i][j] = True if arc i->j.
    Tie-break: smaller speed wins the tie (arc from smaller -> larger),
    matching 'tie broken by smaller speed' = smaller speed is treated as
    'deeper/safer' default. We make it deterministic & antisymmetric.
    """
    n = len(speeds)
    prof = [depth_profile(speeds[i]) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            gt = sum(1 for k in range(6) if prof[i][k] > prof[j][k])
            lt = sum(1 for k in range(6) if prof[i][k] < prof[j][k])
            margin = gt - lt  # margin(i,j)
            if margin > 0:
                adj[i][j] = True            # i -> j
            elif margin < 0:
                adj[j][i] = True            # j -> i
            else:
                # tie: smaller speed -> arc from smaller speed to larger
                if speeds[i] < speeds[j]:
                    adj[i][j] = True
                else:
                    adj[j][i] = True
    return adj

def scores(adj):
    n = len(adj)
    return tuple(sorted(sum(1 for j in range(n) if adj[i][j]) for i in range(n)))

def count_c3(adj):
    n = len(adj)
    c = 0
    for a, b, cc in combinations(range(n), 3):
        # count directed 3-cycles among the triple
        # cycle a->b->c->a or a->c->b->a
        if adj[a][b] and adj[b][cc] and adj[cc][a]:
            c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]:
            c += 1
    return c

def count_H(adj):
    """Number of Hamiltonian paths (Redei). n small."""
    n = len(adj)
    cnt = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if not adj[perm[i]][perm[i + 1]]:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt

def is_primitive(S):
    return reduce(gcd, S) == 1

def lonely_on_grid(speeds):
    """True if exists a in U14 with no runner in section 0 (observer lonely at a/14)."""
    for a in U14:
        if all((v * a) % 14 != 0 for v in speeds):
            return True
    return False

# ---------------------------------------------------------------------------
# PART 0: sanity / baseline
# ---------------------------------------------------------------------------
def part0():
    print("=" * 70)
    print("PART 0: baseline sanity")
    S = [1, 2, 3, 4, 5]
    adj = m3_tournament(S)
    sc = scores(adj)
    print(f"S=(1,2,3,4,5): score={sc}, c3={count_c3(adj)}, H={count_H(adj)}")
    # print depth profiles & tier structure
    print("depth profiles for residues 0..13:")
    for r in range(14):
        dp = depth_profile(r)
        print(f"  res {r:2d}: profile {dp}  sum={sum(dp)}")

# ---------------------------------------------------------------------------
# PART 1: RESIDUE-LEVEL exhaustive search at k=5.
# The M3 tournament depends ONLY on residues mod 14. Enumerate all
# residue 5-tuples (with distinct residues -> 'clean' SDR-like, and also
# with repeats) and check realized scores. Also require the residue set to
# come from a LONELY (loneliness on grid) config to stay LRC-meaningful.
# ---------------------------------------------------------------------------
def part1_residue_clean():
    print("=" * 70)
    print("PART 1a: residue-level, DISTINCT residues 1..13 (clean SDR), k=5")
    print("(distinct nonzero residues -> tie-break rarely needed; honest residue test)")
    realized_scores = {}
    target = (2, 2, 2, 2, 2)
    found = []
    for combo in combinations(range(1, 14), 5):  # distinct nonzero residues
        # treat residues as 'speeds' for M3 (depth depends only on residue);
        # but tie-break uses speed ordering = residue ordering here.
        adj = m3_tournament(list(combo))
        sc = scores(adj)
        realized_scores[sc] = realized_scores.get(sc, 0) + 1
        if sc == target:
            found.append(combo)
    print(f"  total distinct-residue 5-sets: {sum(realized_scores.values())}")
    print(f"  realized score sequences: {sorted(realized_scores.keys())}")
    print(f"  score (2,2,2,2,2) realized? {len(found)>0}  count={len(found)}")
    if found:
        print(f"  WITNESS residues: {found[:5]}")
    # secondary claimed-forbidden scores
    for s in [(1,1,2,3,3),(1,2,2,2,3)]:
        print(f"  score {s} realized? {s in realized_scores} count={realized_scores.get(s,0)}")
    return found, realized_scores

def part1_residue_all():
    print("=" * 70)
    print("PART 1b: residue-level, ALL residue 5-tuples 0..13 WITH repeats allowed")
    print("(including section-0 residues; tie-break by an injected speed index)")
    # We must give M3 actual speeds for tie-break. To explore tie-break freedom,
    # we test residue multisets and, for each, try to realize (2,2,2,2,2) by
    # choosing a speed lift that ASSIGNS the residues to speeds with various
    # orderings. We test by enumerating residue tuples and using the index
    # order as the tie-break speed proxy AND the reverse, to probe tie-break.
    target = (2, 2, 2, 2, 2)
    found = []
    realized = set()
    count = 0
    for combo in product(range(14), repeat=5):
        # require it to be a 'genuine' multiset that is lonely on grid if we
        # think of these as residues (loneliness uses residue==0 detection).
        # To stay LRC-meaningful, require lonely_on_grid via residues:
        if not any(all(r != 0 for r in [(rr * a) % 14 for rr in combo]) for a in U14):
            # this residue-set is a forced-non-lonely (covering) config on grid;
            # still LRC-relevant (off-grid lonely) so we DON'T skip, but mark.
            pass
        count += 1
        # speeds = residues + 14*index to make them distinct & preserve residue,
        # tie-break order = the index order.
        speeds = [combo[i] + 14 * i if combo[i] != 0 else 14 * (i + 1) for i in range(5)]
        # ensure residues preserved:
        adj = m3_tournament(speeds)
        sc = scores(adj)
        realized.add(sc)
        if sc == target:
            found.append((combo, speeds))
    print(f"  total residue 5-tuples (with repeats): {count}")
    print(f"  distinct realized scores: {len(realized)}")
    print(f"  score (2,2,2,2,2) realized? {len(found)>0} count={len(found)}")
    if found:
        for c, s in found[:8]:
            adj = m3_tournament(s)
            print(f"    WITNESS residues={c} speeds={s} c3={count_c3(adj)} H={count_H(adj)}")
    return found, realized

def part1_residue_tiebreak_freedom():
    print("=" * 70)
    print("PART 1c: residue multiset WITH ALL tie-break orderings")
    print("(probe whether (2,2,2,2,2) is reachable via tie-break choice on repeats)")
    # For residue multisets that have ties (margin 0 pairs), enumerate ALL
    # antisymmetric tie-break orientations of the tied pairs and see if any
    # gives (2,2,2,2,2). This is the MOST GENEROUS reading: 'free tie-breaks'.
    target = (2, 2, 2, 2, 2)
    realized_with_free_tiebreak = set()
    witnesses = []
    n = 5
    # enumerate residue multisets (sorted to dedupe)
    seen_multisets = set()
    total = 0
    for combo in combinations_with_replacement(range(14), 5):
        if combo in seen_multisets:
            continue
        seen_multisets.add(combo)
        prof = [depth_profile(r) for r in combo]
        # compute fixed (non-tie) arcs and list of tied pairs
        fixed = [[None] * n for _ in range(n)]
        tied = []
        for i in range(n):
            for j in range(i + 1, n):
                gt = sum(1 for k in range(6) if prof[i][k] > prof[j][k])
                lt = sum(1 for k in range(6) if prof[i][k] < prof[j][k])
                mar = gt - lt
                if mar > 0:
                    fixed[i][j] = True; fixed[j][i] = False
                elif mar < 0:
                    fixed[i][j] = False; fixed[j][i] = True
                else:
                    tied.append((i, j))
        if len(tied) > 12:
            # too many tie orientations to brute (2^12=4096 ok, cap higher)
            pass
        total += 1
        for bits in range(1 << len(tied)):
            adj = [[fixed[i][j] for j in range(n)] for i in range(n)]
            for idx, (i, j) in enumerate(tied):
                if (bits >> idx) & 1:
                    adj[i][j] = True; adj[j][i] = False
                else:
                    adj[i][j] = False; adj[j][i] = True
            sc = scores(adj)
            realized_with_free_tiebreak.add(sc)
            if sc == target:
                witnesses.append((combo, tied, bits))
                break  # one orientation suffices to mark this multiset reachable
    print(f"  residue multisets tested: {total}")
    print(f"  (2,2,2,2,2) reachable with SOME tie-break? {len(witnesses)>0}")
    print(f"  number of multisets that can reach it: {len(witnesses)}")
    if witnesses:
        for combo, tied, bits in witnesses[:8]:
            print(f"    multiset {combo}, #tied-pairs={len(tied)}")
    print(f"  distinct scores reachable (free tie-break): {sorted(realized_with_free_tiebreak)}")
    return witnesses

from itertools import combinations_with_replacement

# ---------------------------------------------------------------------------
# PART 2: GENUINE PRIMITIVE SPEED SETS (honest empirical), broad search.
# Larger vmax than the claim's 40, plus targeted covering/tight/sporadic sets.
# ---------------------------------------------------------------------------
def part2_primitive(vmax, require_lonely=False):
    print("=" * 70)
    print(f"PART 2: genuine primitive 5-speed sets, vmax={vmax}, require_lonely={require_lonely}")
    target = (2, 2, 2, 2, 2)
    realized = {}
    found = []
    total = 0
    for S in combinations(range(1, vmax + 1), 5):
        if not is_primitive(S):
            continue
        if require_lonely and not lonely_on_grid(S):
            continue
        total += 1
        adj = m3_tournament(list(S))
        sc = scores(adj)
        realized[sc] = realized.get(sc, 0) + 1
        if sc == target:
            found.append(S)
            if len(found) <= 10:
                print(f"  *** REGULAR WITNESS: {S} c3={count_c3(adj)} H={count_H(adj)}")
    print(f"  primitive sets tested: {total}")
    print(f"  (2,2,2,2,2) realized? {len(found)>0} count={len(found)}")
    print(f"  realized scores: {sorted(realized.keys())}")
    for s in [(1,1,2,3,3),(1,2,2,2,3)]:
        print(f"  score {s}: count={realized.get(s,0)}")
    return found, realized

# ---------------------------------------------------------------------------
# PART 3: targeted hard / covering / tight 5-subsets and residue patterns
# including the residue-7 ('seven' tier) and section-0 (multiple of 14).
# ---------------------------------------------------------------------------
def part3_targeted():
    print("=" * 70)
    print("PART 3: targeted families (covering, residue-7, section-0, large spread)")
    target = (2, 2, 2, 2, 2)
    fams = []
    # include speeds that are multiples of 14 (section 0 forever on grid)
    fams.append([7, 14, 21, 28, 1])     # multiples of 7/14
    fams.append([14, 28, 42, 56, 1])
    fams.append([7, 21, 35, 49, 1])     # all residue 7
    fams.append([1, 13, 3, 11, 5])      # all odd-tier residues
    fams.append([2, 4, 6, 8, 10])       # all even-tier
    fams.append([7, 2, 1, 4, 3])
    # large spread random-ish deterministic primitive sets
    import random
    random.seed(1234)
    big = []
    for _ in range(200000):
        S = tuple(sorted(random.sample(range(1, 2000), 5)))
        if is_primitive(S):
            big.append(S)
    # dedupe
    big = list(set(big))
    realized = {}
    found = []
    for S in fams:
        adj = m3_tournament(list(S))
        sc = scores(adj)
        print(f"  family {S}: score={sc} c3={count_c3(adj)} H={count_H(adj)}")
        if sc == target:
            found.append(tuple(S))
    for S in big:
        adj = m3_tournament(list(S))
        sc = scores(adj)
        realized[sc] = realized.get(sc, 0) + 1
        if sc == target:
            found.append(S)
            if len([f for f in found]) <= 10:
                print(f"  *** REGULAR WITNESS (large): {S} H={count_H(adj)}")
    print(f"  large random primitive sets tested: {len(big)}")
    print(f"  (2,2,2,2,2) realized in large search? count={sum(1 for s in found if max(s)>=14)}")
    print(f"  realized scores (large): {sorted(realized.keys())}")
    return found

# ---------------------------------------------------------------------------
# PART 4: N-dependence control (claim says regular IS realized at N=15,20).
# Reproduce with general modulus N to validate our map matches the claim.
# ---------------------------------------------------------------------------
def units_mod(N):
    return [a for a in range(1, N) if gcd(a, N) == 1]

def depth_profile_N(v, N, units):
    return tuple(min((v * a) % N, N - ((v * a) % N)) for a in units)

def m3_tournament_N(speeds, N):
    units = units_mod(N)
    n = len(speeds)
    prof = [depth_profile_N(speeds[i], N, units) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    L = len(units)
    for i in range(n):
        for j in range(i + 1, n):
            gt = sum(1 for k in range(L) if prof[i][k] > prof[j][k])
            lt = sum(1 for k in range(L) if prof[i][k] < prof[j][k])
            mar = gt - lt
            if mar > 0:
                adj[i][j] = True
            elif mar < 0:
                adj[j][i] = True
            else:
                if speeds[i] < speeds[j]:
                    adj[i][j] = True
                else:
                    adj[j][i] = True
    return adj

def part4_Ncontrol():
    print("=" * 70)
    print("PART 4: N-dependence control. Does regular (2,2,2,2,2) appear at N!=14?")
    target = (2, 2, 2, 2, 2)
    for N in [13, 14, 15, 20, 28, 22]:
        found = 0
        total = 0
        for S in combinations(range(1, 26), 5):
            if not is_primitive(S):
                continue
            total += 1
            adj = m3_tournament_N(list(S), N)
            if scores(adj) == target:
                found += 1
        print(f"  N={N}: regular realized count={found} / {total} primitive sets (vmax=25)")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    part0()
    f1a, r1a = part1_residue_clean()
    f1b, r1b = part1_residue_all()
    f1c = part1_residue_tiebreak_freedom()
    f2_25, _ = part2_primitive(25)
    f2_40, _ = part2_primitive(40)
    f2_60, _ = part2_primitive(60)
    f3 = part3_targeted()
    part4_Ncontrol()

    print("=" * 70)
    print("SUMMARY")
    print(f"  PART1a (distinct residues): regular found = {len(f1a)}")
    print(f"  PART1b (residues w/ repeats, index tie-break): regular found = {len(f1b)}")
    print(f"  PART1c (free tie-break): regular reachable multisets = {len(f1c)}")
    print(f"  PART2 vmax25: {len(f2_25)}, vmax40: {len(f2_40)}, vmax60: {len(f2_60)}")
    print(f"  PART3 targeted/large: {len(f3)}")
