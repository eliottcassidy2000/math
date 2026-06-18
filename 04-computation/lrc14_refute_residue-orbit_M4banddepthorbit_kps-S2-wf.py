"""
lrc14_refute_residue-orbit_M4banddepthorbit_kps-S2-wf.py

ADVERSARIAL REFUTATION of the M4 band-depth orbit-majority forbidden-class claim.

CLAIM under attack (theme residue-orbit, METHOD M4):
  Map: vertices = runners (n=5 -> 5 vertices). For each unit a in (Z/14)* = {1,3,5,9,11,13},
       section r_i(a)=v_i*a mod 14, depth d_i(a)=min(r_i(a),14-r_i(a)) in {0..7}.
       Arc i->j iff #{a: d_i(a)>d_j(a)} > #{a: d_j(a)>d_i(a)} (6-unit depth-majority).
       Ties (win=0) broken by raw speed (i->j iff S[i]<S[j], i.e. smaller speed -> arc out;
       we keep EXACT explorer convention from method4_band_majority).
  Of the 12 iso classes on 5 vertices (A000568(5)=12), EXACTLY 4 are claimed PROVED
  UNREACHABLE:
     (H=9,  #3cyc=3, score(1,1,2,3,3))
     (H=13, #3cyc=4, score(1,2,2,2,3))
     (H=15, #3cyc=4, score(1,2,2,2,3))
     (H=15, #3cyc=5, score(2,2,2,2,2))   <-- the REGULAR tournament (Paley T_5)
  Reason given: M4 depends ONLY on residues mod 14, so residue-exhaustive => airtight.

STRATEGY:
  M4 is a function of the residue vector (v_i mod 14) ONLY. So the FULL set of iso
  classes realizable by M4 over ANY runner-vertex assignment equals the set realizable
  over all 5-tuples of residues mod 14 (with repetition allowed: two runners CAN share a
  residue, since SDR is strictly stronger than loneliness; loneliness only needs section 0
  empty). We therefore enumerate:
   (A) ALL residue 5-tuples mod 14 (14^5 with repetition) -> the absolute ceiling of M4.
   (B) restrict to LRC-LONELY residue patterns (no runner in section 0 at the witnessing
       grid time) and to ACTUAL primitive integer speed sets over a broad window, including
       covering / sporadic / tight cases.
   (C) off-grid does NOT change M4 (M4 is defined on the grid orbit of units), but we still
       feed genuinely LRC-feasible (primitive, gap>=1/14 if we want, or just primitive)
       speed sets across a large window to confirm.

  If ANY residue 5-tuple (even with collisions) yields one of the 4 forbidden classes,
  the "PROVED unreachable" claim is REFUTED at the residue level. If additionally that
  residue pattern is realizable by genuine distinct integer speeds, it is refuted at the
  speed level too.

All exact integer arithmetic. No floats in decisions.
"""

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

MOD = 14
UNITS = [a for a in range(1, MOD) if gcd(a, MOD) == 1]  # {1,3,5,9,11,13}

# ---------------- EXACT M TOOL ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t):
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
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

# ---------------- M4 MAP (verbatim semantics from method4_band_majority, mod 14) ----------------
def depth(r):
    return min(r % MOD, (MOD - (r % MOD)) % MOD)  # min(r,14-r), with r=0 -> 0

def m4_adj_from_residues(res):
    """res: list of m residues mod 14 (one per runner; repetition allowed).
       Returns adjacency m x m. Tie-break by raw INDEX-tagged speed is supplied separately
       when needed; here for pure-residue test we tie-break by residue value as a stand-in
       only when speeds unavailable -- BUT to match the exact map we REQUIRE speeds for
       tie-break. This function takes an optional speeds list for tie-break."""
    raise NotImplementedError

def m4_adj(S):
    """EXACT M4: input integer speeds S (length m). Uses residues mod 14 and raw-speed
       tie-break exactly as method4_band_majority."""
    m = len(S)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            wi = wj = 0
            for a in UNITS:
                di = depth(S[i]*a)
                dj = depth(S[j]*a)
                if di > dj: wi += 1
                elif dj > di: wj += 1
            if wi > wj: adj[i][j] = 1
            elif wj > wi: adj[j][i] = 1
            else:
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

def m4_adj_residues_tb(res, tb):
    """M4 from residues `res` with explicit tie-break key list `tb` (smaller tb -> arc out).
       This lets us explore ALL tie-break orderings, not just raw-speed."""
    m = len(res)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            wi = wj = 0
            for a in UNITS:
                di = depth(res[i]*a)
                dj = depth(res[j]*a)
                if di > dj: wi += 1
                elif dj > di: wj += 1
            if wi > wj: adj[i][j] = 1
            elif wj > wi: adj[j][i] = 1
            else:
                if tb[i] < tb[j]: adj[i][j] = 1
                elif tb[j] < tb[i]: adj[j][i] = 1
                else: adj[i][j] = 1  # equal tb -> arbitrary; will be canonicalized over both
    return adj

# ---------------- TOURNAMENT INVARIANTS ----------------
def canon_key(adj, m):
    best = None
    for perm in permutations(range(m)):
        bits = []
        for i in range(m):
            for j in range(m):
                if i != j: bits.append(adj[perm[i]][perm[j]])
        t = tuple(bits)
        if best is None or t < best: best = t
    return best
def h_count(adj, m):
    cnt = 0
    for perm in permutations(range(m)):
        ok = True
        for k in range(m-1):
            if adj[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok: cnt += 1
    return cnt
def num_3cycles(adj, m):
    c = 0
    for a, b, cc in combinations(range(m), 3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def score_seq(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def is_valid(adj, m):
    for i in range(m):
        for j in range(i+1, m):
            if adj[i][j] + adj[j][i] != 1: return False
    return True

# ---------------- FREE SET on 5 vertices ----------------
def all_iso_classes(m):
    seen = {}
    pairs = list(combinations(range(m), 2))
    for bits in product([0,1], repeat=len(pairs)):
        adj = [[0]*m for _ in range(m)]
        for (i,j), b in zip(pairs, bits):
            if b: adj[i][j] = 1
            else: adj[j][i] = 1
        k = canon_key(adj, m)
        if k not in seen:
            seen[k] = (h_count(adj, m), num_3cycles(adj, m), score_seq(adj, m))
    return seen

FORBIDDEN = {
    (9, 3, (1,1,2,3,3)),
    (13,4, (1,2,2,2,3)),
    (15,4, (1,2,2,2,3)),
    (15,5, (2,2,2,2,2)),
}

def signature(adj, m):
    return (h_count(adj,m), num_3cycles(adj,m), score_seq(adj,m))

# ================================================================
def main():
    m = 5
    free = all_iso_classes(m)
    print(f"FREE set on {m} vertices: {len(free)} iso classes (A000568(5)=12)")
    free_sigs = {}
    for k,(h,c,s) in free.items():
        free_sigs.setdefault((h,c,s), set()).add(k)
    print("Forbidden signatures under attack:")
    for f in sorted(FORBIDDEN): print("   ", f)
    print()

    # ----------------------------------------------------------------
    # (A) ABSOLUTE CEILING: ALL residue 5-tuples mod 14 WITH repetition.
    #     Tie-break: residues can collide -> need a tie-break. We test EVERY
    #     tie-break ordering by trying all permutations of distinct tb ranks AND
    #     the degenerate equal case. To be thorough we enumerate residue multisets
    #     and for each we try all consistent strict tie-break orders.
    # ----------------------------------------------------------------
    print("="*72)
    print("(A) ABSOLUTE RESIDUE CEILING: all residue 5-tuples mod 14 (repetition allowed),")
    print("    all tie-break orders. This is the MAXIMUM M4 can ever produce.")
    print("="*72)
    realized_keys = set()
    realized_sigs = {}
    hit_forbidden = []  # witnesses
    # iterate residue vectors with repetition; to cut symmetry, residues as sorted tuple
    # but tie-break order matters, so we still need order info. We enumerate residue
    # MULTISETS and assign distinct tb ranks via all permutations of positions.
    count = 0
    for res in product(range(MOD), repeat=m):
        count += 1
        # use distinct tb = the position index permutations to cover all tie orders.
        # but most pairs won't tie; only colliding-depth pairs tie. We just need ONE
        # representative per residue vector with a generic strict tb, PLUS we sweep tb
        # only when there is a tie. To keep it exact and complete, sweep all m! tb orders
        # ONLY for residue vectors that have >=2 equal depths (potential ties).
        depths = [depth(res[i]) for i in range(m)]
        # quick: does any pair tie in the 6-unit majority? approximate by trying generic
        # tb first; then if any pair had wi==wj, sweep.
        adj = m4_adj_residues_tb(list(res), tb=list(range(m)))
        if not is_valid(adj, m):
            continue
        k = canon_key(adj, m); realized_keys.add(k)
        sig = signature(adj, m); realized_sigs.setdefault(sig, list())
        if len(realized_sigs[sig]) < 1: realized_sigs[sig].append((res, tuple(range(m))))
        if sig in FORBIDDEN:
            hit_forbidden.append(("residue-ceiling", res, tuple(range(m)), sig))
        # detect ties and sweep tb
        tie = False
        for i in range(m):
            for j in range(i+1, m):
                wi=wj=0
                for a in UNITS:
                    di=depth(res[i]*a); dj=depth(res[j]*a)
                    if di>dj: wi+=1
                    elif dj>di: wj+=1
                if wi==wj: tie=True
        if tie:
            for tb in permutations(range(m)):
                adj2 = m4_adj_residues_tb(list(res), tb=list(tb))
                if not is_valid(adj2, m): continue
                k2 = canon_key(adj2, m); realized_keys.add(k2)
                sig2 = signature(adj2, m)
                realized_sigs.setdefault(sig2, list())
                if sig2 in FORBIDDEN:
                    hit_forbidden.append(("residue-ceiling-tb", res, tb, sig2))
    print(f"  enumerated {count} residue vectors (14^5={MOD**m})")
    print(f"  M4 realizes {len(realized_keys)} iso KEYS over residue ceiling")
    real_sig_set = set(realized_sigs.keys())
    print(f"  realized SIGNATURES: {len(real_sig_set)}")
    forb_hit = [s for s in real_sig_set if s in FORBIDDEN]
    print(f"  forbidden signatures HIT at residue ceiling: {sorted(forb_hit)}")
    # which of the 12 free are reachable / not
    reach_sigs = set()
    for sig in real_sig_set:
        reach_sigs.add(sig)
    all_free_sigs = set(free_sigs.keys())
    unreached = all_free_sigs - reach_sigs
    print(f"  free signatures NEVER reached even at ceiling: {sorted(unreached)}")
    print(f"  (count reached={len(reach_sigs & all_free_sigs)} / {len(all_free_sigs)})")
    print()
    if hit_forbidden:
        print("  *** FORBIDDEN HITS (REFUTATION at residue level) ***")
        for w in hit_forbidden[:10]:
            print("    ", w)
    else:
        print("  No forbidden signature reachable at the residue ceiling.")
    print()

    # ----------------------------------------------------------------
    # (B) GENUINE LRC SPEED SETS, broad window. Distinct positive primitive speeds.
    #     This is the honest physical test. We sweep speeds 1..VMAX, include covering
    #     and sporadic configs. Record realized signatures and any forbidden hit.
    # ----------------------------------------------------------------
    print("="*72)
    print("(B) GENUINE primitive integer speed sets, broad windows.")
    print("="*72)
    def primitive(S):
        g=0
        for v in S: g=gcd(g,v)
        return g==1
    for VMAX in [20, 30, 42]:
        realized = {}
        forb_hits = []
        n_used = 0
        for S in combinations(range(1, VMAX+1), m):
            if not primitive(S): continue
            adj = m4_adj(list(S))
            if not is_valid(adj, m): continue
            n_used += 1
            sig = signature(adj, m)
            realized.setdefault(sig, S)
            if sig in FORBIDDEN:
                forb_hits.append((S, sig))
        forb_seen = [s for s in realized if s in FORBIDDEN]
        print(f"  VMAX={VMAX}: used {n_used} primitive 5-sets; realized {len(realized)} signatures.")
        print(f"     forbidden hit: {sorted(forb_seen) if forb_seen else 'NONE'}")
        if forb_hits:
            for w in forb_hits[:5]: print("       witness:", w)
    print()

    # ----------------------------------------------------------------
    # (C) LRC-FEASIBLE / sporadic / covering style residue patterns explicitly:
    #     residues that include a 0 (parked runner) and patterns with collisions, to
    #     see if loneliness-allowed (section-0-empty but residues colliding) patterns
    #     unlock the forbidden classes. We also try ALL ordered residue tuples with
    #     distinct nonzero residues (perfect SDR) -- the cleanest LRC-lonely case.
    # ----------------------------------------------------------------
    print("="*72)
    print("(C) DISTINCT NONZERO residues (perfect SDR, the canonical LRC-lonely grid case)")
    print("="*72)
    realized = set(); forb = []
    for res in combinations(range(1, MOD), m):  # distinct nonzero residues
        for perm in permutations(res):
            adj = m4_adj_residues_tb(list(perm), tb=list(range(m)))
            if not is_valid(adj, m): continue
            sig = signature(adj, m)
            realized.add(sig)
            if sig in FORBIDDEN:
                forb.append((perm, sig))
    print(f"  distinct-nonzero-residue signatures: {len(realized)}")
    print(f"  forbidden hit: {sorted(set(s for _,s in forb)) if forb else 'NONE'}")
    free_reached = realized & set(free_sigs.keys())
    print(f"  free signatures reached (perfect SDR): {len(free_reached)} / {len(free_sigs)}")
    print(f"  free signatures NOT reached: {sorted(set(free_sigs)-realized)}")
    print()

    # ----------------------------------------------------------------
    # FINAL VERDICT
    # ----------------------------------------------------------------
    print("="*72)
    print("VERDICT")
    print("="*72)
    any_hit = bool(hit_forbidden) or any(False for _ in [])
    # gather all forbidden hits across A,B,C
    print(f"  (A) residue ceiling forbidden hits: {len(hit_forbidden)}")
    print(f"  forbidden classes confirmed unreachable across A/B/C: "
          f"{FORBIDDEN - set(forb_hit)}")
    if hit_forbidden:
        print("  RESULT: REFUTED -- at least one forbidden class realized. See witnesses above.")
    else:
        print("  RESULT: CONFIRMED at residue ceiling (14^5 exhaustive incl. collisions &")
        print("          all tie-break orders). The 4 forbidden classes are unreachable by M4.")

if __name__ == "__main__":
    main()
