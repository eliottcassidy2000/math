"""
Adversarial refutation of the M3b "section-parity character vote" forbidden-class claim.

CLAIM (theme character-spectral, method M3b):
 Map: vertices = n runners with speeds S (positive ints).
 Weight w(a) = (-1)^{section of v0 at a/14} = +1 if v0*a mod 14 is EVEN else -1,
   for a in (Z/14)* = {1,3,5,9,11,13}.  v0 = reference runner (smallest? first? -- see below).
 d0(v,a) = min(v*a mod 14, 14 - (v*a mod 14))  -- cyclic section-distance from 0.
 Arc x -> y iff  sum_a w(a)*sign(d0(x,a) - d0(y,a)) > 0 ; tie-break x>y (by speed).
 CLAIMED FORBIDDEN at n=5 (exhaustive over primitive S subset {1..13}, |S|=5):
   the iso classes:
     (9,3,(1,1,2,3,3))   -- (H, c3, score-sequence)
     (13,4,(1,2,2,2,3))
     (15,4,(1,2,2,2,3))
     regular (15,5,(2,2,2,2,2))
 i.e. these tournament iso classes are NEVER produced by the map over the search domain.

 We adversarially try to PRODUCE any of these forbidden classes with the EXACT map
 over a far broader set of LRC-constrained inputs.

We identify a tournament iso class by a canonical invariant tuple. To be safe we use
the FULL canonical form (canonical adjacency under vertex relabeling) for small n,
and also report (H, c3, sorted score sequence) which is what the claim lists.

Reference v0: the claim says "Reference v0". Ambiguous. We test BOTH conventions:
  v0 = min(S)  and  v0 = the runner listed first / the largest. We test several.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
import sys, random, itertools

UNITS = [1,3,5,9,11,13]  # (Z/14)*

def sec(v,a):  # section = v*a mod 14
    return (v*a) % 14

def d0(v,a):
    s = (v*a) % 14
    return min(s, 14 - s)

def weight_for_ref(v0):
    # w(a) = +1 if v0*a mod 14 even else -1
    return {a: (1 if (v0*a) % 14 % 2 == 0 else -1) for a in UNITS}

def build_tournament(S, ref):
    """S = list of distinct positive ints (the speeds = vertices). ref = reference speed v0.
    Returns adjacency matrix adj[i][j]=1 if i->j. Order vertices by S sorted ascending.
    Tie-break x>y means: if vote sum == 0, larger speed -> smaller? We use 'x>y' as in
    'arc x->y iff vote>0; tie-break x>y' meaning when tie, the larger speed beats smaller,
    i.e. x->y when x's speed > y's speed.
    """
    V = sorted(S)
    n = len(V)
    w = weight_for_ref(ref)
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            x = V[i]; y = V[j]
            vote = 0
            for a in UNITS:
                dx = d0(x,a); dy = d0(y,a)
                if dx > dy: vote += w[a]
                elif dx < dy: vote -= w[a]
            if vote > 0:
                adj[i][j] = 1
            elif vote < 0:
                adj[i][j] = 0
            else:
                # tie: x>y means larger speed wins. x speed vs y speed.
                adj[i][j] = 1 if x > y else 0
    # sanity: tournament (exactly one of adj[i][j],adj[j][i])
    for i in range(n):
        for j in range(i+1,n):
            assert adj[i][j] + adj[j][i] == 1, (S,ref,i,j,adj[i][j],adj[j][i])
    return adj, V

def score_seq(adj):
    n=len(adj)
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def count_3cycles(adj):
    n=len(adj); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not adj[i][j]: continue
            for k in range(n):
                if k in (i,j): continue
                if adj[j][k] and adj[k][i]: c+=1
    return c//3

def ham_paths(adj):
    """count Hamiltonian paths (directed) = H, Redei (always odd)."""
    n=len(adj)
    # DP over subsets
    from functools import lru_cache
    full=(1<<n)-1
    # dp[mask][last] = number of paths covering mask ending at last
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            cur=dp[mask][last]
            if not cur: continue
            for nxt in range(n):
                if mask&(1<<nxt): continue
                if adj[last][nxt]:
                    dp[mask|(1<<nxt)][nxt]+=cur
    return sum(dp[full][v] for v in range(n))

def canon(adj):
    """canonical form under vertex relabeling: min over permutations of the
    upper-triangular bit signature. Only for small n (<=7)."""
    n=len(adj)
    best=None
    for p in permutations(range(n)):
        bits=0
        b=0
        for i in range(n):
            for j in range(n):
                if i!=j:
                    bits = (bits<<1) | adj[p[i]][p[j]]
        if best is None or bits<best:
            best=bits
    return best

def iso_key(adj):
    return (ham_paths(adj), count_3cycles(adj), score_seq(adj), canon(adj))

def short_key(adj):
    return (ham_paths(adj), count_3cycles(adj), score_seq(adj))

# ---- Forbidden targets (from claim) given as (H, c3, score) ----
FORBIDDEN_SHORT = {
    (9,3,(1,1,2,3,3)),
    (13,4,(1,2,2,2,3)),
    (15,4,(1,2,2,2,3)),
    (15,5,(2,2,2,2,2)),
}

def is_primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

def main():
    print("=== M3b section-parity character vote: adversarial refutation ===")
    print("UNITS (Z/14)* =", UNITS)
    # Reference conventions to test
    REF_MODES = ["min","max","first_listed","each"]

    hits = {f:[] for f in FORBIDDEN_SHORT}
    seen_short = set()

    def process(S, ref, src):
        adj,V = build_tournament(S, ref)
        sk = short_key(adj)
        seen_short.add(sk)
        if sk in FORBIDDEN_SHORT:
            hits[sk].append((tuple(V), ref, src))

    # === PHASE 1: exhaustive n=5 over primitive subsets of {1..13} (the claim's domain) ===
    print("\n[Phase 1] n=5 exhaustive primitive subsets of {1..13}, ref in {min,max}...")
    cnt=0
    for S in combinations(range(1,14),5):
        if not is_primitive(S): continue
        cnt+=1
        for ref in (min(S), max(S)):
            process(S, ref, "ph1")
    print(f"  primitive 5-subsets of 1..13: {cnt}")
    print("  short-keys seen so far:", len(seen_short))
    for f in FORBIDDEN_SHORT:
        print(f"   forbidden {f}: {'HIT '+str(hits[f][:3]) if hits[f] else 'not produced'}")

    # === PHASE 2: broaden speed range. n=5, primitive subsets of {1..30}, multiple refs ===
    print("\n[Phase 2] n=5 primitive subsets of {1..30}, ref in {min,max,median}...")
    cnt=0
    for S in combinations(range(1,31),5):
        if not is_primitive(S): continue
        cnt+=1
        Ss=sorted(S)
        refs = {Ss[0], Ss[-1], Ss[2]}
        # also test ref = any speed congruent to convenient residues
        for ref in refs:
            process(S, ref, "ph2")
    print(f"  primitive 5-subsets of 1..30: {cnt}")
    for f in FORBIDDEN_SHORT:
        print(f"   forbidden {f}: {'HIT '+str(hits[f][:3]) if hits[f] else 'not produced'}")

    # === PHASE 3: ref = EACH member (all reference choices), n=5 subsets of {1..40} sampled ===
    print("\n[Phase 3] n=5 random primitive subsets of {1..60}, ref = EVERY member...")
    rng=random.Random(20260617)
    tested=0
    for _ in range(60000):
        S=tuple(sorted(rng.sample(range(1,61),5)))
        if not is_primitive(S): continue
        tested+=1
        for ref in S:
            process(S, ref, "ph3")
    print(f"  random primitive 5-subsets tested: {tested}")
    for f in FORBIDDEN_SHORT:
        print(f"   forbidden {f}: {'HIT '+str(hits[f][:3]) if hits[f] else 'not produced'}")

    # === PHASE 4: residue-targeted. The map only depends on speeds via (v mod 14) actually?
    # d0(v,a) depends on v*a mod 14, i.e. on v mod 14. And vote depends only on residues mod 14.
    # So tournament depends ONLY on the multiset of residues mod 14 (and tie-break on actual sizes).
    # Enumerate ALL residue-vectors mod 14 for n=5 exhaustively (with distinct residues and with repeats).
    print("\n[Phase 4] n=5 EXHAUSTIVE over residues mod 14 (the map depends only on v mod 14)...")
    # residues 0..13; tie-break uses actual speed but if residues equal AND we pick speeds in
    # residue order, tie-break is determined. We test residue vectors; pick representative speeds
    # = residue if nonzero else 14. To exercise tie-break both ways we also try reversed sizes.
    resfound=0
    res_short=set()
    for rv in itertools.combinations_with_replacement(range(0,14),5):
        # build speeds: distinct positive ints with these residues mod 14
        # use r + 14*idx to keep distinct; vary ordering for tie-break coverage
        for orient in (0,1):
            base=[]
            for idx,r in enumerate(rv):
                val = r if r!=0 else 14
                val = val + 14*idx
                base.append(val)
            if orient==1:
                base=base[::-1]
            S=tuple(base)
            if not is_primitive(S):
                # primitivity not required for residue coverage, but keep flag
                pass
            for ref in set(S):
                adj,V=build_tournament(S,ref)
                sk=short_key(adj)
                res_short.add(sk)
                if sk in FORBIDDEN_SHORT:
                    hits[sk].append((tuple(V),ref,"ph4"))
                    resfound+=1
    print(f"  residue-vectors processed; distinct short-keys via residues: {len(res_short)}")
    for f in FORBIDDEN_SHORT:
        print(f"   forbidden {f}: {'HIT '+str(hits[f][:3]) if hits[f] else 'not produced'}")

    # === PHASE 5: LARGER n (the forbidden iso classes are n=5, but a 5-vertex INDUCED
    # subtournament of a larger LRC config could realize them). Take LRC covering/tight/random
    # configs at n up to 13, build the FULL tournament, then check ALL 5-subsets of vertices.
    print("\n[Phase 5] larger LRC configs -> check ALL 5-vertex induced subtournaments...")
    def induced(adj, idxs):
        m=len(idxs)
        sub=[[0]*m for _ in range(m)]
        for a in range(m):
            for b in range(m):
                if a!=b: sub[a][b]=adj[idxs[a]][idxs[b]]
        return sub
    big_configs=[]
    # covering family {1..11,13} + {84m}
    for m in range(1,4):
        big_configs.append(tuple(sorted(list(range(1,12))+[13,84*m])))
    # AP {1..13}
    big_configs.append(tuple(range(1,14)))
    # Goddyn-Wong-ish {1..11,13,24}
    big_configs.append(tuple(sorted(list(range(1,12))+[13,24])))
    # sporadic / random primitive n in 6..13
    for n in range(6,14):
        for _ in range(4000):
            S=tuple(sorted(rng.sample(range(1,80),n)))
            if not is_primitive(S): continue
            big_configs.append(S)
    print(f"  big configs: {len(big_configs)}")
    sub_found=0
    for S in big_configs:
        Vsort=sorted(S)
        for ref in (min(S),max(S)):
            adj,V=build_tournament(S,ref)
            n=len(V)
            for idxs in combinations(range(n),5):
                sub=induced(adj,list(idxs))
                sk=short_key(sub)
                if sk in FORBIDDEN_SHORT:
                    hits[sk].append((tuple(V[k] for k in idxs),ref,"ph5-induced"))
                    sub_found+=1
    print(f"  induced 5-subtournaments matching forbidden: {sub_found}")
    for f in FORBIDDEN_SHORT:
        print(f"   forbidden {f}: {'HIT '+str(hits[f][:3]) if hits[f] else 'not produced'}")

    print("\n=== SUMMARY ===")
    total_seen = len(seen_short | res_short)
    print(f"distinct n=5 short-keys produced as TOP-LEVEL tournaments: {len(seen_short)}")
    print(f"distinct n=5 short-keys via residue-exhaust: {len(res_short)}")
    any_hit=False
    for f in FORBIDDEN_SHORT:
        if hits[f]:
            any_hit=True
            print(f"REFUTED target {f}: {len(hits[f])} witnesses, e.g. {hits[f][:5]}")
        else:
            print(f"NOT realized: {f}")
    if not any_hit:
        print("NO forbidden class realized in any phase -> claim holds across this search.")
    return hits, seen_short, res_short

if __name__=="__main__":
    main()
