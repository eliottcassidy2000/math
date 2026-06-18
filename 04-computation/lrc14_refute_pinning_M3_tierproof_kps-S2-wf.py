"""
Part 2 of refutation: TIER STRUCTURE PROOF + adversarial escape analysis.

Goal: turn the empirical 'forbidden (2,2,2,2,2)' into an (almost) airtight
residue-level statement, and verify the only escapes are exactly the
degenerate all-equal-residue multisets the claim concedes.
"""
from itertools import combinations, permutations, combinations_with_replacement, product
from math import gcd
from functools import reduce

U14 = [1,3,5,9,11,13]
def depth(r): r%=14; return min(r,14-r)
def prof(v): return tuple(depth((v*a)%14) for a in U14)

# Establish the 3-tier classification of residues by depth profile.
print("="*70)
print("TIER STRUCTURE: classify residues 0..13 by depth profile")
profiles = {}
for r in range(14):
    profiles.setdefault(prof(r), []).append(r)
for p,rs in sorted(profiles.items(), key=lambda kv: -sum(kv[0])):
    print(f"  profile {p} sum={sum(p)}: residues {rs}")

# The PAIRWISE margin between two residues r,s = sum sign(prof(r)[k]-prof(s)[k]).
# Compute the full margin matrix between the 4 'tier representatives':
#   SEVEN = res 7 (profile sum 42)
#   EVEN  = {2,4,6,8,10,12} (profile sum 24) -- but profiles differ within!
#   ODD   = {1,3,5,9,11,13} (profile sum 18)
#   ZERO  = res 0 (profile sum 0, the parked/section-0 runner)
print("="*70)
print("PAIRWISE MARGIN between every ordered residue pair (margin(r over s))")
def margin(r,s):
    pr,ps=prof(r),prof(s)
    gt=sum(1 for k in range(6) if pr[k]>ps[k])
    lt=sum(1 for k in range(6) if pr[k]<ps[k])
    return gt-lt
# show margin matrix
hdr = "      " + "".join(f"{s:4d}" for s in range(14))
print(hdr)
for r in range(14):
    row="".join(f"{margin(r,s):4d}" for s in range(14))
    print(f"  {r:2d}: {row}")

# KEY CLAIM to verify: among DISTINCT residues, every pair has NONZERO margin
# EXCEPT pairs within the SAME profile class with identical profiles.
print("="*70)
print("Which distinct-residue pairs are TIED (margin 0)?")
tied_pairs=[(r,s) for r in range(14) for s in range(r+1,14) if margin(r,s)==0]
print(f"  tied distinct pairs: {tied_pairs}")
# group: these should be exactly same-profile pairs
for r,s in tied_pairs:
    print(f"    ({r},{s}): prof(r)={prof(r)} prof(s)={prof(s)} same={prof(r)==prof(s)}")

# So the M3 tournament on distinct residues is a STRICT order refinement:
# SEVEN beats everything; among EVEN vs ODD, EVEN beats ODD; within tiers,
# margins resolve by sub-profile, and only TRUE profile-twins tie.
# Score (2,2,2,2,2) requires a REGULAR tournament = no dominating vertex.
# But a SEVEN runner (res 7) beats ALL others (margin +6 vs everyone) -> out-deg 4.
# And ZERO (res 0) loses to ALL others -> out-deg 0. So any set containing
# res 7 OR res 0 with 4 distinct others is NOT regular.

print("="*70)
print("DOMINANCE CHECK: res 7 vs all, res 0 vs all")
print(f"  margin(7 over r) for r!=7: {[margin(7,r) for r in range(14) if r!=7]}")
print(f"  margin(r over 0) for r!=0: {[margin(r,0) for r in range(14) if r!=0]}")

# EXHAUSTIVE residue-multiset proof that (2,2,2,2,2) needs all-equal residues.
print("="*70)
print("EXHAUSTIVE: over ALL 5-residue MULTISETS, with FREE tie-break, which")
print("reach (2,2,2,2,2)? And are they ALL all-equal-residue (degenerate)?")
def scores(adj):
    n=len(adj); return tuple(sorted(sum(1 for j in range(n) if adj[i][j]) for i in range(n)))
target=(2,2,2,2,2)
reach=[]
for combo in combinations_with_replacement(range(14),5):
    P=[prof(r) for r in combo]
    fixed=[[None]*5 for _ in range(5)]; tied=[]
    for i in range(5):
        for j in range(i+1,5):
            gt=sum(1 for k in range(6) if P[i][k]>P[j][k]); lt=sum(1 for k in range(6) if P[i][k]<P[j][k])
            m=gt-lt
            if m>0: fixed[i][j]=True; fixed[j][i]=False
            elif m<0: fixed[i][j]=False; fixed[j][i]=True
            else: tied.append((i,j))
    ok=False
    for bits in range(1<<len(tied)):
        adj=[[fixed[i][j] for j in range(5)] for i in range(5)]
        for idx,(i,j) in enumerate(tied):
            if (bits>>idx)&1: adj[i][j]=True; adj[j][i]=False
            else: adj[i][j]=False; adj[j][i]=True
        if scores(adj)==target: ok=True; break
    if ok: reach.append(combo)
# classify reachable multisets
def all_same_profile(combo):
    return len(set(prof(r) for r in combo))==1
all_twin=[c for c in reach if all_same_profile(c)]
not_twin=[c for c in reach if not all_same_profile(c)]
print(f"  multisets reaching (2,2,2,2,2) with free tie-break: {len(reach)}")
print(f"  of these, ALL-same-profile (degenerate twins): {len(all_twin)}")
print(f"  of these, NOT all-same-profile (GENUINE escape?): {len(not_twin)}")
if not_twin:
    print("  !!! GENUINE NON-DEGENERATE ESCAPE (would refute):")
    for c in not_twin[:20]: print(f"     {c}  profiles={[prof(r) for r in c]}")
else:
    print("  -> every (2,2,2,2,2)-reachable multiset is all-same-profile.")
    print("     Same profile => all 10 pairs tied => tournament is decided")
    print("     ENTIRELY by tie-break. Under the smaller-speed tie-break this is")
    print("     the TRANSITIVE order (H=1), NOT regular. Regular only via")
    print("     arbitrary tie orientation, which the claim already concedes.")
    # show the profile classes that admit it
    profs_seen=set(prof(c[0]) for c in all_twin)
    print(f"     profile classes admitting it: {sorted(profs_seen)}")

# Now: do all-same-profile multisets EVER arise from a LONELY / LRC-valid
# DISTINCT-speed primitive set? all-same-profile with DISTINCT residues is
# impossible for the ODD tier? No: odd tier {1,3,5,9,11,13} has TWO sub-profiles.
print("="*70)
print("Can a DISTINCT-residue 5-set be all-same-profile? (needed for regular)")
for p,rs in profiles.items():
    print(f"  profile {p}: {len(rs)} residues {rs} -> can pick 5 distinct? {len(rs)>=5}")
print("  -> only EVEN tier {2,4,6,8,10,12} has >=5 residues with a COMMON profile?")
# check: do all 6 even residues share ONE profile?
print("  even-tier profiles:")
for r in [2,4,6,8,10,12]: print(f"    res {r}: {prof(r)}")
print("  odd-tier profiles:")
for r in [1,3,5,9,11,13]: print(f"    res {r}: {prof(r)}")
