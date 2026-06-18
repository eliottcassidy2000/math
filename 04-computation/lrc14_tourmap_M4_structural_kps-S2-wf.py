# Deep structural analysis of METHOD 4 (section rotational tournament) and M6.
# M4: vertex=runner; at grid time a in (Z/14)*, section r_i=(v_i*a) mod 14;
#     arc i->j iff (r_i - r_j) mod 14 in {1..6}.
# Claim to test: the realized iso classes are EXACTLY the "rotational/circulant"
# tournaments that embed into the 14-cycle circular tournament, and the forbidden
# classes are non-circular-embeddable ones.
#
# We:
#  (1) PROVE forbiddenness is structural by characterizing M4 abstractly: M4 is the
#      restriction of the C_14 ROTATIONAL TOURNAMENT (Z/14 with i->j iff (i-j)%14 in 1..6)
#      to the 13 USED residues {1,...,13} (residue 0 = lonely-forbidden since observer).
#      Actually any tournament realized by M4 is an induced subtournament of R_14
#      (rotational tournament on Z/14, restricted to nonzero distinct residues) up to the
#      tie-break. So the realized set = {induced subtournaments of R_13'} basically.
#  (2) Enumerate the full set of induced subtournaments of the rotational tournament on
#      Z/14 \ {0} = {1..13} on any m-subset of DISTINCT residues, for m=3,4,5,6, and
#      compare to M4's realized set. If they match, the M4 realized set is fully
#      characterized and the forbidden set is PROVED (modulo loneliness allowing all
#      distinct-residue subsets, which it does via perfect SDR).
#  (3) Test much larger speed windows to confirm no new classes appear.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import importlib.util, os

spec = importlib.util.spec_from_file_location(
    "main", os.path.join(os.path.dirname(__file__), "lrc14_tourmap_danger-interval_kps-S2-wf.py"))
mn = importlib.util.module_from_spec(spec); spec.loader.exec_module(mn)

def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

# ---------- The rotational tournament R_14 restricted to residue subsets ----------
def rot_adj_from_residues(res):
    # res: list of DISTINCT residues mod 14 (the "sections")
    m=len(res)
    adj=[[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            d=(res[i]-res[j])%14
            if 1<=d<=6: adj[i][j]=True
            elif 7<=d<=13: adj[i][j]=False
            else:
                # d==0 impossible (distinct); but keep tie-break
                adj[i][j]=(res[i]<res[j])
    return adj

def induced_subtournaments_of_R14(m, allow_zero=False, allow_collision=False):
    """All iso classes arising from choosing m residues from Z/14.
       allow_zero: include residue 0 (the lonely-forbidden section).
       Distinct residues only (SDR / perfect)."""
    pool = list(range(0,14)) if allow_zero else list(range(1,14))
    seen={}
    for res in combinations(pool, m):
        adj=rot_adj_from_residues(list(res))
        cn=mn.canon(adj,m)
        if cn not in seen:
            seen[cn]=(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),tuple(res))
    return seen

# ---------- M4 realized set (recompute, large window) ----------
def m4_realized(m, vmax):
    sets=[c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1]
    real={}
    for S in sets:
        S=list(S)
        for a in range(1,14):
            if gcd(a,14)!=1: continue
            adj=mn.method4_adj(S,a); cn=mn.canon(adj,m)
            if cn not in real:
                real[cn]=(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),a))
    return real

print("="*72)
print("M4 structural characterization: realized = induced subtournaments of R_14?")
print("="*72)
for m in [3,4,5,6]:
    free=mn.FREE[m]
    # induced subtournaments of rotational R_14 on nonzero residues (distinct)
    rot_nonzero = induced_subtournaments_of_R14(m, allow_zero=False)
    rot_withzero = induced_subtournaments_of_R14(m, allow_zero=True)
    vmax = {3:20,4:16,5:13,6:11}[m]
    real = m4_realized(m, vmax)
    rkeys=set(real); rotk=set(rot_nonzero); rotz=set(rot_withzero)
    print(f"\n m={m} (free A000568={free}):")
    print(f"   M4 realized (vmax={vmax}):            {len(rkeys)}")
    print(f"   induced sub R_14 nonzero residues:   {len(rotk)}")
    print(f"   induced sub R_14 incl residue 0:     {len(rotz)}")
    print(f"   M4 == R_14(nonzero) ? {rkeys==rotk}    M4 subset of R_14(all)? {rkeys<=rotz}")
    extra = rkeys - rotz
    missing_from_m4 = rotk - rkeys
    if extra: print(f"   M4 has classes NOT in any rotational embedding: {len(extra)} (unexpected)")
    if missing_from_m4: print(f"   rotational classes M4 never hits: {len(missing_from_m4)}")
    # forbidden = free - realized
    print(f"   FORBIDDEN by M4: {free-len(rkeys)} classes")

# ---------- Characterize WHICH classes are rotational-embeddable in general ----------
print()
print("="*72)
print("Rotational embeddability count vs free set (this is the THEORETICAL ceiling)")
print("="*72)
for m in [3,4,5,6]:
    rotk=induced_subtournaments_of_R14(m, allow_zero=False)
    print(f" m={m}: {len(rotk)}/{mn.FREE[m]} iso classes embed as induced sub of R_14 (nonzero residues)")
