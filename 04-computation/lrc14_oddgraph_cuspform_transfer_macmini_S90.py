#!/usr/bin/env python3
"""mac-mini-S90: BUILD ATTEMPT -- the tournament odd-graph -> cusp-form f_14 transfer.
Target (HYP-3768 honest negative): the ι-EVEN Dedekind sum s(n,Phi6) is BLIND to tightness (AP & GW
share it). The tight/non-tight distinction = the last bit = the ι-ODD cusp form f_14 (curve 14a).
TEST whether the tournament ι-odd object (the covering-min phase-cloud rotational tournament, and the
tight-config units tournament) (1) SEPARATES tight from covering where Dedekind fails, and (2) matches
f_14's arithmetic (rank 0, torsion Z/6, conductor 14, Klein-four cusps, a_p)."""
from math import gcd
from itertools import combinations
Phi6=lambda n: n*n-n+1

def rot_tournament(points, mod):
    """rotational tournament on residues `points` mod `mod` (odd): i->j iff (p_j-p_i) mod mod in first half."""
    P=sorted(set(x%mod for x in points)); m=len(P); half=mod//2
    adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            d=(P[j]-P[i])%mod
            adj[i][j]=1 if 1<=d<=half else 0
    return P,adj
def scores(adj): return [sum(r) for r in adj]
def c3(adj):
    m=len(adj); C=sum(1 for a,b,c in combinations(range(m),3)
        if (adj[a][b]+adj[b][c]+adj[c][a]) in (0,3))
    # cyclic triangles = C(m,3) - transitive; count cyclic directly:
    cyc=0
    for a,b,c in combinations(range(m),3):
        e=adj[a][b]+adj[b][c]+adj[c][a]+adj[b][a]+adj[c][b]+adj[a][c]  # =3 always
        s=adj[a][b]+adj[b][c]+adj[c][a]
        if s==1 or s==2: cyc+=1
    return cyc
def ham_paths(adj):
    m=len(adj)
    if m>12: return None
    from functools import lru_cache
    full=(1<<m)-1
    import sys; sys.setrecursionlimit(10000)
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def f(mask,last):
        if mask==full: return 1
        t=0
        for nx in range(m):
            if not (mask>>nx)&1 and adj[last][nx]:
                t+=f(mask|(1<<nx),nx)
        return t
    return sum(f(1<<s,s) for s in range(m))

n=14; D=Phi6(n)   # 183
# (1) covering-min phase cloud: {n*k mod D : k=1..n-2} + killer(-n) + observer(0)
cloud=[(n*k)%D for k in range(1,n-1)]+[(-n)%D]+[0]
P,adj=rot_tournament(cloud,D)
sc=scores(adj); cyc=c3(adj); H=ham_paths(adj)
print(f"COVERING-MIN phase-cloud tournament (n=14, mod {D}, {len(P)} vertices):")
print(f"  scores={sorted(sc)}  (HYP-3813 claim: mostly 6, one 7 one 5)")
print(f"  3-cycles c3={cyc}  (HYP-3813 claim: 90);  H (Ham paths)={H}  (odd by Redei: {H%2==1 if H else '?'})")
print(f"  score variance={sum((s-sum(sc)/len(sc))**2 for s in sc)/len(sc):.3f} (0=regular=tight extremal)")

# (2) tight AP config: phases {i/14} = residues {1..13} mod 14; rotational on the UNITS (Z/14)* and on all
units=[a for a in range(1,14) if gcd(a,14)==1]  # (Z/14)* order 6
print(f"\nTIGHT AP config: (Z/14)* = {units} (order phi(14)=6); antipodal pairs {[(a,14-a) for a in units if a<7]}")
print(f"  (Z/14)* is cyclic Z/6 (gen 3: {[pow(3,k,14) for k in range(6)]}); curve 14a torsion = Z/6  -> MATCH?")

# (3) f_14 / curve 14a arithmetic
def ap14(p):
    if p in(2,7): return None
    cnt=1
    for x in range(p):
        for y in range(p):
            if (y*y+x*y+y-(x*x*x+4*x-6))%p==0: cnt+=1
    return p+1-cnt
print("\nf_14 (14a) data: conductor 14=2*7, rank 0, torsion Z/6Z, X_0(14) cusps = Klein four V_4")
print("  a_p:", {p:ap14(p) for p in [3,5,11,13,17,19]})

# (4) SEPARATION TEST: does a tournament ι-odd invariant separate tight (AP,GW non-covering) from covering?
print("\nSEPARATION TEST (the ι-odd detector the Dedekind sum LACKS):")
for name,S,ts_num,ts_den in [("AP {1..13} TIGHT",list(range(1,14)),1,14),
                             ("GW {1..11,13,24} TIGHT",[*range(1,12),13,24],1,14),
                             ("deep well COVERING",list(range(1,13))+[182],14,183),
                             ("near-AP resid COVERING",[*range(1,12),13,84],37,89)]:
    # difference-tournament at tight time t*=ts_num/ts_den: i->j iff frac((v_i-v_j)t*) in (0,1/2)
    m=len(S); adj2=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            fr=((S[i]-S[j])*ts_num % ts_den)/ts_den
            adj2[i][j]=1 if 0<fr<0.5 else 0
    sc2=scores(adj2); var=sum((s-(m-1)/2)**2 for s in sc2)/m
    print(f"  {name:26s}: score var={var:6.3f}, c3={c3(adj2)}, regular={var<1e-9}")
print("\n=> if tight families give REGULAR (var~0) and covering give IRREGULAR, the tournament SEPARATES")
print("   what the Dedekind sum cannot -- the ι-odd tightness detector = the cusp-form's combinatorial shadow.")
