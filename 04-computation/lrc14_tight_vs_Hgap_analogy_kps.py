"""kps-2026-07-04: analogy {12,24} (tight LRC coverers=AP,GW) <-> {7,21} (forbidden H, PROVED
THM-029/079). Test the hard bridge: H(tournament) of the tight families; + structural analysis."""
from itertools import combinations
def H_count(n, adj):  # Hamiltonian PATH count via subset DP. adj[i][j]=True if i->j
    from functools import lru_cache
    full=(1<<n)-1
    # dp[mask][v] = # ham paths on 'mask' ending at v
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            if not (mask>>v)&1 or dp[mask][v]==0: continue
            for u in range(n):
                if (mask>>u)&1 or not adj[v][u]: continue
                dp[mask|(1<<u)][u]+=dp[mask][v]
    return sum(dp[full][v] for v in range(n))
def rot_tournament(residues):
    """runners at 'residues' mod 14; i->j if (r_j - r_i) mod 14 in {1..6}; tie(=7) by index order."""
    n=len(residues); adj=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=(residues[j]-residues[i])%14
            if 1<=d<=6: adj[i][j]=True
            elif d==7: adj[i][j]= (i<j)   # break diameter tie by order
            # d in 8..13 => j->i (handled when i,j swap)
    return adj
def is_tournament(adj,n):
    return all((adj[i][j] ^ adj[j][i]) for i in range(n) for j in range(n) if i!=j)

# AP {1..13}: residues 1..13; GW {1..11,13,24}: residues {1..11,13, 24%14=10}
AP=[i%14 for i in range(1,14)]
GW=[i%14 for i in [1,2,3,4,5,6,7,8,9,10,11,13,24]]
print("AP residues mod14:", AP)
print("GW residues mod14:", GW, "(collision at 10: runners 10 & 24; missing 12)")
for name,res in [("AP",AP),("GW",GW)]:
    adj=rot_tournament(res); n=len(res)
    tour=is_tournament(adj,n)
    if tour:
        h=H_count(n,adj)
        print(f"  {name}: H(rotational tournament, n={n}) = {h}  ; H mod 7 = {h%7}, mod 21 = {h%21}, /7={h/7 if h%7==0 else '-'}")
    else:
        print(f"  {name}: NOT a valid tournament (collision => not simple); H undefined at the collision config")
print()
# STRUCTURAL analogy
print("STRUCTURAL ANALOGY:")
print(f"  {{12,24}} = 12*{{1,2}}   base=12=n-2 (n=14), mult set {{1,2}}, ratio 2 (GW=double)")
print(f"  {{7,21}}  =  7*{{1,3}}   base=7=n/2  (n=14), mult set {{1,3}}, ratio 3")
print(f"  MEETING POINT: 84 = 12*7 = 4*21  = the FIRST COVERING residue-liar {{1..11,13,84}}")
print(f"    84 = 6*14 (blocks 14-grid); lonely at 37/89, 89=F11; 84=12*7=(tight base)*(heptagon)=4*21")
print(f"  H-formula H=1+2a1+4a2: H=7<=>a1=3 (forbidden); H=21<=>a1+2a2=10 (six-way block)")
print(f"  M-formula M=k/(12k+5): tight (1/14) at k=1(X=12,AP),2(X=24,GW); loose k>=3")
print(f"  => tight k in {{1,2}} <-> forbidden cycle-count {{3,10}} (a1 or a1+2a2)?  3,10 vs 1,2: ?")
print(f"  ratios: 24/12=2, 21/7=3; 2*3=6=ord(14 mod Phi6=183) [Eisenstein]; golden ord=16 [F49]")
# is H(loose residue-liars) related? compute for a couple
print()
for X in [36,84]:
    res=[i%14 for i in list(range(1,12))+[13,X]]
    adj=rot_tournament(res); n=13
    if is_tournament(adj,n):
        h=H_count(n,adj); print(f"  {{1..11,13,{X}}} residues {sorted(set(res))}: H={h} (mod7={h%7})")
    else:
        print(f"  {{1..11,13,{X}}} residues {res}: collision (X%14={X%14}) => degenerate")
