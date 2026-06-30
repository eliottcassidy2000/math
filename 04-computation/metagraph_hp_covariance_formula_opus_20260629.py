"""
PROVE on the finite metagraph: H = sum_pi 1[pi is a Ham path], a sum of dependent indicators.
Claim (pair-covariance):  Cov(1[pi HP],1[pi' HP]) = 2^{-2(n-1)}(2^c - 1)  if COMPATIBLE
                                                   = -2^{-2(n-1)}          if CONFLICT,
where c = # common directed consecutive-arcs, conflict = some pair consecutive-opposite in both.
Then Var(H) = 2^{-2(n-1)} [ sum_{compatible(pi,pi')} 2^{c} - (n!)^2 ].
Verify against direct Var(H)=E[H^2]-E[H]^2. Also: does FKG hold (all Cov>=0)? NO -> mixed sign.
"""
from itertools import permutations
from fractions import Fraction as F
def arcs(p):  # directed consecutive arcs as frozenset of ordered pairs
    return set((p[i],p[i+1]) for i in range(len(p)-1))
def Hcount(n,adj):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            av=adj[last]&~mask
            while av:
                nb=av&-av; nx=nb.bit_length()-1; dp[mask|nb][nx]+=c; av^=nb
    return sum(dp[(1<<n)-1])
def direct_var(n):
    E=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(E); s1=0;s2=0
    for bits in range(1<<m):
        adj=[0]*n
        for k,(i,j) in enumerate(E):
            if (bits>>k)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        h=Hcount(n,adj); s1+=h; s2+=h*h
    N=1<<m; return F(s2,N)-F(s1,N)**2
for n in range(3,7):
    P=list(permutations(range(n)))
    base=F(1,4**(n-1))
    var_formula=F(0); npos=0;nneg=0; sum2c=0; ncompat=0
    for p in P:
        ap=arcs(p)
        for q in P:
            aq=arcs(q)
            # conflict: some ordered pair (a,b) in ap with (b,a) in aq
            conflict=any((b,a) in aq for (a,b) in ap)
            if conflict:
                cov=-base; nneg+=1
            else:
                c=len(ap & aq)
                cov=base*(2**c-1); 
                sum2c+=2**c; ncompat+=1
                if cov>0: npos+=1
            var_formula+=cov
    vf=var_formula; vd=direct_var(n)
    # closed reduction
    red=base*(sum2c - len(P)**2)
    print(f"n={n}: Var(H) direct={vd}  formula={vf}  reduction=2^-2(n-1)[Sum_compat 2^c - (n!)^2]={red}  match={vd==vf==red}")
    print(f"      #ordered pairs={len(P)**2}  compatible={ncompat} (cov>=0)  conflict={nneg} (cov<0)  => FKG FAILS (mixed sign): {nneg>0}")
