"""
THREAD A closed form. W(n)=Sum_{compatible pi,pi'} 2^c = n! * G(n), where (left-regular symmetry)
G(n)=Sum_{sigma in S_n} g(sigma), g(sigma)=2^{#successions(sigma_{k+1}=sigma_k+1)} if NO value-descent-
adjacency (sigma_{k+1}=sigma_k-1), else 0.
CLAIM (inclusion-exclusion on bonds): G(n)=Sum_{T subset of [n-1], ALL maximal runs EVEN length} (n-|T|)! * 2^{#runs(T)}.
Then Var(H)=(n!/4^{n-1})(G(n)-n!).
Verify all three: direct g-sum, even-run formula, and the original Var(H).
"""
from itertools import permutations
from fractions import Fraction as F
def g(sigma):
    n=len(sigma); succ=0
    for k in range(n-1):
        if sigma[k+1]==sigma[k]+1: succ+=1
        if sigma[k+1]==sigma[k]-1: return 0
    return 2**succ
def G_direct(n):
    return sum(g(s) for s in permutations(range(n)))
def runs(T,L):
    # T set of bond positions in [0..L-1]; return list of run lengths (maximal consecutive)
    rs=[]; i=0; T=sorted(T)
    if not T: return rs
    cur=1
    for k in range(1,len(T)):
        if T[k]==T[k-1]+1: cur+=1
        else: rs.append(cur); cur=1
    rs.append(cur); return rs
def G_evenrun(n):
    import itertools
    L=n-1; tot=0
    for r in range(L+1):
        for T in itertools.combinations(range(L),r):
            rl=runs(set(T),L)
            if all(x%2==0 for x in rl):
                tot+=__import__('math').factorial(n-len(T))*2**len(rl)
    return tot
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
def Var_direct(n):
    E=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(E); s1=0;s2=0
    for bits in range(1<<m):
        adj=[0]*n
        for k,(i,j) in enumerate(E):
            if (bits>>k)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        h=Hcount(n,adj); s1+=h; s2+=h*h
    N=1<<m; return F(s2,N)-F(s1,N)**2
import math
print(f"{'n':>2} {'G_direct':>10} {'G_evenrun':>10} {'match':>6} {'Var(H) formula':>16} {'Var(H) direct':>16} {'ok':>4}")
seq=[]
for n in range(3,8):
    gd=G_direct(n); ge=G_evenrun(n); seq.append(gd)
    vf=F(math.factorial(n),4**(n-1))*(gd-math.factorial(n))
    vd=Var_direct(n) if n<=6 else None
    print(f"{n:>2} {gd:>10} {ge:>10} {str(gd==ge):>6} {str(vf):>16} {str(vd):>16} {str(vf==vd) if vd is not None else '--':>4}")
print(f"\nG(n) = Sum_sigma g = {seq}  (n=3..7)  [= W(n)/n!];  search OEIS")
# extend G(n) via even-run formula (fast) to n=12
print("extended G(n) via even-run closed form:")
gg=[G_evenrun(n) for n in range(3,13)]
print(f"  n=3..12: {gg}")
