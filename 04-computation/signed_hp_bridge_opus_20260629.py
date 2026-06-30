"""
Bridge the palindromic-polynomial world (THM-084/088, reversal involution) to the
dihedral decomposition (a,b,f). Compute for Paley T_p:
  H        = #Hamiltonian paths (all-forward orderings) = [x^{n-1}] F
  S        = signed HP count = sum_{HP P} sgn(P)        = [x^{n-1}] SF
  f, b     = dihedral reflection-fixed count / sign multiplicity (from prior work)
and look for a clean relation linking S to f or b -- the "reversal sign" unification.

sgn of an ordering: append at the end; new inversions = #{already-placed v > nxt}.
"""
from math import comb

def QR(p): return set((x*x)%p for x in range(1,p))

def H_and_S(p):
    """H = #all-forward orderings (HPs); S = signed sum over them. DP with sign."""
    qr=QR(p)
    adj=[[ (j-i)%p in qr for j in range(p)] for i in range(p)]
    # dp[mask][last] = (count, signed_count) over all-forward prefixes covering mask ending last
    Hc=[[0]*p for _ in range(1<<p)]
    Sc=[[0]*p for _ in range(1<<p)]
    for v in range(p):
        Hc[1<<v][v]=1; Sc[1<<v][v]=1   # single vertex: sign +1
    higher=[ (~((1<<(v+1))-1)) for v in range(p)]  # bits > v
    for mask in range(1<<p):
        hr=Hc[mask]; sr=Sc[mask]
        for last in range(p):
            hc=hr[last]
            if not hc and not sr[last]: continue
            if not (mask>>last)&1: continue
            al=adj[last]; sc=sr[last]
            for nxt in range(p):
                if (mask>>nxt)&1 or not al[nxt]: continue
                inv=bin(mask & higher[nxt]).count("1")   # placed vertices > nxt
                sign = -1 if inv&1 else 1
                nm=mask|(1<<nxt)
                Hc[nm][nxt]+=hc
                Sc[nm][nxt]+=sc*sign
    full=(1<<p)-1
    return sum(Hc[full]), sum(Sc[full])

# dihedral f (anti-palindromic for p=3mod4; 0 for p=1mod4) and b
def f_count(p):
    qr=QR(p); m=(p-1)//2
    if p%4==1: return 0
    pairs=[(a,(p-a)%p) for a in range(1,m+1)]; seq=[]; used=[False]*m; cnt=0
    def arc(i,j): return (j-i)%p in qr
    def dfs(pos):
        nonlocal cnt
        if pos==m:
            full=seq+[0]+[(-x)%p for x in reversed(seq)]
            if all(arc(full[i],full[i+1]) for i in range(p-1)): cnt+=1
            return
        for pi in range(m):
            if used[pi]: continue
            for val in pairs[pi]:
                if pos>0 and not arc(seq[-1],val): continue
                used[pi]=True; seq.append(val); dfs(pos+1); seq.pop(); used[pi]=False
    dfs(0); return cnt

print(f"{'p':>3} {'p%4':>4} {'H':>14} {'S=signed HP':>14} {'f':>8} {'b':>12} {'S vs ?':>20}")
for p in [3,5,7,11,13]:
    H,S=H_and_S(p); f=f_count(p)
    b=(H-p*f)//(2*p)
    # candidate relations
    rels=[]
    if f and S%f==0: rels.append(f"S={S//f}*f")
    if b and S%b==0: rels.append(f"S={S//b}*b")
    if S==0: rels.append("S=0")
    if H and S%H==0: rels.append(f"S={S//H}*H")
    print(f"{p:>3} {p%4:>4} {H:>14} {S:>14} {f:>8} {b:>12} {'; '.join(rels):>20}")
print("\nSF sign (THM-088): C(p,2) parity -> p=1mod4 palindromic(even), p=3mod4 anti-palindromic(odd).")
for p in [3,5,7,11,13]:
    print(f"  p={p}: C(p,2)={comb(p,2)} ({'even->palindromic' if comb(p,2)%2==0 else 'odd->ANTI-palindromic'})")
