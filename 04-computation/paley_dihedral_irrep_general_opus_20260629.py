"""
THM-131 GENERALIZATION (finished): D_{2p} irrep decomposition of the Hamiltonian-path
permutation representation of the Paley tournament T_p, p = 3 (mod 4) prime.
opus 2026-06-29 (dihedral<->tournament thread).

PROOF (rigorous). G=D_{2p}=<r,s | r^p=s^2=1, srs=r^-1>, p odd prime. Conjugacy classes:
{e}; {r^j,r^-j}, j=1..(p-1)/2; and ONE class of all p reflections (p odd). Irreps:
trivial rho0, sign rho1 (1 on rotations, -1 on reflections), and (p-1)/2 two-dim rho_k
(chi(r^j)=2cos(2pi jk/p), chi(reflection)=0).

V = permutation rep on the H Hamiltonian paths (directed). Action: r:(v_i)->(v_i+1);
s:(v_0..v_{p-1})->(-v_{p-1}..-v_0) (reverse+negate; uses T_p ~ T_p^op, THM-127).
  chi_V(e)=H.
  chi_V(r^j)=0 for j!=0: a directed HP fixed by +j forces source v_0+j=v_0 => j=0.
  chi_V(reflection)=f (same for all reflections; one conjugacy class), f=#anti-palindromic
     HPs (v_i = -v_{p-1-i}, middle vertex 0).
Multiplicities m(rho)=(1/2p)[H*chi_rho(e) + f * sum_{p reflections} chi_rho(refl)]:
  a := m(triv) = (H + p f)/(2p)
  b := m(sign) = (H - p f)/(2p)        # the chirality / sign-isotypic multiplicity
  c := m(2-dim)= H/p   (each of the (p-1)/2)
Identities: a+b = H/p (= #Z_p-orbits); a-b = f; a = #D_{2p}-orbits of HPs (Burnside);
dim: a+b+2c*(p-1)/2 = H. Chirality fraction b/(a+b) = (1 - p f / H)/2 -> 1/2 (f=o(H/p)).
"""
from math import factorial

def QR(p): return set((x*x)%p for x in range(1,p))

def H_count(p):
    qr=QR(p); adj=[[ (j-i)%p in qr for j in range(p)] for i in range(p)]
    dp=[[0]*p for _ in range(1<<p)]
    for v in range(p): dp[1<<v][v]=1
    for mask in range(1<<p):
        row=dp[mask]
        for last in range(p):
            c=row[last]
            if not c or not (mask>>last)&1: continue
            al=adj[last]
            for nxt in range(p):
                if (mask>>nxt)&1 or not al[nxt]: continue
                dp[mask|(1<<nxt)][nxt]+=c
    full=(1<<p)-1
    return sum(dp[full][v] for v in range(p))

def f_count(p):
    qr=QR(p); m=(p-1)//2
    pairs=[(a,(p-a)%p) for a in range(1,m+1)]
    seq=[]; used=[False]*m; count=0
    def arc(i,j): return (j-i)%p in qr
    def dfs(pos):
        nonlocal count
        if pos==m:
            full=seq+[0]+[(-x)%p for x in reversed(seq)]
            if all(arc(full[i],full[i+1]) for i in range(p-1)): count+=1
            return
        for pi in range(m):
            if used[pi]: continue
            for val in pairs[pi]:
                if pos>0 and not arc(seq[-1],val): continue
                used[pi]=True; seq.append(val); dfs(pos+1); seq.pop(); used[pi]=False
    dfs(0); return count

print(f"{'p':>3} {'H(T_p)':>16} {'f':>8} {'a':>14} {'b=sign':>14} {'c':>14} {'b/(a+b)':>8}")
for p in [3,7,11,19]:
    H=H_count(p); f=f_count(p)
    a=(H+p*f)//(2*p); b=(H-p*f)//(2*p); c=H//p
    assert (H+p*f)%(2*p)==0 and (H-p*f)%(2*p)==0 and H%p==0
    assert a+b+2*c*((p-1)//2)==H and a-b==f
    print(f"{p:>3} {H:>16} {f:>8} {a:>14} {b:>14} {c:>14} {(1-p*f/H)/2:>8.4f}")
# extra f data point (H(23) too large for this DP, but f is a pruned DFS):
f23=f_count(23)
print(f"\nf(23)={f23} (extends the reflection-fixed/self-reverse-complement HP sequence 1,9,185,573057,...)")
print("chirality fraction b/(a+b)=(1-pf/H)/2 -> 1/2 since f=o(H/p): Paley HP-rep is asymptotically maximally chiral.")
