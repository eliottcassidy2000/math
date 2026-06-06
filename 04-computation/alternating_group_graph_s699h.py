"""The ALTERNATING GROUP GRAPH AG_n (Cayley(A_n, {(1 2 i),(1 i 2)})) and the full 3-cycle Cayley
graph, as members of the forbidden-distance Cayley family (HN/LRC/UD unification, S699g). Compute
χ, independence α, spectrum (= REPRESENTATION THEORY = the non-abelian Fourier transform; the
character ratios are the Bessel/Dirichlet analog), Hoffman bounds, triangles, as n varies.
opus-2026-06-06-S699h."""
from itertools import permutations
try:
    import numpy as np; HAVE=True
except Exception: HAVE=False
def parity(p):
    p=list(p); seen=[False]*len(p); par=0
    for i in range(len(p)):
        if seen[i]: continue
        j=i; c=0
        while not seen[j]: seen[j]=True; j=p[j]; c+=1
        par+=c-1
    return par%2
def An(n): return [p for p in permutations(range(n)) if parity(p)==0]
def mul(p,q): return tuple(p[q[i]] for i in range(len(q)))   # (p∘q)
def inv(p):
    r=[0]*len(p)
    for i,x in enumerate(p): r[x]=i
    return tuple(r)
def threecycle(a,b,c,n):  # cycle a->b->c->a as a permutation
    p=list(range(n)); p[a]=b; p[b]=c; p[c]=a; return tuple(p)
def cayley(n, gens):
    V=An(n); idx={g:i for i,g in enumerate(V)}; N=len(V)
    adj=[[0]*N for _ in range(N)]
    for g in V:
        for s in gens:
            h=mul(s,g)
            if h in idx and idx[h]!=idx[g]: adj[idx[g]][idx[h]]=1; adj[idx[h]][idx[g]]=1
    return V,adj
def greedy_chi(adj):
    N=len(adj); color=[-1]*N
    order=sorted(range(N), key=lambda v:-sum(adj[v]))
    for v in order:
        used={color[w] for w in range(N) if adj[v][w] and color[w]>=0}
        c=0
        while c in used: c+=1
        color[v]=c
    return max(color)+1
def greedy_indep(adj):
    N=len(adj); order=sorted(range(N), key=lambda v:sum(adj[v])); S=[]; banned=[False]*N
    for v in order:
        if not banned[v]:
            S.append(v); banned[v]=True
            for w in range(N):
                if adj[v][w]: banned[w]=True
    return len(S)
def triangles(adj):
    N=len(adj); t=0
    for i in range(N):
        for j in range(i+1,N):
            if not adj[i][j]: continue
            for k in range(j+1,N):
                if adj[i][k] and adj[j][k]: t+=1
    return t
def main():
    print("AG_n = Cayley(A_n, {(1 2 i),(1 i 2): 3<=i<=n}):  n | |A_n| | deg | triangles | χ(greedy) | α | spectrum bounds")
    for n in range(3,7):
        gens=[threecycle(0,1,i,n) for i in range(2,n)]+[inv(threecycle(0,1,i,n)) for i in range(2,n)]
        V,adj=cayley(n,gens); N=len(V); deg=sum(adj[0])
        if N<=360:
            tri=triangles(adj) if N<=120 else None
            chi=greedy_chi(adj); al=greedy_indep(adj)
            if HAVE:
                A=np.array(adj,float); ev=np.linalg.eigvalsh(A); lmin=ev[0]; lmax=ev[-1]
                hoff_chi=1-lmax/lmin; hoff_alpha=N*(-lmin)/(deg-lmin)
                print(f"  {n} | {N:3d} | {deg} | tri={tri} | χ≈{chi} | α≈{al} | λmin={lmin:.3f} λmax={lmax:.3f}; Hoffman χ≥{hoff_chi:.2f}, α≤{hoff_alpha:.1f}")
            else:
                print(f"  {n} | {N:3d} | {deg} | tri={tri} | χ≈{chi} | α≈{al} | (no numpy)")
    print()
    print("FULL 3-CYCLE Cayley graph of A_n (g~h iff gh^{-1} a 3-cycle): spectrum = REPRESENTATION THEORY")
    print("  n | |A_n| | deg=2C(n,3) | λmin | Hoffman χ≥ | (eigenvalues = |C|·χ_ρ(3cyc)/dim ρ)")
    for n in range(3,7):
        allg=[threecycle(a,b,c,n) for a in range(n) for b in range(n) for c in range(n) if len({a,b,c})==3]
        # dedupe (each 3-cycle once) — threecycle(a,b,c) for ordered distinct triples gives each cycle once per rotation*... keep set
        gens=list(set(allg))
        V,adj=cayley(n,gens); N=len(V); deg=sum(adj[0])
        if HAVE and N<=360:
            A=np.array(adj,float); ev=np.linalg.eigvalsh(A); lmin=ev[0]
            print(f"  {n} | {N:3d} | {deg} | λmin={lmin:.3f} | χ≥{1-deg/lmin:.2f} | distinct eigs={sorted(set(round(x,2) for x in ev))[:8]}")
if __name__=='__main__': main()
