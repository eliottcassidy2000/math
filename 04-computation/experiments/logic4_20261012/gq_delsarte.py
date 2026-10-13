"""Delsarte LP bound (= Schrijver's theta') for alpha(G_q), G_q = Cay(F_{q^2}, mu_{q+1}), q prime:
maximize sum f over f: F_{q^2} -> R_{>=0}, f(0) = 1, f = 0 on mu_{q+1}, Fourier transform f^ >= 0.
By symmetry f is constant on the norm circles C_n = {N(z) = n}, n in F_q^x (Frobenius fixes F_q for prime q).
Characters psi_a(z) = exp(2 pi i Tr(az)/q); f^(a) = sum_n f_n K(a,n), K(a,n) = sum_{z in C_n} psi_a(z) depends on N(a).
Also: explicit independent sets that are unions of parallel lines in a non-square-norm direction."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
def run(q):
    adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
    elems2,idx2,mul,add2,pw=field(p,2*e)
    # norm N(z) = z^(q+1) in F_q = prime field, trace Tr(z) = z + z^q
    def N(z): return pw(z,q+1)[0]   # element of F_q = coefficient tuple (c,0)
    def Tr(z): return add2(z,pw(z,q))[0]
    norms=[N(z) for z in elems]
    circles={m:[i for i,z in enumerate(elems) if norms[i]==m] for m in range(1,q)}
    assert all(len(c)==q+1 for c in circles.values())
    # representative a per norm class (nonzero), plus a = 0
    reps={}
    for i,z in enumerate(elems):
        if norms[i] not in reps and i!=0: reps[norms[i]]=z
    # K(a,n) real parts
    def K(a,m):
        return sum(np.cos(2*np.pi*Tr(mul(a,elems[i]))/q) for i in circles[m])
    ns=list(range(1,q))
    Kmat={na: {m: K(reps[na],m) for m in ns} for na in ns}
    # LP: variables f_n for n in ns; f_1 = 0 (unit circle); maximize 1 + (q+1) sum f_n
    # constraints: for each character class na: 1 + sum_m f_m K(na,m) >= 0  ; and a = 0: 1 + (q+1) sum f_m >= 0 trivial
    c=-np.array([(q+1) for m in ns],dtype=float)
    A=[];b=[]
    for na in ns:
        A.append([-Kmat[na][m] for m in ns]); b.append(1.0)
    bounds=[(0,0) if m==1 else (0,None) for m in ns]
    res=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method='highs')
    lp=1-res.fun
    print("q=%d: Delsarte LP bound alpha <= %.4f  (n/5 = %.2f, n/6 = %.2f); optimal f by circle: %s"%(q,lp,n/5,n/6,{m:round(v,4) for m,v in zip(ns,res.x) if v>1e-9}))
    # explicit line-union independent sets: direction v with N(v) non-square; cosets a + F_q v; compatibility
    sq={(t*t)%q for t in range(1,q)}
    v=next(z for i,z in enumerate(elems) if i!=0 and norms[i] not in sq)
    line=[idx[mul(z,v)] for z in elems if pw(z,q)==z]  # F_q multiples: elements fixed by Frobenius
    line=sorted(set(idx[mul(elems[i],v)] for i in range(n) if pw(elems[i],q)==elems[i]))
    assert len(line)==q
    cosets=[]; seen=set()
    for i in range(n):
        if i in seen: continue
        cs=sorted(idx[add(elems[i],elems[j])] for j in line); cosets.append(cs); seen.update(cs)
    assert len(cosets)==q
    def compatible(c1,c2): return all(j not in adj[i] for i in c1 for j in c2)
    import networkx as nx
    H=nx.Graph(); H.add_nodes_from(range(q))
    for i in range(q):
        for j in range(i+1,q):
            if compatible(cosets[i],cosets[j]): H.add_edge(i,j)
    cl=max(nx.find_cliques(H),key=len)
    S=[v_ for i in cl for v_ in cosets[i]]
    assert all(j not in adj[i] for i in S for j in S)
    print("   line-union independent set: %d compatible parallel lines -> alpha >= %d"%(len(cl),len(S)))
for q in [int(a) for a in sys.argv[1:]] or [5,7,11,13,17,19,23]:
    run(q)
