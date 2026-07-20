#!/usr/bin/env python3
"""
death-star-2026-07-20-S59y (HYP-8235) -- compute h(G_8/Z_2) and track growth vs chi.
Build G_8/Z_2 (merged tournament metagraph n=8) via BFS on the iso-class single-
arc-flip graph, WL-refinement canonical form, complement-merge. Then Hadwiger number.
"""
import random, sys
from itertools import permutations, product
from collections import defaultdict, deque
random.seed(8)
N=8; P=[(i,j) for i in range(N) for j in range(i+1,N)]; M=len(P)
def flip(a,arc):
    b=[r[:] for r in a]; i,j=P[arc]; b[i][j],b[j][i]=b[j][i],b[i][j]; return b
def comp(a): return [[a[j][i] if i!=j else 0 for j in range(N)] for i in range(N)]
def refine(a):
    color=[sum(a[i]) for i in range(N)]
    for _ in range(N):
        sig=[(color[i],tuple(sorted(color[j] for j in range(N) if a[i][j])),
              tuple(sorted(color[j] for j in range(N) if a[j][i]))) for i in range(N)]
        order={c:k for k,c in enumerate(sorted(set(sig)))}
        nc=[order[sig[i]] for i in range(N)]
        if nc==color: break
        color=nc
    return color
def canon(a):
    color=refine(a); cells=defaultdict(list)
    for i in range(N): cells[color[i]].append(i)
    groups=[cells[c] for c in sorted(cells)]
    # fast path: all singletons -> unique order
    if all(len(g)==1 for g in groups):
        perm=[g[0] for g in groups]
        return tuple(a[perm[i]][perm[j]] for i in range(N) for j in range(N))
    best=None
    for combo in product(*[list(permutations(g)) for g in groups]):
        perm=[v for pg in combo for v in pg]
        f=tuple(a[perm[i]][perm[j]] for i in range(N) for j in range(N))
        if best is None or f<best: best=f
    return best
def mkey(a): return min(canon(a),canon(comp(a)))
print("building G_8/Z_2 (BFS flip graph, n=8)...", flush=True)
transitive=[[1 if i<j else 0 for j in range(N)] for i in range(N)]
seen={}; rep={}; edges=set(); k0=mkey(transitive); seen[k0]=0; rep[0]=transitive; q=deque([0]); done=0
while q:
    mid=q.popleft(); a=rep[mid]; done+=1
    if done%400==0: print(f"  ...{done} merged classes processed, |seen|={len(seen)}", flush=True)
    for arc in range(M):
        a2=flip(a,arc); k2=mkey(a2)
        if k2 not in seen: nid=len(seen); seen[k2]=nid; rep[nid]=a2; q.append(nid)
        j=seen[k2]
        if j!=mid: edges.add((min(mid,j),max(mid,j)))
V=len(seen); adj=defaultdict(set)
for u,w in edges: adj[u].add(w); adj[w].add(u)
print(f"G_8/Z_2: V={V}, E={len(edges)} (A000568(8)=6880; V=(6880+SC_8)/2 => SC_8={2*V-6880})", flush=True)
degs=sorted(len(adj[v]) for v in range(V))
print(f"  degrees {degs[0]}..{degs[-1]}, mean {2*len(edges)/V:.1f}", flush=True)
def maxclique(nodes,neigh):
    best=[]
    def bk(R,Pp,X):
        nonlocal best
        if not Pp and not X:
            if len(R)>len(best): best=list(R)
            return
        if len(R)+len(Pp)<=len(best): return
        for v in list(Pp):
            bk(R+[v],[w for w in Pp if w in neigh[v]],[w for w in X if w in neigh[v]])
            Pp=[w for w in Pp if w!=v]; X=X+[v]
    bk([],list(nodes),[]); return best
om=maxclique(range(V),{v:adj[v] for v in range(V)})
print(f"  omega = {len(om)}", flush=True)
def kt_minor(t,tries):
    for _ in range(tries):
        neigh={v:set(adj[v]) for v in range(V)}; members={v:{v} for v in range(V)}; active=set(range(V))
        while len(active)>t:
            us=[u for u in active if neigh[u]]
            if not us: break
            u=random.choice(us); w=random.choice(list(neigh[u])); members[u]|=members[w]
            for x in list(neigh[w]):
                neigh[x].discard(w)
                if x!=u: neigh[x].add(u); neigh[u].add(x)
            neigh[u].discard(u); del neigh[w]; active.discard(w)
            if len(active)<=t+2:
                cl=maxclique(active,neigh)
                if len(cl)>=t: return [members[c] for c in cl[:t]]
        cl=maxclique(active,neigh)
        if len(cl)>=t: return [members[c] for c in cl[:t]]
    return None
def verify(branch):
    allv=set()
    for B in branch:
        B=set(B)
        if allv&B: return False
        allv|=B; s=next(iter(B)); st=[s]; vis={s}
        while st:
            x=st.pop()
            for y in adj[x]&B:
                if y not in vis: vis.add(y); st.append(y)
        if vis!=B: return False
    for i in range(len(branch)):
        for j in range(i+1,len(branch)):
            if not any(adj[x]&set(branch[j]) for x in branch[i]): return False
    return True
print("\n=== Hadwiger number of G_8/Z_2 (max t with certified K_t minor) ===", flush=True)
had=0
for t in range(4,26):
    mn=kt_minor(t, 150 if t<=10 else 400)
    if mn and verify(mn):
        had=t; print(f"  K_{t}: YES", flush=True)
    else:
        print(f"  K_{t}: none in budget -> h(G_8/Z_2) >= {had}", flush=True); break
print(f"\n  RESULT: h(G_8/Z_2) >= {had}; chi=7 (=n-1); omega={len(om)}", flush=True)
print(f"  GROWTH: n=7 (chi=6, h>=12); n=8 (chi=7, h>={had}) -> h/chi ~ {had/7:.2f}", flush=True)
