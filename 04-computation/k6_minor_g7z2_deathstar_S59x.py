#!/usr/bin/env python3
"""
death-star-2026-07-20-S59x (HYP-8225) -- build G_7/Z_2 (merged tournament
metagraph, n=7, wiggly single-arc-flip edges, complement-merged) and EXHIBIT a
K_6 minor (Hadwiger t=6 is a theorem, so chi=6 => it exists; we produce explicit
branch sets + verify the certificate), then probe K_7.
Edge def matches chromatic_n7_deep_s314 (single arc flip between iso classes).
"""
import sys, random
from itertools import permutations, product
from collections import defaultdict, deque
random.seed(59)
N=7
P=[(i,j) for i in range(N) for j in range(i+1,N)]   # 21 arcs
M=len(P)

def flip(a, arc):
    b=[row[:] for row in a]; (i,j)=P[arc]
    b[i][j],b[j][i]=b[j][i],b[i][j]
    return b
def comp(a):
    return [[a[j][i] if i!=j else 0 for j in range(N)] for i in range(N)]
def refine(a):
    color=[sum(a[i]) for i in range(N)]
    for _ in range(N):
        sig=[]
        for i in range(N):
            out=tuple(sorted(color[j] for j in range(N) if a[i][j]))
            inn=tuple(sorted(color[j] for j in range(N) if a[j][i]))
            sig.append((color[i],out,inn))
        order=sorted(set(sig)); idx={c:k for k,c in enumerate(order)}
        nc=[idx[sig[i]] for i in range(N)]
        if nc==color: break
        color=nc
    return color
def canon(a):
    color=refine(a)
    cells=defaultdict(list)
    for i in range(N): cells[color[i]].append(i)
    groups=[cells[c] for c in sorted(cells)]
    best=None
    # permute within each cell
    ranges=[list(permutations(g)) for g in groups]
    for combo in product(*ranges):
        perm=[0]*N; pos=0
        # place cells in sorted order of their color
        idx=0
        for g,pg in zip(groups,combo):
            for v in pg:
                perm[idx]=v; idx+=1
        f=tuple(a[perm[i]][perm[j]] for i in range(N) for j in range(N))
        if best is None or f<best: best=f
    return best
def mkey(a):
    return min(canon(a),canon(comp(a)))

print("building G_7/Z_2 via BFS on the iso-class flip graph ...", flush=True)
transitive=[[1 if i<j else 0 for j in range(N)] for i in range(N)]
seen={}; rep={}; edges=set()
k0=mkey(transitive); seen[k0]=0; rep[0]=transitive; q=deque([0])
while q:
    mid=q.popleft(); a=rep[mid]
    for arc in range(M):
        a2=flip(a,arc); k2=mkey(a2)
        if k2 not in seen:
            nid=len(seen); seen[k2]=nid; rep[nid]=a2; q.append(nid)
        j=seen[k2]
        if j!=mid: edges.add((min(mid,j),max(mid,j)))
V=len(seen)
adj=defaultdict(set)
for (u,w) in edges: adj[u].add(w); adj[w].add(u)
print(f"  V = {V} (expect 272), |E| = {len(edges)}", flush=True)
degs=sorted(len(adj[v]) for v in range(V))
print(f"  degree range {degs[0]}..{degs[-1]}, mean {sum(degs)/V:.1f}", flush=True)

# ---- omega (max clique) via Bron-Kerbosch ----
def bk_maxclique(adj,V):
    best=[]
    def bk(R,Pp,X):
        nonlocal best
        if not Pp and not X:
            if len(R)>len(best): best=list(R)
            return
        if len(R)+len(Pp)<=len(best): return
        u=max(Pp|X,key=lambda z:len(adj[z]&Pp)) if (Pp|X) else None
        for v in list(Pp-(adj[u] if u is not None else set())):
            bk(R|{v}, Pp&adj[v], X&adj[v]); Pp=Pp-{v}; X=X|{v}
    bk(set(),set(range(V)),set())
    return best
omega=bk_maxclique(adj,V)
print(f"  omega (max clique) = {len(omega)}: {sorted(omega)[:8]}", flush=True)

# ---- greedy chromatic upper bound (DSATUR) + note ----
def dsatur(adj,V):
    color={}; sat={v:set() for v in range(V)}
    for _ in range(V):
        u=max((v for v in range(V) if v not in color),
              key=lambda v:(len(sat[v]),len(adj[v])))
        c=0
        while c in sat[u]: c+=1
        color[u]=c
        for w in adj[u]: sat[w].add(c)
    return max(color.values())+1
print(f"  DSATUR chromatic upper bound = {dsatur(adj,V)} (repo: chi=6)", flush=True)

# ---- K_t minor via randomized edge contraction + clique check ----
def maxclique_small(nodes,neigh):
    nodes=list(nodes); best=[]
    def bk(R,Pp,X):
        nonlocal best
        if not Pp and not X:
            if len(R)>len(best): best=list(R); return
        if len(R)+len(Pp)<=len(best): return
        for v in list(Pp):
            bk(R+[v],[w for w in Pp if w in neigh[v]],[w for w in X if w in neigh[v]])
            Pp=[w for w in Pp if w!=v]; X=X+[v]
    bk([],nodes,[]); return best
def find_kt_minor(V,adj,t,tries=400):
    for _ in range(tries):
        neigh={v:set(adj[v]) for v in range(V)}
        members={v:{v} for v in range(V)}
        active=set(range(V))
        while len(active)>t:
            us=[u for u in active if neigh[u]]
            if not us: break
            u=random.choice(us); w=random.choice(list(neigh[u]))
            members[u]|=members[w]
            for x in list(neigh[w]):
                neigh[x].discard(w)
                if x!=u: neigh[x].add(u); neigh[u].add(x)
            neigh[u].discard(u); del neigh[w]; active.discard(w)
            if len(active)<=t+3:
                cl=maxclique_small(active,neigh)
                if len(cl)>=t: return [members[c] for c in cl[:t]]
        cl=maxclique_small(active,neigh)
        if len(cl)>=t: return [members[c] for c in cl[:t]]
    return None

def verify_minor(branch,adj):
    # each branch connected in original graph, pairwise adjacent, disjoint
    allv=set()
    for B in branch:
        if allv & B: return False,"not disjoint"
        allv|=B
        # connected
        B=set(B); start=next(iter(B)); st=[start]; vis={start}
        while st:
            x=st.pop()
            for y in adj[x]&B:
                if y not in vis: vis.add(y); st.append(y)
        if vis!=B: return False,"branch not connected"
    for i in range(len(branch)):
        for j in range(i+1,len(branch)):
            if not any((adj[x]&branch[j]) for x in branch[i]):
                return False,f"branch {i},{j} not adjacent"
    return True,"OK"

print("\nsearching for a K_6 minor (randomized contraction)...", flush=True)
k6=find_kt_minor(V,adj,6)
if k6:
    ok,msg=verify_minor(k6,adj)
    sizes=[len(B) for B in k6]
    print(f"  K_6 MINOR FOUND. branch-set sizes {sizes}. certificate: {ok} ({msg})")
    for i,B in enumerate(k6): print(f"    B{i+1} = {sorted(B)}")
else:
    print("  no K_6 minor found in 400 tries (unexpected if chi=6)")

print("\nprobing for a K_7 minor (beyond what chi=6 forces)...", flush=True)
k7=find_kt_minor(V,adj,7,tries=600)
if k7:
    ok,_=verify_minor(k7,adj)
    print(f"  K_7 MINOR FOUND (Hadwiger number >= 7). sizes {[len(B) for B in k7]}, cert {ok}")
else:
    print("  no K_7 minor found in 600 tries (Hadwiger number likely = 6, matching chi)")
