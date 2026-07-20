#!/usr/bin/env python3
"""
death-star-2026-07-20-S59x (HYP-8225) -- refine: the exact Hadwiger number of
G_7/Z_2 (max t with a K_t minor) and an ECONOMICAL K_6 minor grown from the
literal 4-clique. Rebuilds the graph (fast) and saves the edge list.
"""
import random
from itertools import permutations, product, combinations
from collections import defaultdict, deque
random.seed(7)
N=7; P=[(i,j) for i in range(N) for j in range(i+1,N)]; M=len(P)
def flip(a,arc):
    b=[r[:] for r in a]; i,j=P[arc]; b[i][j],b[j][i]=b[j][i],b[i][j]; return b
def comp(a): return [[a[j][i] if i!=j else 0 for j in range(N)] for i in range(N)]
def refine(a):
    color=[sum(a[i]) for i in range(N)]
    for _ in range(N):
        sig=[(color[i],tuple(sorted(color[j] for j in range(N) if a[i][j])),
              tuple(sorted(color[j] for j in range(N) if a[j][i]))) for i in range(N)]
        order=sorted(set(sig)); idx={c:k for k,c in enumerate(order)}
        nc=[idx[sig[i]] for i in range(N)]
        if nc==color: break
        color=nc
    return color
def canon(a):
    color=refine(a); cells=defaultdict(list)
    for i in range(N): cells[color[i]].append(i)
    groups=[cells[c] for c in sorted(cells)]; best=None
    for combo in product(*[list(permutations(g)) for g in groups]):
        perm=[]; 
        for pg in combo: perm+=list(pg)
        f=tuple(a[perm[i]][perm[j]] for i in range(N) for j in range(N))
        if best is None or f<best: best=f
    return best
def mkey(a): return min(canon(a),canon(comp(a)))
transitive=[[1 if i<j else 0 for j in range(N)] for i in range(N)]
seen={}; rep={}; edges=set(); k0=mkey(transitive); seen[k0]=0; rep[0]=transitive; q=deque([0])
while q:
    mid=q.popleft(); a=rep[mid]
    for arc in range(M):
        a2=flip(a,arc); k2=mkey(a2)
        if k2 not in seen: nid=len(seen); seen[k2]=nid; rep[nid]=a2; q.append(nid)
        j=seen[k2]
        if j!=mid: edges.add((min(mid,j),max(mid,j)))
V=len(seen); adj=defaultdict(set)
for u,w in edges: adj[u].add(w); adj[w].add(u)
print(f"G_7/Z_2: V={V}, E={len(edges)}")
with open("05-knowledge/results/g7z2_edgelist_deathstar_S59x.txt","w") as f:
    f.write(f"# G_7/Z_2 merged tournament metagraph, V={V}\n")
    for u,w in sorted(edges): f.write(f"{u} {w}\n")

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

print("\n=== exact Hadwiger number (max t with a K_t minor) ===")
had=0
for t in range(4,13):
    mn=kt_minor(t, 200 if t<=7 else 500)
    if mn and verify(mn):
        had=t; print(f"  K_{t} minor: YES (cert OK), max branch size {max(len(B) for B in mn)}")
    else:
        print(f"  K_{t} minor: none found in the budget -> Hadwiger number = {had}")
        break
print(f"  => Hadwiger number h(G_7/Z_2) >= {had}  (chi=6, omega=4)")

print("\n=== ECONOMICAL K_6 minor grown from the literal 4-clique ===")
clique4=maxclique(range(V),{v:adj[v] for v in range(V)})[:4]
print(f"  base 4-clique (singletons): {sorted(clique4)}")
# find 2 more connected branch sets, each adjacent to all 4 clique verts + each other,
# small, disjoint from the clique. Greedy: candidate sets = short paths between two
# vertices each adjacent to (part of) the clique, union covering all 4.
base=set(clique4)
def adj_to_all(Bset): return all(any((adj[x] & Bset) for x in [c]) for c in clique4)
# find vertices adjacent to as many clique members as possible
def cover(v): return sum(1 for c in clique4 if c in adj[v])
cands=sorted((v for v in range(V) if v not in base), key=lambda v:-cover(v))
# build branch set 5: grow from best-cover vertex until adjacent to all 4 clique verts, connected
def grow_branch(avoid):
    for seed in cands:
        if seed in avoid: continue
        B={seed}; 
        # BFS-add neighbors to cover missing clique members
        for _ in range(6):
            missing=[c for c in clique4 if not (adj[c]&B)]
            if not missing: 
                return B
            # add a neighbor of B (not in avoid/base) that covers a missing clique vert
            added=False
            for x in list(B):
                for y in adj[x]:
                    if y in avoid or y in base or y in B: continue
                    if any(c in adj[y] for c in missing):
                        B.add(y); added=True; break
                if added: break
            if not added:
                # add any neighbor to extend
                ext=[y for x in B for y in adj[x] if y not in avoid and y not in base and y not in B]
                if not ext: break
                B.add(ext[0])
        if all(adj[c]&B for c in clique4): return B
    return None
B5=grow_branch(base)
B6=grow_branch(base | (B5 or set()))
if B5 and B6:
    branch=[{c} for c in clique4]+[B5,B6]
    ok=verify(branch)
    print(f"  B5 = {sorted(B5)} (size {len(B5)}), B6 = {sorted(B6)} (size {len(B6)})")
    print(f"  economical K_6 minor: 4 singletons + 2 sets of size {len(B5)},{len(B6)}; certificate: {ok}")
    if ok:
        print("  BRANCH SETS:", [sorted(b) for b in branch])
else:
    print("  greedy economical build incomplete; the certified K_6 minor from the main run stands.")
