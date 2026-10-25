#!/usr/bin/env python3
"""Exact finite controls, standard library only. No repository writes."""
import itertools, json
from functools import lru_cache
from pathlib import Path

def cross(a,b): return bool(a&b) and not (a<=b or b<=a)
def edges_of(F): return {(i,j) for i in range(len(F)) for j in range(i+1,len(F)) if cross(F[i],F[j])}
def chromatic(n,E):
    if n==0:return 0,[]
    adj=[set() for _ in range(n)]
    for a,b in E:adj[a].add(b);adj[b].add(a)
    order=sorted(range(n),key=lambda a:-len(adj[a]))
    for k in range(1,n+1):
        col=[-1]*n
        def dfs(pos):
            if pos==n:return True
            v=order[pos]; bad={col[w] for w in adj[v]}
            for c in range(k):
                if c not in bad:
                    col[v]=c
                    if dfs(pos+1):return True
            col[v]=-1
            return False
        if dfs(0):return k,col

def direct_cover(F):
    # Independent subset-cover DP uses inclusion/disjointness, no graph calls.
    n=len(F); good=[]
    for mask in range(1,1<<n):
        ids=[i for i in range(n) if mask>>i&1]
        if all(F[i].isdisjoint(F[j]) or F[i].issubset(F[j]) or F[j].issubset(F[i]) for i,j in itertools.combinations(ids,2)):
            good.append(mask)
    @lru_cache(None)
    def dp(mask):
        if not mask:return 0
        low=mask&-mask
        return 1+min(dp(mask^sub) for sub in good if sub&low and sub&mask==sub)
    return dp((1<<n)-1)

def hierarchy(U,F):
    H=set(map(frozenset,F))|{frozenset(U)}|{frozenset([x]) for x in U}
    assert all(A and A<=U for A in H)
    parent={A:min((B for B in H if A<B),key=len) for A in H if A!=U}
    children={B:[A for A in H if parent.get(A)==B] for B in H}
    for B,C in children.items():
        if len(B)>1:
            assert set().union(*C)==B
            assert all(A.isdisjoint(D) for A,D in itertools.combinations(C,2))
            assert len(C)>=2
        else:assert not C
    return len(H)

def check(F,U):
    k,col=chromatic(len(F),edges_of(F))
    assert direct_cover(F)==k
    for c in range(k):hierarchy(U,[F[i] for i in range(len(F)) if col[i]==c])
    return {'clusters':len(F),'crossing_edges':len(edges_of(F)),'minimum_trees':k}

out={}
U=frozenset(itertools.product(range(3),repeat=2))
classes=[[frozenset((x,y) for x,y in U if x==b) for b in range(3)]]
classes += [[frozenset((x,y) for x,y in U if (y-m*x)%3==b) for b in range(3)] for m in range(3)]
F=sum(classes,[])
assert all(cross(F[i],F[j])==(i//3!=j//3) for i,j in itertools.combinations(range(12),2))
out['affine_plane_F3']=check(F,U)
U=frozenset(range(6));F=[frozenset(x for x in U if x%p==a) for p in (2,3) for a in range(p)]
out['CRT_2_by_3']=check(F,U)
P2=F[:2];P3=F[2:]
diagonal={A&B for A in P2 for B in P3 if A&B}
assert diagonal=={frozenset([x]) for x in U}
assert all(A not in diagonal and A!=U for A in F)
out['diagonal_CRT']={'diagonal_level_one_blocks':len(diagonal),'independent_nontrivial_blocks_absent':len(F),'minimum_diagonal_trees':check(list(diagonal)+[U],U)['minimum_trees']}
counts={}
for n in range(6):
    possible=list(itertools.combinations(range(n),2));count=0
    for mask in range(1<<len(possible)):
        E={e for i,e in enumerate(possible) if mask>>i&1}
        F=[frozenset([('v',v)]+[('e',a,b) for a,b in E if v in (a,b)]) for v in range(n)]
        assert edges_of(F)==E
        assert all(not A<=B for i,A in enumerate(F) for j,B in enumerate(F) if i!=j)
        U=frozenset().union(*F) if F else frozenset([('dummy',)])
        check(F,U);count+=1
    counts[str(n)]=count
out['all_labelled_graph_realizations_n_le_5']=counts
U=frozenset(range(3));subsets=[frozenset(x for x in U if mask>>x&1) for mask in range(1,8)]
for mask in range(128):check([S for i,S in enumerate(subsets) if mask>>i&1],U)
out['all_cluster_families_on_3_labels']=128
# Clique bound is not always tight: C5 realizes omega=2, chi=3.
E={(0,1),(1,2),(2,3),(3,4),(0,4)}
F=[frozenset([('v',v)]+[('e',a,b) for a,b in E if v in (a,b)]) for v in range(5)]
out['odd_cycle_hostile']=check(F,frozenset().union(*F))
assert out['odd_cycle_hostile']['minimum_trees']==3
out['status']='PASS'
print(json.dumps(out,indent=2,sort_keys=True))
