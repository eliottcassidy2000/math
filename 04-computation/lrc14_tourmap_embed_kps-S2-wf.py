# Just the rotational-embeddability structural answer (no M(S) calls -> fast).
from math import gcd
from itertools import combinations
import importlib.util, os
spec = importlib.util.spec_from_file_location(
    "main", os.path.join(os.path.dirname(__file__), "lrc14_tourmap_danger-interval_kps-S2-wf.py"))
mn = importlib.util.module_from_spec(spec); spec.loader.exec_module(mn)

m=5
def all_classes(m):
    seen={}
    pairs=list(combinations(range(m),2))
    for bits in range(2**len(pairs)):
        adj=[[False]*m for _ in range(m)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: adj[a][b]=True
            else: adj[b][a]=True
        cn=mn.canon(adj,m)
        seen.setdefault(cn,(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m)))
    return seen
allc=all_classes(m)

def rot_adj(res,N):
    mm=len(res); adj=[[False]*mm for _ in range(mm)]
    for i in range(mm):
        for j in range(mm):
            if i==j: continue
            d=(res[i]-res[j])%N
            adj[i][j]=(1<=d<=(N-1)//2)
    return adj

# precompute rotational-embeddable canon keys for odd N up to 21
emb_keys=set()
for N in range(5,22,2):
    for res in combinations(range(N),m):
        emb_keys.add(mn.canon(rot_adj(list(res),N),m))

print("Rotational-embeddability of all 12 n=5 classes (odd N in 5..21):")
for cn,v in sorted(allc.items(),key=lambda kv:kv[1]):
    print(f"  {v}  rotational_embeddable={cn in emb_keys}")
print()
print(f"Total rotationally-embeddable classes: {sum(1 for cn in allc if cn in emb_keys)}/12")
forb=[v for cn,v in allc.items() if cn not in emb_keys]
print("NON-rotational classes:")
for v in sorted(forb): print("  ",v)
