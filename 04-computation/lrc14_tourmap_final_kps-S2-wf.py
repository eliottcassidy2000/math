# Final confirmation: identify the structural nature of the forbidden class and
# check M3 robustness. Also test M2 (control: realizes everything) to confirm contrast.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import importlib.util, os

spec = importlib.util.spec_from_file_location(
    "main", os.path.join(os.path.dirname(__file__), "lrc14_tourmap_danger-interval_kps-S2-wf.py"))
mn = importlib.util.module_from_spec(spec); spec.loader.exec_module(mn)

def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

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

m=5
allc=all_classes(m)

# Is the forbidden class (15,4,(1,2,2,2,3)) ROTATIONALLY EMBEDDABLE in Z/14? in Z/p any p?
# A tournament is "rotational"/"locally transitive" iff it is a circulant tournament.
# At n=5 the only circulant (rotational) tournament is the QR tournament = regular (2,2,2,2,2),H=15,c3=5.
# Check: is (15,4) class an induced subtournament of ANY rotational tournament Z/N i->j iff (i-j)%N in 1..(N-1)/2 ?
def rot_adj(res,N):
    m=len(res); adj=[[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            d=(res[i]-res[j])%N
            adj[i][j] = (1<=d<= (N-1)//2)
    return adj

# find canon key of forbidden class
forb_key=None
for cn,v in allc.items():
    if v==(15,4,(1,2,2,2,3)): forb_key=cn
print("Forbidden class canon found:", forb_key is not None)

print("\nIs (15,4,(1,2,2,2,3)) an induced subtournament of rotational Z/N for odd N up to 17?")
embeddable=False
for N in range(5,18,2):
    for res in combinations(range(N),m):
        adj=rot_adj(list(res),N)
        if mn.canon(adj,m)==forb_key:
            embeddable=True
            print(f"  EMBEDS in Z/{N} at residues {res}")
            break
    if embeddable: break
if not embeddable:
    print("  NOT embeddable in any rotational Z/N (N odd <=17): it is a NON-rotational tournament.")

# Which n=5 classes ARE rotationally embeddable (in any odd Z/N)?
print("\nRotational-embeddability of all 12 classes (odd N<=17):")
for cn,v in sorted(allc.items(),key=lambda kv:kv[1]):
    emb=False
    for N in range(5,18,2):
        for res in combinations(range(N),m):
            if mn.canon(rot_adj(list(res),N),m)==cn: emb=True;break
        if emb:break
    print(f"  {v}  rotational_embeddable={emb}")

# M3 robustness
print("\n[M3 signed danger-arrival @ optimum] window scan:")
for vmax in [12,18,24,30]:
    real={}
    for S in (c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1):
        S=list(S); gap,tau0=mn.M(S)
        adj=mn.method3_adj(S,tau0); real.setdefault(mn.canon(adj,m),tau0)
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}")
    if vmax==30:
        for cn in sorted(forb,key=lambda c:allc[c]): print("     FORBIDDEN",allc[cn])
