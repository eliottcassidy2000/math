# lrc14_tourmap_M3_characterize_kps-S2-wf.py
# Characterize the M3 lap-count-parity map precisely.
# arc(i,j) = (vi>vj) XOR (floor(|vi-vj|*tau) is odd).
# Questions:
#  1. Confirm realized H-set at n=3,4,5 over ALL free tau (intrinsic).
#  2. Is M3 EVER non-transitive in a way loneliness changes? (already: no.)
#  3. WHY is the forbidden set what it is? Describe the realized classes structurally.
#  4. Sanity: M3 with the parity twist REMOVED reduces to M1a (transitive)?  yes by construction.
#  5. Does M3 depend on tau only through the parity vector p_{ij}=floor(|vi-vj|*tau) mod 2?
#     If so, M3 is a function of the speeds and a Z/2 "parity pattern" -> characterize which
#     parity patterns are achievable, and which tournaments they yield.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from functools import reduce, lru_cache

def primitive(S): return reduce(gcd,S)==1
def gen_speed_sets(n,maxspeed):
    return [tuple(sorted(S)) for S in combinations(range(1,maxspeed+1),n) if primitive(S)]

_PERMS={n:list(permutations(range(n))) for n in (3,4,5)}
@lru_cache(maxsize=None)
def analyze_up(n,up):
    pos={}; idx=0
    for a in range(n):
        for b in range(a+1,n): pos[(a,b)]=up[idx]; idx+=1
    def arc(x,y): return pos[(x,y)]==1 if x<y else pos[(y,x)]==0
    n3=0
    for a,b,c in combinations(range(n),3):
        if arc(a,b) and arc(b,c) and arc(c,a): n3+=1
        if arc(a,c) and arc(c,b) and arc(b,a): n3+=1
    H=0
    for perm in _PERMS[n]:
        if all(arc(perm[k],perm[k+1]) for k in range(n-1)): H+=1
    best=None
    for perm in _PERMS[n]:
        bits=[]
        for a in range(n):
            for b in range(a+1,n):
                pa,pb=perm[a],perm[b]
                bits.append(pos[(pa,pb)] if pa<pb else 1-pos[(pb,pa)])
        key=tuple(bits)
        if best is None or key<best: best=key
    sc=[0]*n
    for a in range(n):
        for b in range(a+1,n):
            if pos[(a,b)]==1: sc[a]+=1
            else: sc[b]+=1
    return best,n3,H,tuple(sorted(sc))
def classify_up(n,up): return analyze_up(n,up)

def M3_up(verts,tau):
    n=len(verts); up=[]
    for a in range(n):
        for b in range(a+1,n):
            vi,vj=verts[a],verts[b]
            d=abs(vi-vj); k=int(d*tau)
            bit = (1 if vi>vj else 0)^(k%2)  # but verts sorted ascending so vi<vj always -> vi>vj false
            up.append(bit)
    return tuple(up)

def all_classes(n):
    res={}
    for bits in range(2**(n*(n-1)//2)):
        up=tuple((bits>>k)&1 for k in range(n*(n-1)//2))
        c,n3,H,sc=analyze_up(n,up)
        if c not in res: res[c]=(n3,H,sc)
    return res
FREE={n:all_classes(n) for n in (3,4,5)}

print("M3 characterization  (kind-pasteur-S2-wf)")
print("="*70)

# Realized H-set over ALL free tau, broad
for n in (3,4,5):
    maxspeed={3:20,4:16,5:13}[n]
    SS=gen_speed_sets(n,maxspeed)
    D=60
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
    realized={}
    for S in SS:
        for tau in free_taus:
            up=M3_up(S,tau)
            c,n3,H,sc=analyze_up(n,up)
            if c not in realized: realized[c]=(n3,H,sc)
    Hs=sorted(set(v[1] for v in realized.values()))
    print(f"\nn={n}: realized {len(realized)}/{len(FREE[n])}; realized H-values={Hs}")
    # list realized classes
    for c in sorted(realized,key=lambda c:(realized[c][1],realized[c][0])):
        n3,H,sc=realized[c]; print(f"   realized: H={H} c3={n3} score={sc}")
    print("   FORBIDDEN:")
    for c in sorted(set(FREE[n])-set(realized),key=lambda c:(FREE[n][c][1],FREE[n][c][0])):
        n3,H,sc=FREE[n][c]; print(f"     H={H} c3={n3} score={sc}")

# Does M3 depend on tau only via parity pattern? Check: for n=5, enumerate all achievable
# parity vectors p_{ab}=floor(|v_a-v_b|*tau) mod 2 over free tau, for a few speed sets, and
# show the resulting tournament depends only on (sorted speeds, parity vector).
print("\n" + "="*70)
print("M3 parity-pattern structure (n=5 example S=(1,2,3,4,5))")
S=(1,2,3,4,5)
D=80
free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
seen=set()
patterns={}
for tau in free_taus:
    n=len(S)
    pv=[]
    for a in range(n):
        for b in range(a+1,n):
            d=abs(S[a]-S[b]); pv.append(int(d*tau)%2)
    pv=tuple(pv)
    up=M3_up(S,tau)
    c,n3,H,sc=analyze_up(n,up)
    if pv not in patterns: patterns[pv]=(c,n3,H,sc)
    seen.add((c,H))
print(f"  #distinct parity vectors over free tau: {len(patterns)}")
print(f"  #distinct (class,H) realized for this S: {len(seen)}")
hh=sorted(set(h for _,h in seen))
print(f"  H values for S=(1,2,3,4,5): {hh}")
print("\nDONE")
