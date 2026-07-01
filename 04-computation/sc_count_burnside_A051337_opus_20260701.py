"""#SC(n) via Burnside over anti-automorphisms. fix_comp(pi)=#tournaments T with pi(T)=complement(T)."""
from itertools import combinations
from math import factorial
from collections import Counter
def fix_comp(perm, n):
    pairs=list(combinations(range(n),2)); pidx={p:i for i,p in enumerate(pairs)}; P=len(pairs)
    inv=[0]*n
    for i,pi in enumerate(perm): inv[pi]=i
    parent=list(range(P)); rank=[0]*P; par=[0]*P   # par[x] = parity from x to parent[x]
    def find(x):
        if parent[x]==x: return x,0
        r,p=find(parent[x]); parent[x]=r; par[x]^=p; return r,par[x]
    def union(a,b,rel):
        ra,pa=find(a); rb,pb=find(b)
        if ra==rb: return (pa^pb)==rel
        parent[ra]=rb; par[ra]=pa^pb^rel; return True
    for i,(u,v) in enumerate(pairs):
        a,c=inv[u],inv[v]
        if a<c: flip=1; qq=pidx[(a,c)]
        else: flip=0; qq=pidx[(c,a)]
        if not union(i,qq,flip): return 0
    return 2**len(set(find(i)[0] for i in range(P)))
def partitions(n,m=None):
    if m is None: m=n
    if n==0: yield []; return
    for k in range(min(n,m),0,-1):
        for rest in partitions(n-k,k): yield [k]+rest
def rep_perm(part,n):
    perm=[0]*n; s=0
    for c in part:
        for i in range(c): perm[s+i]=s+(i+1)%c
        s+=c
    return perm
def class_size(part,n):
    cnt=Counter(part); denom=1
    for c,mult in cnt.items(): denom*=(c**mult)*factorial(mult)
    return factorial(n)//denom
def SC(n):
    return sum(class_size(p,n)*fix_comp(rep_perm(p,n),n) for p in partitions(n))//factorial(n)
A={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:9733056,11:903753248,12:154108311168}
scv={}
print(f"{'n':>3} {'#SC':>12} {'A(n-1)':>10} {'even:=A(n-1)':>13}")
for n in range(3,15):
    s=SC(n); scv[n]=s
    tag=('YES' if (n%2==0 and s==A.get(n-1)) else '')
    print(f"{n:>3} {s:>12} {str(A.get(n-1)):>10} {tag:>13}")
print("\nODD  #SC:", {n:scv[n] for n in scv if n%2==1})
print("EVEN #SC:", {n:scv[n] for n in scv if n%2==0})
# odd rule probes
print("\nOdd-n probes (n=5,7,9,11,13):")
for n in [5,7,9,11,13]:
    if n in scv:
        half=(n-1)//2
        print(f"  n={n}: #SC={scv[n]}  ; A(n-2)={A.get(n-2)} #SC/#SC(n-2)={scv[n]/scv.get(n-2,1):.4f} ; #SC(n)-2*A(n-2)={scv[n]-2*(A.get(n-2) or 0)}")

print("\n=== CORRECTION check: even rule #SC(2k)=A000568(2k-1)? ===")
for n in [4,6,8,10]:
    print(f"  n={n}: #SC={scv.get(n)} vs A000568({n-1})={A.get(n-1)}  equal? {scv.get(n)==A.get(n-1)}")
print("\n=== contributing cycle types (fix_comp>0) -- the valid anti-automorphisms ===")
for n in [4,5,6,7,8,9]:
    contrib=[]
    for p in partitions(n):
        fc=fix_comp(rep_perm(p,n),n)
        if fc>0: contrib.append((tuple(p),fc))
    print(f"  n={n}: {len(contrib)} contributing types; e.g. {contrib[:5]}")
    # characterize: cycle lengths mod 4
    modok=all(all((c%4 in (0,)) or (c in (1,2)) for c in p) for p,_ in contrib)
    print(f"      all cycles length divisible by 4 OR in {{1,2}}? {modok}")
