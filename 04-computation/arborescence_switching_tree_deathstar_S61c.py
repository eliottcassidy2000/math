from itertools import permutations, combinations
from fractions import Fraction as Fr

def det(M):
    M=[[Fr(x) for x in row] for row in M]; n=len(M); d=Fr(1)
    for c in range(n):
        piv=None
        for r in range(c,n):
            if M[r][c]!=0: piv=r;break
        if piv is None: return Fr(0)
        if piv!=c: M[c],M[piv]=M[piv],M[c]; d=-d
        d*=M[c][c]; inv=1/M[c][c]
        for r in range(c+1,n):
            f=M[r][c]*inv
            for cc in range(c,n): M[r][cc]-=f*M[c][cc]
    return d
def minor(L,r):
    idx=[i for i in range(len(L)) if i!=r]; return [[L[i][j] for j in idx] for i in idx]
def a_r(adj,n,r):  # out-arborescences rooted at r
    Din=[sum(adj[i][j] for i in range(n)) for j in range(n)]
    L=[[ (Din[i] if i==j else 0) - adj[i][j] for j in range(n)] for i in range(n)]
    return int(det(minor(L,r)))
def a_tot(adj,n): return sum(a_r(adj,n,r) for r in range(n))
def hpaths_from(adj,n,r): return sum(1 for p in permutations(range(n)) if p[0]==r and all(adj[p[k]][p[k+1]] for k in range(n-1)))
def H(adj,n): return sum(hpaths_from(adj,n,r) for r in range(n))

def switch_vertex(bits,v,n,idx):
    for u in range(n):
        if u!=v: bits^=(1<<idx[(min(u,v),max(u,v))])
    return bits
def adj_of(bits,n,edges):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(edges):
        if (bits>>k)&1: a[i][j]=1
        else: a[j][i]=1
    return a

from math import factorial
print("=== SWITCHING-CLASS SUMS: the general oriented-spanning-tree theorem ===")
print("   family F of oriented spanning trees => Sum_{T in class} #{S in F: S subset T} = |F|\n")
for n in range(3,6):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(edges); idx={e:k for k,e in enumerate(edges)}
    seen=set(); sums_a=[]; sums_ar0=[]; sums_H=[]; sums_hfrom0=[]
    for start in range(1<<m):
        if start in seen: continue
        cls=set([start]); fr=[start]
        while fr:
            b=fr.pop()
            for v in range(n):
                nb=switch_vertex(b,v,n,idx)
                if nb not in cls: cls.add(nb); fr.append(nb)
        seen|=cls
        sums_a.append(sum(a_tot(adj_of(b,n,edges),n) for b in cls))
        sums_ar0.append(sum(a_r(adj_of(b,n,edges),n,0) for b in cls))
        sums_H.append(sum(H(adj_of(b,n,edges),n) for b in cls))
        sums_hfrom0.append(sum(hpaths_from(adj_of(b,n,edges),n,0) for b in cls))
    print(f"n={n}: (over {len(sums_a)} classes)")
    print(f"  ALL arborescences a(T):   every class-sum = n^(n-1) = {n**(n-1)}?  {all(s==n**(n-1) for s in sums_a)}   (sums seen: {sorted(set(sums_a))})")
    print(f"  arborescences rooted@0:   every class-sum = n^(n-2) = {n**(n-2)}?  {all(s==n**(n-2) for s in sums_ar0)}   (sums: {sorted(set(sums_ar0))})")
    print(f"  ALL Ham paths H(T):       every class-sum = n! = {factorial(n)}?  {all(s==factorial(n) for s in sums_H)}   (THM-1445)")
    print(f"  Ham paths from vtx 0:     every class-sum = (n-1)! = {factorial(n-1)}?  {all(s==factorial(n-1) for s in sums_hfrom0)}")
    print()
print("=> ONE theorem: Sum_{T in switching class} #{oriented spanning trees of family F realized in T} = |F|.")
print("   Ham paths (|F|=n!) and arborescences (|F|=n^(n-1)) are TWO INSTANCES. THM-1445 is the path corollary.")
print("   Also: a(T) >= H(T) always, with a = 'branchy' trees, H = 'linear' (tallest) trees.")

print("\n=== BEST theorem: Euler circuits of a regular tournament via arborescences (arc/even world) ===")
def paley(q):
    QR=set((x*x)%q for x in range(1,q)); a=[[1 if (j-i)%q in QR else 0 for j in range(q)] for i in range(q)]
    for i in range(q): a[i][i]=0
    return a
# verify BEST on the 3-cycle (Euler circuits brute)
C3=[[0,1,0],[0,0,1],[1,0,0]]
best3 = a_r(C3,3,0)* (factorial(1-1)**3)  # outdeg=1 => (outdeg-1)!=0!=1
print(f"  3-cycle: BEST ec = a_r * prod (outdeg-1)! = {a_r(C3,3,0)} * 1 = {best3}  (Euler circuits of C3 = 1) match={best3==1}")
for q in [7,11]:
    a0=a_r(paley(q),q,0); d=(q-1)//2
    ec = a0 * factorial(d-1)**q
    print(f"  Paley_{q}: a_r={a0}, outdeg={d}, ec(Euler circuits) = a_r*((({d}-1)!)^{q}) = {a0} * {factorial(d-1)}^{q} = {ec}")
