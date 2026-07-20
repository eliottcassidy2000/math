from itertools import permutations, combinations, product
from fractions import Fraction as Fr

def ham_paths_from(adj, n, r=None):
    cnt=0
    for p in permutations(range(n)):
        if r is not None and p[0]!=r: continue
        if all(adj[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
    return cnt
def H(adj,n): return ham_paths_from(adj,n)

def det(M):
    import copy
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
    idx=[i for i in range(len(L)) if i!=r]
    return [[L[i][j] for j in idx] for i in idx]

def matrix_tree_in(adj,n,r):  # minor of (D_out - A)
    Dout=[sum(adj[i]) for i in range(n)]
    L=[[ (Dout[i] if i==j else 0) - adj[i][j] for j in range(n)] for i in range(n)]
    return int(det(minor(L,r)))
def matrix_tree_out(adj,n,r): # minor of (D_in - A)
    Din=[sum(adj[i][j] for i in range(n)) for j in range(n)]  # in-deg of j
    # careful: Din[j] = number of i with adj[i][j]=1
    Din=[sum(adj[i][j] for i in range(n)) for j in range(n)]
    L=[[ (Din[i] if i==j else 0) - adj[i][j] for j in range(n)] for i in range(n)]
    return int(det(minor(L,r)))

def brute_out_arb(adj,n,r):
    # count spanning trees (edge sets) that are out-arborescences rooted at r
    edges=[(i,j) for i in range(n) for j in range(n) if i!=j and adj[i][j]]
    cnt=0
    for treearcs in combinations(edges,n-1):
        # must be arborescence rooted at r: every non-root has in-deg 1, root in-deg 0, connected/acyclic
        indeg={v:0 for v in range(n)}
        for (i,j) in treearcs: indeg[j]+=1
        if indeg[r]!=0: continue
        if any(indeg[v]!=1 for v in range(n) if v!=r): continue
        # check reachability from r (then it's a tree)
        adj2={v:[] for v in range(n)}
        for (i,j) in treearcs: adj2[i].append(j)
        seen={r}; st=[r]
        while st:
            x=st.pop()
            for y in adj2[x]:
                if y not in seen: seen.add(y); st.append(y)
        if len(seen)==n: cnt+=1
    return cnt

def tournaments(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(edges)
    for bits in range(1<<m):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj,bits

print("=== calibrate Matrix-Tree vs brute-force out-arborescences (n=4, sample) ===")
n=4; ok_out=True; ok_in_is_out=True
for adj,bits in list(tournaments(n))[:20]:
    for r in range(n):
        b=brute_out_arb(adj,n,r)
        mo=matrix_tree_out(adj,n,r)
        mi=matrix_tree_in(adj,n,r)
        if b!=mo: ok_out=False
    if not ok_out: break
print("  matrix_tree_out == brute out-arborescence:", ok_out)

print("\n=== a_r(T), H(T), Ham-paths-from-r for n=3,4 (a few) ===")
for n in [3,4]:
    print(f" n={n}:")
    for adj,bits in list(tournaments(n))[:6]:
        ar=[matrix_tree_out(adj,n,r) for r in range(n)]
        hpr=[ham_paths_from(adj,n,r) for r in range(n)]
        print(f"   bits={bits:2d}  a_r={ar} tot={sum(ar)}  Hpaths_from_r={hpr} H={sum(hpr)}  (a_r>=Hfrom_r: {all(ar[r]>=hpr[r] for r in range(n))})")

print("\n=== Paley closed form a_r(Paley_q) = (1/q)[q(q+1)/4]^((q-1)/2) ===")
def paley(q):  # q prime, q=3 mod 4
    QR=set((x*x)%q for x in range(1,q))
    adj=[[1 if (j-i)%q in QR else 0 for j in range(q)] for i in range(q)]
    for i in range(q): adj[i][i]=0
    return adj
for q in [3,7,11]:
    adj=paley(q)
    ar0=matrix_tree_out(adj,q,0)
    closed=Fr(1,q)*Fr(q*(q+1),4)**((q-1)//2)
    print(f"   q={q}: a_0(Paley) = {ar0}   closed form (1/q)[q(q+1)/4]^((q-1)/2) = {closed}   match={ar0==closed}")
