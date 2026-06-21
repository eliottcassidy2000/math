"""
Third pass: is P PREDICTIVE relative to the OCF invariants the repo actually uses?
- Compute H (Hamiltonian path count = Redei) via DP, alpha_1=c3, and check:
  (Q1) Does P distinguish tournaments that H + c3 do NOT? (adds discriminating power)
  (Q2) Conversely is P redundant given H? (is P a function of H?)
  (Q3) Confirm specific claimed witnesses exist literally.
- Also test 'parallel-oriented' interpretation invariance: the value P depends on the
  convex embedding choice (vertices placed at 1..n in label order). Is P a graph invariant
  (relabeling-invariant) or an artifact of the fixed staircase labeling? Test by relabeling.
"""
from itertools import combinations, product, permutations
from collections import defaultdict

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A
def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t
def Pinv(A,n):
    P=0
    for w,x,y,z in combinations(range(1,n+1),4):
        if A[w][y]==A[x][z]: P+=1
    return P
def H_redei(A,n):
    # count Hamiltonian paths (directed) = sum over all perms with all consecutive arcs present
    # DP over subsets: paths[mask][last]
    full=(1<<n)-1
    paths=[[0]*(n) for _ in range(1<<n)]
    for v in range(n): paths[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            cur=paths[mask][last]
            if not cur: continue
            for nxt in range(n):
                if mask&(1<<nxt): continue
                if A[last+1][nxt+1]:  # arc last->nxt
                    paths[mask|(1<<nxt)][nxt]+=cur
    return sum(paths[full][v] for v in range(n))

out=[]
def log(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)

# Q1/Q2: joint distribution of (H, c3, P) over full enumeration n<=6
def joint(n):
    T=tiles(n); F=len(T)
    triples=defaultdict(int)
    pairs_Hc3=set(); full=set()
    H_to_P=defaultdict(set)
    Hc3_to_P=defaultdict(set)
    for bits in product([0,1],repeat=F):
        A=adj(n,list(bits),T)
        h=H_redei(A,n); c=c3(A,n); p=Pinv(A,n)
        pairs_Hc3.add((h,c)); full.add((h,c,p))
        H_to_P[h].add(p)
        Hc3_to_P[(h,c)].add(p)
    P_is_func_of_H = all(len(v)==1 for v in H_to_P.values())
    P_is_func_of_Hc3 = all(len(v)==1 for v in Hc3_to_P.values())
    P_adds_beyond_Hc3 = len(full) > len(pairs_Hc3)
    return P_is_func_of_H, P_is_func_of_Hc3, P_adds_beyond_Hc3, len(pairs_Hc3), len(full)

log("=== (H,c3) vs (H,c3,P) discriminating power ===")
for n in (4,5,6):
    a,b,c,np_,nf=joint(n)
    log(f"n={n}: P func of H? {a} ; P func of (H,c3)? {b} ; "
        f"P adds beyond (H,c3)? {c} ; #(H,c3)={np_} #(H,c3,P)={nf}")

# Q3: relabeling invariance — is P invariant under vertex relabeling (graph invariant)
# or an artifact of the convex order 1..n? Take a fixed tournament, permute labels, recompute P.
def relabel_test(n, bits):
    T=tiles(n); A=adj(n,list(bits),T)
    base=Pinv(A,n)
    vals=set()
    cnt=0
    for perm in permutations(range(1,n+1)):
        # relabel: new A'[perm[i-1]][perm[j-1]] = A[i][j]
        B=[[0]*(n+1) for _ in range(n+1)]
        for i in range(1,n+1):
            for j in range(1,n+1):
                if i!=j: B[perm[i-1]][perm[j-1]]=A[i][j]
        vals.add(Pinv(B,n))
        cnt+=1
        if cnt>=200: break
    return base, sorted(vals)

log("\n=== P relabeling-(in)variance (is P a graph invariant or embedding artifact?) ===")
import random; random.seed(1)
for n in (5,6):
    F=len(tiles(n))
    bits=[random.randint(0,1) for _ in range(F)]
    base,vals=relabel_test(n,bits)
    log(f"n={n}: sample tournament P(convex order)={base} ; P over relabelings = {vals}  "
        f"({'INVARIANT' if len(vals)==1 else 'DEPENDS ON EMBEDDING -> NOT a graph invariant'})")

with open("05-knowledge/results/conn_verify3_crossingroad_kps-Sx-wf.out","w") as f:
    f.write("\n".join(out))
