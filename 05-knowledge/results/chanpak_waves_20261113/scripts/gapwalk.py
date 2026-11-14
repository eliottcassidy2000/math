"""Wave 24, P64: second-largest eigenvalue modulus of the block walk on F_p^3 minus 0.
Blocks: the 16 invertible blocks J w with |w| <= 5 (w empty or O v Z), uniform weights."""
import itertools, numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, sys
def step(s,m):
    H,E,L,h=s
    if m=='O': return (H,0,2*L,h)
    if m=='Z': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def blockmat(w):
    cols=[]
    for e in ((1,0,0),(0,1,0),(0,0,1)):
        s=step((e[0],e[1],e[2],0),'J')
        for ch in w: s=step(s,ch)
        cols.append(s[:3])
    return np.array([[cols[c][r] for c in range(3)] for r in range(3)],dtype=np.int64)
ws=['']+['O'+''.join(v)+'Z' for n in range(0,4) for v in itertools.product('OZ',repeat=n)]
Ms=[blockmat(w) for w in ws]; assert len(Ms)==16
def gap(p):
    idx=np.arange(1,p**3)                       # nonzero vectors, encoded x + p y + p^2 z
    V=np.stack([idx%p,(idx//p)%p,idx//(p*p)])  # 3 x n
    n=len(idx); rows=[];cols=[]
    for M in Ms:
        W=(M%p)@V%p; code=W[0]+p*W[1]+p*p*W[2]
        rows.append(np.arange(n)); cols.append(code-1)
    P=sps.csr_matrix((np.full(16*n,1/16),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n))
    k=min(6,n-2)
    vals=spl.eigs(P,k=k,which='LM',return_eigenvectors=False,maxiter=5000,tol=1e-9)
    mods=sorted(np.abs(vals),reverse=True)
    return mods
print("p   |lambda_1|  |lambda_2|  |lambda_3|")
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43]:
    m=gap(p); print("%-3d %.6f   %.6f   %.6f"%(p,m[0],m[1],m[2])); sys.stdout.flush()
