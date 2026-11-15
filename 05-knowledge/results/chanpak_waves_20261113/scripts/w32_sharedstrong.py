exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'w32_ev1strong.py')).read().split("ev1=0; ev1_strong=[]")[0])
exec(open('fixedvec.py').read().split("groups={}")[0].split("I=np.identity(dim,dtype=np.int64)\np=PR[0]")[1]) if False else None
def kernel_mod(A,p):
    A=A.copy()%p; rows,cols=A.shape; piv=[]; r=0
    for c in range(cols):
        pr=None
        for i in range(r,rows):
            if A[i,c]: pr=i;break
        if pr is None: continue
        A[[r,pr]]=A[[pr,r]]; inv=pow(int(A[r,c]),p-2,p); A[r]=(A[r]*inv)%p
        for i in np.nonzero(A[:,c])[0]:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%p
        piv.append(c); r+=1
    free=[c for c in range(cols) if c not in piv]; basis=[]
    for f in free:
        v=np.zeros(cols,dtype=np.int64); v[f]=1
        for i,c in enumerate(piv): v[c]=(-A[i,f])%p
        basis.append(v)
    return basis
def prim(w):
    l=len(w)
    for d in range(1,l+1):
        if l%d==0 and w==w[:d]*(l//d):
            c=w[:d]; return min(c[i:]+c[:i] for i in range(d))
groups={}
for L in range(1,5):
    for w in itertools.product(masks,repeat=L):
        M=I.copy()
        for li in w: M=(MI[li]@M)%p
        K=kernel_mod((M-I)%p,p)
        if not K or not bulk_strong(list(w)): continue
        for v in K:
            nzi=np.nonzero(v)[0][0]; inv=pow(int(v[nzi]),p-2,p); key=tuple(((v*inv)%p).tolist())
            groups.setdefault(key,set()).add(prim(w))
real={k:v for k,v in groups.items() if len(v)>1}
print('strong eigenvalue-1 fixed vectors:',len(groups),' shared by >=2 distinct primitive strong blocks:',len(real))
for k,v in list(real.items())[:8]: print('  ',sorted(v))
