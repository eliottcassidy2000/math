exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'w32_unipotent.py')).read().split("I=np.identity(dim,dtype=np.int64)")[0])
I=np.identity(dim,dtype=np.int64)
p=PR[0]
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
    free=[c for c in range(cols) if c not in piv]
    basis=[]
    for f in free:
        v=np.zeros(cols,dtype=np.int64); v[f]=1
        for i,c in enumerate(piv): v[c]=(-A[i,f])%p
        basis.append(v)
    return basis
groups={}
cnt=0
for L in range(1,5):
    for w in itertools.product(masks,repeat=L):
        M=I.copy()
        for li in w: M=(MI[li]@M)%p
        K=kernel_mod((M-I)%p,p)
        if not K: continue
        cnt+=1
        for v in K:
            nzi=np.nonzero(v)[0][0]; inv=pow(int(v[nzi]),p-2,p); key=tuple(((v*inv)%p).tolist())
            groups.setdefault(key,[]).append(w)
multi={k:v for k,v in groups.items() if len(v)>1}
print('blocks with eigenvalue 1:',cnt,' distinct fixed vectors:',len(groups),' fixed vectors shared by >=2 blocks:',len(multi))
for k,v in list(multi.items())[:6]:
    print('  shared by',len(v),'blocks, e.g.',v[:4], ' support size',sum(1 for x in k if x))
