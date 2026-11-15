import numpy as np, glob
def load(f):
    L=[l for l in open(f).read().split('\n') if len(l.split())==3]
    return np.array([list(map(int,l.split()[:2])) for l in L],dtype=np.int64).reshape(-1,2)
V={n:np.loadtxt('pairs_%d.txt'%n,dtype=np.int64).reshape(-1,2) for n in range(2,11)}; V[1]=np.array([[1,0]])
toHF=lambda P: np.stack([P[:,0],2*P[:,0]+P[:,1]],1)
R=np.load('reducible_hf.npy',allow_pickle=True)
S11=toHF(np.unique(np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s11_*.txt'))]),axis=0))
S12=toHF(np.unique(np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s12_*.txt'))]+[load(f) for f in sorted(glob.glob('j12_*.txt'))]),axis=0))
small=[toHF(V[n]) for n in (1,2,3,4)]
def prod(A,B):
    return np.stack([(A[:,0][:,None]*B[:,0][None,:]).ravel(),(A[:,1][:,None]*B[:,1][None,:]).ravel()],1)
exactHF=np.concatenate([toHF(V[n]) for n in V]+list(R))
witHF=np.concatenate([S11,S12]+[prod(S,Q) for S in (S11,S12) for Q in small])
allHF=np.concatenate([exactHF,witHF])
print('pairs: exact %d, witnessed (incl. products with blocks of order <= 4) %d'%(len(exactHF),len(witHF)))
for name,stat in (('K1: N = H + F2',lambda P:P[:,0]+P[:,1]),('TT2: N = 3F2 - H',lambda P:3*P[:,1]-P[:,0])):
    for lab,P in (('exact',exactHF),('exact+witnessed',allHF)):
        U=np.unique(stat(P)); start=36835 if name.startswith('K1') else 188269
        M=np.arange(start,4000001); M=M[M%2==1]; miss=M[~np.isin(M,U)]
        cert=(miss[0]-2) if len(miss) else 4000000
        print('%-18s %-16s every odd N in [%d, %d] covered; misses up to 4e6: %d (by class mod 3: %s)'
              %(name,lab,start,cert,len(miss),[int((miss%3==r).sum()) for r in (0,1,2)]))
