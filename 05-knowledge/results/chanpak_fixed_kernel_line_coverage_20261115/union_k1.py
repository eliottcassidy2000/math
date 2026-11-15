import numpy as np, glob
def load(f):
    L=open(f).read().split('\n')
    L=[l for l in L if len(l.split())==3]
    return np.array([list(map(int,l.split()[:2])) for l in L],dtype=np.int64).reshape(-1,2)
V={n:np.loadtxt('pairs_%d.txt'%n,dtype=np.int64).reshape(-1,2) for n in range(2,11)}; V[1]=np.array([[1,0]])
R=np.load('reducible_hf.npy',allow_pickle=True)          # (h, F2) at orders 11..14
S11=np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s11_*.txt'))])
S12=np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s12_*.txt'))]+[load(f) for f in sorted(glob.glob('j12_*.txt'))])
def k1(P): return P[:,0]*3+P[:,1]                        # from (H, P2)
exact=np.unique(np.concatenate([k1(V[n]) for n in V]+[r[:,0]+r[:,1] for r in R]))      # h + F2 for reducible
witn=np.unique(np.concatenate([k1(S11),k1(S12), S11[:,0]+2*(2*S11[:,0]+S11[:,1])]))  # sampled, and sampled-11 => K1 : (h, 2F2)
for name,U in (('exact only (census <= 10, reducible 11-14)',exact),('exact + witnessed samples',np.unique(np.concatenate([exact,witn])))):
    M=np.arange(36835,1400001); M=M[M%2==1]; miss=M[~np.isin(M,U)]
    print('%s: every odd N in [36835, %d] covered; misses in (that, 1.4e6]: %d'%(name,(miss[0]-2) if len(miss) else 1400000,len(miss)))
    for r in (0,1,2):
        mr=miss[miss%3==r]
        print('     N = %d mod 3: first miss %s, misses up to 1.4e6: %d'%(r,mr[0] if len(mr) else None,len(mr)))
