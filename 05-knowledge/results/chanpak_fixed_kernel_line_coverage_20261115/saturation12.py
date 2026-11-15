# Saturation of the K1 junction at order 12: empty lines (N odd, N not 0 mod 3) above the order-10 window top,
# as targeted samples accumulate; then the certified continuous range for all odd N.
import numpy as np, glob
V10=np.loadtxt('pairs_10.txt',dtype=np.int64); v10=np.unique(3*V10[:,0]+V10[:,1])
S11=np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s11_*.txt'))]); v11=np.unique(3*S11[:,0]+S11[:,1])
base=[np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s12_*.txt'))]
targ=[np.loadtxt(f,dtype=np.int64,usecols=(0,1)).reshape(-1,2) for f in sorted(glob.glob('j12_*.txt'))]
lo,hi=241521,400000
N=np.arange(lo,hi+1); N=N[(N%2==1)&(N%3!=0)]
def empties(arrs):
    s=np.unique(np.concatenate([3*a[:,0]+a[:,1] for a in arrs]))
    return int((~np.isin(N,s)).sum()), s
print('lines N odd, N not 0 mod 3, in [%d, %d]: %d'%(lo,hi,len(N)))
for frac in (0,0.25,0.5,0.75,1.0):
    arrs=base+[a[:int(len(a)*frac)] for a in targ]
    e,_=empties(arrs)
    print('  base 4.8M + %3d%% of targeted run (%d in-window-filtered samples): empty lines %d'%(100*frac,sum(int(len(a)*frac) for a in targ),e))
e,s12=empties(base+targ)
U=np.unique(np.concatenate([v10,v11,s12]))
M=np.arange(36835,U[-1]+1); M=M[M%2==1]; miss=M[~np.isin(M,U)]
print('certified (K1, orders 10 exact + 11, 12 witnessed): every odd N in [36835, %d] is 3H(A)+P2(A) for some A of order 10-12'%(miss[0]-2 if len(miss) else U[-1]))
