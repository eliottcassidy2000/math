import numpy as np, sys
from math import gcd
V={n:np.loadtxt('pairs_%d.txt'%n,dtype=np.int64).reshape(-1,2) for n in range(2,11)}
V[1]=np.array([[1,0]],dtype=np.int64)
# readers B: orders 1..4 distinct pairs -> kernels (alpha,beta) = (X+2Y, X), X = h+P2, Y = h
readers=[]
for n in range(1,5):
    for h,P in V[n]:
        X,Y=h+P,h; readers.append((n,int(h),int(P),int(X+2*Y),int(X)))
print('readers (order, h, P2) -> kernel (alpha, beta):')
for r in readers: print('   ',r)
def full_window(vals, mod_ok):
    """vals: sorted unique ints. mod_ok: function giving admissible mask. Return longest run of consecutive admissible N all covered."""
    lo,hi=vals[0],vals[-1]
    N=np.arange(lo,hi+1,dtype=np.int64)
    adm=mod_ok(N)
    hit=np.zeros(len(N),dtype=bool); hit[vals-lo]=True
    bad=adm & ~hit
    # longest stretch without admissible misses
    idx=np.nonzero(bad)[0]
    edges=np.concatenate(([-1],idx,[len(N)]))
    gaps=np.diff(edges)-1
    k=int(np.argmax(gaps)); a=edges[k]+1; b=edges[k+1]-1
    return int(N[a]),int(N[b]),int(adm.sum()),int(bad.sum())
for (ro,rh,rP,al,be) in readers[:6]:
    print('\nkernel alpha=%d beta=%d  (reader of order %d, h=%d, P2=%d)'%(al,be,ro,rh,rP))
    allv=[]
    for n in range(5,11):
        vals=np.unique(al*V[n][:,0]+be*V[n][:,1])
        res=sorted(set((vals%24).tolist()))
        ok=lambda N,res=res: np.isin(N%24,res)
        a,b,nadm,nbad=full_window(vals,ok)
        print('  order %2d: %8d values in [%d, %d]; residues mod 24 %s; longest fully covered run [%d, %d] (ratio %.2f); admissible misses %d of %d'
              %(n,len(vals),vals[0],vals[-1],res,a,b,b/max(a,1),nbad,nadm))
        allv.append(vals)
    U=np.unique(np.concatenate(allv+[np.unique(al*V[n][:,0]+be*V[n][:,1]) for n in range(1,5)]))
    odd=lambda N: (N%2==1)
    a,b,nadm,nbad=full_window(U,odd)
    print('  union orders <= 10: longest run of odd N all covered [%d, %d]; odd misses in [1,%d]: %d'%(a,b,U[-1],nbad))
