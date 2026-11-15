# count cross pairs (H0 top bit 1, H1 top bit 0, lowest bits differ) of a-bit numbers with
#   h = H0-H1 = rev(H0)-rev(H1) = d   and   p(H0) = p(H1)   (anti-diagonal collision lines)
import numpy as np, sys
amax=int(sys.argv[1])
N=1<<amax
p=np.zeros(N,dtype=np.int64)
for u in range(1,N):
    w=u>>1; p[u]=p[w] if u&1 else p[w]+w+(w&1)
for a in range(2,amax+1):
    M=1<<a
    H=np.arange(M,dtype=np.int64)
    r=np.zeros(M,dtype=np.int64)
    for i in range(a): r|=((H>>i)&1)<<(a-1-i)
    key=H-r                       # h=d  <=>  H0-rev(H0) = H1-rev(H1)
    pa=p[:M]
    top=H>=M//2
    from collections import defaultdict
    A=defaultdict(lambda:[0,0,0,0])   # key,(p) -> counts: [top&odd, top&even, low&odd, low&even]
    for h_,k_,pp,t,o in zip(H,key,pa,top,H&1):
        A[(int(k_),int(pp))][(0 if t else 2)+(0 if o else 1)]+=1
    hd_p=sum(v[0]*v[3]+v[1]*v[2] for v in A.values())
    # h=d only (any p)
    B=defaultdict(lambda:[0,0,0,0])
    for h_,k_,t,o in zip(H,key,top,H&1):
        B[int(k_)][(0 if t else 2)+(0 if o else 1)]+=1
    hd=sum(v[0]*v[3]+v[1]*v[2] for v in B.values())
    print('a=%2d  cross pairs with h=d: %8d (%.3f * 2^a)   with also p(H0)=p(H1): %6d (%.3f * 2^(a/2))'%(a,hd,hd/2**a,hd_p,hd_p/2**(a/2)))
