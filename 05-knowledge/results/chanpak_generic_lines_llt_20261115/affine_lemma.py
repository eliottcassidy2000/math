# Lemma: for a signed difference sigma in {-1,0,1}^a (x - y bitwise) and free bits z on Z = {k: sigma_k = 0},
#   p(x) - p(y) = c_0(sigma) + sum_{k in Z} c_k(sigma) z_k,
#   c_k = - sum_{j>k} sigma_j 2^(j-k-1) - sum_{i<k} sigma_i 2^(k-i-1) - sigma_(k+1) - sigma_(k-1).
import itertools, random
def p(u):
    r=0
    while u>1:
        w=u>>1
        if not (u&1): r+=w+(w&1)
        u=w
    return r
def coeffs(sig):
    a=len(sig); c={}
    for k in range(a):
        if sig[k]!=0: continue
        s=-sum(sig[j]*2**(j-k-1) for j in range(k+1,a))-sum(sig[i]*2**(k-i-1) for i in range(k))
        s-= (sig[k+1] if k+1<a else 0)+(sig[k-1] if k>0 else 0)
        c[k]=s
    return c
bad=0; tested=0
for a in range(1,11):
    sigs=list(itertools.product((-1,0,1),repeat=a))
    for sig in (sigs if a<=7 else random.Random(a).sample(sigs,400)):
        Z=[k for k in range(a) if sig[k]==0]; c=coeffs(sig)
        base=None
        for zb in itertools.product((0,1),repeat=len(Z)) if len(Z)<=7 else [tuple(random.Random(sum(sig)+t).randint(0,1) for _ in Z) for t in range(64)]:
            xb=[0]*a; yb=[0]*a
            for k in range(a):
                if sig[k]==1: xb[k],yb[k]=1,0
                elif sig[k]==-1: xb[k],yb[k]=0,1
            for k,zk in zip(Z,zb): xb[k]=yb[k]=zk
            x=sum(b<<k for k,b in enumerate(xb)); y=sum(b<<k for k,b in enumerate(yb))
            lin=sum(c[k]*zk for k,zk in zip(Z,zb))
            val=p(x)-p(y)-lin
            if base is None: base=val
            elif val!=base: bad+=1
            tested+=1
print('affine-in-free-bits lemma: %d evaluations, %d violations'%(tested,bad))
