import numpy as np, math, random
random.seed(99)
C=105; n1=(C-1)//2
Hf=np.arange(1,n1+1); S=np.sin(2*math.pi*np.outer(Hf,Hf)/C)
def hm(d):
    g=C//d; Kd=set((g*j)%C for j in range(d)); m=np.zeros(n1,bool)
    for x in Kd:
        if 0<x<=n1: m[x-1]=True
    return m
def silent(eps,D):
    Phi=S@eps; B=S[:,D]@eps[D]; A=Phi-B
    return D.any() and np.all(np.abs(A*B)<1e-6)
def search(D,restarts=200,steps=4000):
    best=None
    for r in range(restarts):
        e=np.array([random.choice([-1,1]) for _ in range(n1)]); e[0]=1
        def viol(e_):
            Phi=S@e_; B=S[:,D]@e_[D]; A=Phi-B; return float(np.sum((A*B)**2))
        cur=viol(e); imp=True; st=0
        while imp and st<steps:
            imp=False; st+=1
            order=list(range(1,n1)); random.shuffle(order)
            for k in order:
                e[k]=-e[k]; v=viol(e)
                if v<cur-1e-9: cur=v; imp=True
                else: e[k]=-e[k]
            if cur<1e-7: break
        if best is None or cur<best[0]: best=(cur,e.copy())
        if cur<1e-7: break
    return best
def is_chain(s):
    s=sorted(s); return all(s[i+1]%s[i]==0 for i in range(len(s)-1))
# Coprime-prime combos (predict A=0) vs nested-chain combos (predict A>0)
tests=[((3,5),"coprime"),((3,7),"coprime"),((5,7),"coprime"),
       ((3,15),"chain 3|15"),((5,15),"chain 5|15"),((7,35),"chain 7|35"),
       ((3,21),"chain 3|21"),((5,35),"chain 5|35"),((7,21),"chain 7|21"),
       ((15,21),"non-chain (gcd3)"),((15,35),"non-chain (gcd5)")]
print(f"C=105: H_d (+) H_e silent-ability (search; A_t*B_t lemma verifies):")
for (d,e),tag in tests:
    D=hm(d)^hm(e)
    b=search(D)
    found=silent(b[1],D)
    print(f"  H{d}+H{e:>3} [{tag:>17}] chain={is_chain((d,e))!s:>5}  minviol={b[0]:.2e}  SILENT-CUT={found}")
