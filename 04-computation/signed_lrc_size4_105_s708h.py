import numpy as np, math, random
random.seed(7)
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
# two independent chain moves; try to make BOTH silent at one cut -> size-4 class
pairs=[((3,15),(7,21)),((3,15),(7,35)),((5,15),(7,21)),((3,21),(5,15))]
for (a,b),(c,d) in pairs:
    D1=hm(a)^hm(b); D2=hm(c)^hm(d)
    def viol(e):
        Phi=S@e
        B1=S[:,D1]@e[D1]; A1=Phi-B1
        B2=S[:,D2]@e[D2]; A2=Phi-B2
        return float(np.sum((A1*B1)**2)+np.sum((A2*B2)**2))
    best=None
    for r in range(120):
        e=np.array([random.choice([-1,1]) for _ in range(n1)]); e[0]=1
        cur=viol(e); imp=True; st=0
        while imp and st<3000:
            imp=False; st+=1
            order=list(range(1,n1)); random.shuffle(order)
            for k in order:
                e[k]=-e[k]; v=viol(e)
                if v<cur-1e-9: cur=v; imp=True
                else: e[k]=-e[k]
            if cur<1e-6: break
        if best is None or cur<best[0]: best=(cur,e.copy())
        if cur<1e-6: break
    ok = silent(best[1],D1) and silent(best[1],D2) and (D1^D2).any()
    # also check the product D1^D2 is silent (group closure -> Klein 4)
    prod_sil = silent(best[1], D1^D2)
    print(f"chains H{a}+H{b} & H{c}+H{d}: minviol={best[0]:.2e} both-silent={ok} product-also-silent(=>size4)={prod_sil}")
