import numpy as np, itertools
R='/Users/e/Documents/GitHub/math/04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt'
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=np.array([nxt() for _ in range(dim)],dtype=object); fin=np.array([nxt() for _ in range(dim)],dtype=object)
MI={}
for _ in range(nl):
    li=nxt(); MI[li]=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=object)
def run(masks_seq):
    y=init.copy()
    for q,m in enumerate(masks_seq):
        mm=sum(1<<(s-1) for s in (1,2,3) if (m>>(s-1))&1 and s<=q)  # spans must be <= q
        y=MI[mm].dot(y)
    return int(fin.dot(y))
pairs=[((0,0,4,6),(0,1,4,6)),((0,2,1,5),(0,2,3,5)),((0,4,2,2),(0,4,2,6)),((1,1,2,3),(1,1,3,3)),
       ((0,6,2,4),(0,6,2,5)),((1,1,4,2),(1,1,6,2)),((0,2,4,4),(0,2,5,4)),((2,4,5,5),(2,5,6,3))]
pre=(0,0,2)       # a 3-cycle start
sufs=[(),(1,2),(0,1,3)]
Ld=5
for u,v in pairs:
    best=None
    for suf in sufs:
        H={}
        for d in itertools.product((0,1),repeat=Ld):
            w=list(pre)
            for b in d: w+=list(v if b else u)
            w+=list(suf); H[d]=run(w)
        z=tuple([0]*Ld)
        wts=[H[tuple(1 if j==i else 0 for j in range(Ld))]-H[z] for i in range(Ld)]
        add=all(H[d]==H[z]+sum(wts[i] for i in range(Ld) if d[i]) for d in H)
        ratios=[wts[i]/wts[i+1] if wts[i+1] else None for i in range(Ld-1)]
        if best is None: best=(suf,add,wts,ratios,H[z])
        if add: best=(suf,add,wts,ratios,H[z]); break
    print('pair',u,v,'additive in digits:',best[1],' weights',best[2],' ratios',[round(r,4) if r else r for r in best[3]])
