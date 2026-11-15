import pickle, numpy as np
import os as _o; tok=open(_o.environ.get('SPAN3_INTREP','../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt')).read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); [nx() for _ in range(2*dim)]
M={}
for _ in range(nl):
    li=nx(); M[li]=np.array([[nx() for _ in range(dim)] for _ in range(dim)],dtype=float)
code=[(i,a) for i in range(3) for a in range(1,8)]
need=pickle.load(open('memneed_1.765.pkl','rb'))['need']
def growth(c):
    P=np.identity(dim)
    for u in c:
        i,a=code[u]
        for l in (0,)*i+(a,): P=M[l]@P
    return np.abs(np.linalg.eigvals(P)).max()**(1/sum(code[u][0]+1 for u in c))
g={c:growth(c) for c in need}
name=lambda c:'|'.join('0'*code[u][0]+str(code[u][1]) for u in c)
for m in (1,2,3,4,5,6):
    unc=[(g[c],c) for c,L in need.items() if L>m+1]
    if not unc: print('memory %d: every cheap cycle cut'%m); continue
    gmin,cmin=min(unc)
    print('memory %d code units (edge check on %d-words): cheapest uncut cycle growth %.5f  %s  (bottleneck length %d); uncut cycles %d'
          %(m,m+1,gmin,name(cmin),need[cmin],len(unc)))
