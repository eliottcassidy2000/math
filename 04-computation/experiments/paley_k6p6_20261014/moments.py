import json, numpy as np, sympy
W=json.load(open('k6p6_W_allprimes.json')); W={int(k):v for k,v in W.items()}
def chi_table(q):
    c=np.full(q,-1,dtype=np.int64); c[np.unique((np.arange(1,q)**2)%q)]=1; c[0]=0; return c
rows=[]
for q in sympy.primerange(5,1001):
    c=chi_table(q); x=np.arange(q)
    base=c[x][None,:]*c[(x-1)%q][None,:]*c[(x[None,:]-x[:,None])%q]      # base[a,y]=chi(y)chi(y-1)chi(y-a)
    t=base.sum(axis=1)                                                     # t(a) Legendre trace (= -a_q(E_a))
    K=c[x][:,None]*c[x][None,:]*c[(x-1)%q][None,:]*c[(x[None,:]-x[:,None])%q]   # K[a,x]=chi(a)chi(x)chi(x-1)chi(x-a)
    M=-int(t@K@t)
    ca=c[x]; cam1=c[(x-1)%q]
    feats=dict(M=M, m2=int((t*t).sum()), m2a=int((ca*t*t).sum()), m2b=int((cam1*t*t).sum()), m2ab=int((ca*cam1*t*t).sum()),
               m3=int((t**3).sum()), m3a=int((ca*t**3).sum()), m4=int((t**4).sum()))
    rows.append((q,W[q],feats))
json.dump(rows,open('moments.json','w'))
qs=np.array([r[0] for r in rows],dtype=float); Wv=np.array([r[1] for r in rows],dtype=float)
names=list(rows[0][2].keys())
F={n:np.array([r[2][n] for r in rows],dtype=float) for n in names}
poly=[qs**3,qs**2,qs,np.ones_like(qs)]
def fit(cols,label):
    A=np.column_stack(cols); c,_,rk,_=np.linalg.lstsq(A,Wv,rcond=None); err=np.abs(Wv-A@c).max()
    print(f'{label:40s} rank {rk:2d}  max|res| {err:.4g}  coeffs {np.round(c,5)}')
fit([F['M']]+poly,'W ~ M + poly')
for r in (1,3):
    sel=(qs%4==r)
    A=np.column_stack([F['M'][sel]]+[p[sel] for p in poly]); c,_,_,_=np.linalg.lstsq(A,Wv[sel],rcond=None)
    print(f'  q={r} mod 4 only: max|res| {np.abs(Wv[sel]-A@c).max():.4g}, coeffs {np.round(c,5)}')
fit([F[n] for n in names]+poly,'W ~ all moments + poly')
print('sample q, W, M:',[(r[0],r[1],r[2]['M']) for r in rows[:10]])
