import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, sys, itertools
sys.path.insert(0,RIB)
from sympy import Matrix, symbols, factor
w=int(sys.argv[1]); units=sys.argv[2].split(',')
R=os.path.join(RIB,f'span{w}_intrep.txt')
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=Matrix([nxt() for _ in range(dim)]); fin=Matrix([[nxt() for _ in range(dim)]])
MATS={}
for _ in range(nl):
    li=nxt(); MATS[li]=Matrix([[nxt() for _ in range(dim)] for _ in range(dim)])
def mask(t): return sum(1<<(int(ch)-1) for ch in t if ch in '123')
def Mw(b):
    M=Matrix.eye(dim)
    for t in b.split('.'): M=MATS[mask(t)]*M
    return M
U=[Mw(u) for u in units]
pre=Mw('0')*init
basis=pre; frontier=[pre]
while frontier:
    nf=[]
    for v in frontier:
        for M in U:
            wv=M*v
            if Matrix.hstack(basis,wv).rank()>basis.rank(): basis=Matrix.hstack(basis,wv); nf.append(wv)
    frontier=nf
obs=fin; frontier=[fin]
while frontier:
    nf=[]
    for u in frontier:
        for M in U:
            wv=u*M
            if Matrix.vstack(obs,wv).rank()>obs.rank(): obs=Matrix.vstack(obs,wv); nf.append(wv)
    frontier=nf
Hk=obs*basis; r=Hk.rank()
print('units',units,'reach',basis.shape[1],'obs',obs.shape[0],'minimal dim',r)
rows=[];cur=Matrix.zeros(0,Hk.shape[1])
for i in range(Hk.shape[0]):
    t=Matrix.vstack(cur,Hk[i,:])
    if t.rank()>cur.rank(): cur=t; rows.append(i)
    if len(rows)==r: break
cols=[];cur=Matrix.zeros(Hk.shape[0],0)
for j in range(Hk.shape[1]):
    t=Matrix.hstack(cur,Hk[:,j])
    if t.rank()>cur.rank(): cur=t; cols.append(j)
    if len(cols)==r: break
O=obs.extract(rows,list(range(obs.shape[1]))); Bs=basis.extract(list(range(basis.shape[0])),cols)
Hi=(O*Bs).inv()
red=[O*M*Bs*Hi for M in U]
x=symbols('x')
for u,M in zip(units,red):
    ranks=[]; P=Matrix.eye(r)
    for j in range(1,5): P=M*P; ranks.append(P.rank())
    print(' unit',u,'charpoly',factor(M.charpoly(x).as_expr()),'power ranks',ranks)
# pairwise products ranks
print(' pair ranks:')
for (i,u),(j,v) in itertools.product(enumerate(units),repeat=2):
    print('  ',u,'then',v,': rank',(red[j]*red[i]).rank(), 'charpoly', factor((red[j]*red[i]).charpoly(x).as_expr()))
import pickle
pickle.dump({'units':units,'red':red,'z0':O*pre,'fin':fin*Bs*Hi},open(f'unitrep_w{w}.pkl','wb'))
