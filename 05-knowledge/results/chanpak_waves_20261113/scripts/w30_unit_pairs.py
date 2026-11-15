import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import pickle, itertools
from sympy import symbols, factor, Matrix
from collections import Counter
D=pickle.load(open(os.path.join(HERE,'w30_hankel_13units.pkl'),'rb')); units=D['units']; red=D['red']; x=symbols('x')
print('unit ranks:',{u:M.rank() for u,M in zip(units,red)})
print('unit charpolys:',Counter(str(factor(M.charpoly(x).as_expr())) for M in red))
rk=Counter(); comp=[]; inc=[]
for (i,u),(j,v) in itertools.product(enumerate(units),repeat=2):
    P=red[j]*red[i]; r=P.rank(); rk[r]+=1
    cp=factor(P.charpoly(x).as_expr())
    nz=[ev for ev,m in P.eigenvals().items() if ev!=0]
    if len(nz)==1: comp.append((u,v,r))
    else: inc.append((u,v,r,str(cp)))
print('pair product ranks:',dict(rk))
print('pairs with a single nonzero eigenvalue (compatible):',len(comp),' incompatible:',len(inc))
r1=[(u,v) for u,v,r in comp if r==1]
print('rank-one pairs:',len(r1))
# powers: rank of u^k
for u,M in zip(units,red):
    P=M; rs=[]
    for k in range(1,6): rs.append(P.rank()); P=M*P
    print('  ',u,'power ranks',rs)
