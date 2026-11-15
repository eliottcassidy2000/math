import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import pickle, itertools, sys
from sympy import Matrix, Rational, nsimplify, symbols, factor, sqrt
D=pickle.load(open(sys.argv[1],'rb'))
units=D['units']; red=D['red']; z0=D['z0']; fin=D['fin']
r=red[0].shape[0]
print('units',units,'dim',r)
for u,M in zip(units,red): print(u, M.tolist())
print('z0',z0.T.tolist(),'fin',fin.tolist())
# image plane automaton: state = column space of current product (rank>=2 required), canonical rref of basis
def canon(M):
    C=M.columnspace()
    B=Matrix.hstack(*C)
    R=B.T.rref()[0]
    R=R[:B.shape[1],:]
    return tuple(tuple(x) for x in R.tolist()), B.shape[1]
start=Matrix.eye(r)
states={}; order=[]; trans={}
# begin from images of single units
frontier=[]
for i,M in enumerate(red):
    key,rk=canon(M)
    if key not in states: states[key]=len(order); order.append((key,M,rk)); frontier.append(key)
while frontier:
    nf=[]
    for key in frontier:
        _,M,rk=order[states[key]]
        for j,U in enumerate(red):
            P=U*M
            if P.rank()<2: trans[(states[key],j)]=None; continue
            k2,rk2=canon(P)
            if k2 not in states: states[k2]=len(order); order.append((k2,P,rk2)); nf.append(k2)
            trans[(states[key],j)]=states[k2]
    frontier=nf
    if len(order)>200: print('too many states'); break
print('plane states:',len(order),' ranks',[o[2] for o in order])
for s in range(len(order)):
    print(' state',s,'rank',order[s][2],'->',[ (units[j], trans.get((s,j))) for j in range(len(units))])
