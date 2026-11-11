from check_automaton import *
import sympy as S
from itertools import combinations
k=2;a=[[False,True],[False,False]]
v,l,tr=compile_automaton(k,a);D=len(v)
M=[]
for rows in tr:
 m=S.zeros(D)
 for i,js in enumerate(rows):
  for j in js:m[j,i]+=1
 M.append(m)
L=S.Matrix([int(x) for x in l]).T
for p,q in combinations(range(4),2):
 eq=(M[p]-2*S.eye(D)).col_join(M[q]-2*S.eye(D))
 basis=eq.nullspace()
 seen=[(L*b)[0] for b in basis]
 print(p,q,'common_dim',len(basis),'count_visible',any(seen))
for p,q in __import__('itertools').permutations(range(4),2):
 I=S.eye(D);Z=S.zeros(D)
 E=(M[p]-2*I).row_join(Z).col_join((M[q]-2*I).row_join(Z)).col_join((-I).row_join(M[p]-I)).col_join(I.row_join(M[q]-I)).col_join(L.row_join(S.zeros(1,D))).col_join(S.zeros(1,D).row_join(L))
 y=S.zeros(4*D+2,1);y[4*D]=1
 sol=S.linsolve((E,y))
 print('digits',p,q,'consistent',sol!=S.EmptySet)
