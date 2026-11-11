from check_automaton import *
from itertools import permutations
import sympy as S
k=2;a=[[False,True],[False,False]]
v,acc,tr=compile_automaton(k,a);D=len(v);L=S.Matrix([int(x) for x in acc]).T
M=[]
for rows in tr:
 m=S.zeros(D)
 for i,js in enumerate(rows):
  for j in js:m[j,i]+=1
 M.append(S.diag(m,1))
certs={}
for p,q in permutations(range(4),2):
 B=(L*M[p][:D,:D]-2*L).row_join(S.Matrix([[-1]])).col_join((L*M[q][:D,:D]-2*L).row_join(S.Matrix([[1]])))
 while True:
  old=B.rows
  C=B.col_join(B*M[p]).col_join(B*M[q]);rr,piv=C.T.rref();B=C[list(piv),:]
  if B.rows==old:break
 certs[p,q]=[[int(x) for x in r] for r in B.tolist()]
 print('pair',p,q,'rank',B.rows,flush=True)
wave=[(v,())];found=[]
for n in range(9):
 for v,w in wave:
  vv=v+[1]
  for pair,rows in certs.items():
   if all(sum(x*y for x,y in zip(r,vv))==0 for r in rows):found.append((pair,w,sum(x for x,a in zip(v,acc) if a)))
 print('length',n,'words',len(wave),'found',found[:5],flush=True)
 if n<8:wave=[(step(v,tr[c]),w+(c,)) for v,w in wave for c in range(4)]
