from check_automaton import *
import sympy as S
from math import factorial
for k in range(4):
 A=[[i<j for j in range(k)] for i in range(k)]
 init,acc,tr=compile_automaton(k,A)
 d=2*3**k
 # Single permutation block. Original local order(empty,good,bad).
 C=S.Matrix([[0,0,1],[1,1,0],[-1,0,0]])
 E=S.Matrix([[0,1],[1,0]])
 P=E
 for _ in range(k):P=S.kronecker_product(C,P)
 Pi=P.inv();ms=[]
 for c,rows in enumerate(tr):
  M=S.zeros(d)
  for i,js in enumerate(rows[:d]):
   for j in js:M[j,i]+=1
  B=Pi*M*P
  assert all(B[i,j]==0 for i in range(d) for j in range(i))
  if ms: assert list(B.diagonal())==list(ms[0].diagonal())
  ms.append(B)
 for word in [(0,),((1<<k)-1,0)]:
  B=S.eye(d)
  for c in word:B=ms[c]*B
  L=len(word)
  assert list(B.diagonal())==[x**L for x in ms[0].diagonal()]
  p=S.eye(d)
  for j in range(k+2):p=p*(B-j**L*S.eye(d))
  assert p**(2*k+2)==S.zeros(d)
 print({'k':k,'per_permutation_dimension':d,'alphabet_size':1<<k,'common_triangularization':True,'uniform_annihilator':True})
