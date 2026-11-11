from check_automaton import *
import numpy as np
A=[[False,True],[False,False]]
init,acc,tr=compile_automaton(2,A);D=len(init);phi=[3,1,2,0]
M=[]
for rows in tr:
 a=np.zeros((D,D),dtype=object)
 for i,js in enumerate(rows):
  for j in js:a[j,i]+=1
 M.append(a)
def star(w):return tuple(phi[x] for x in w[::-1])
def mat(w):
 x=np.eye(D,dtype=object)
 for b in w:x=M[b]@x
 return x
rng=random.Random(261109);checks=0
for _ in range(100):
 w=tuple(rng.randrange(4) for _ in range(rng.randrange(9)))
 X=mat(w);Y=mat(star(w))
 assert star(star(w))==w
 assert automaton_count(2,A,w)==automaton_count(2,A,star(w))
 for b in range(4):
  assert np.array_equal(mat(w+(b,)),M[b]@X)
  assert np.array_equal(mat(star(w+(b,))),Y@M[phi[b]])
  assert np.array_equal(mat((b,)+w),X@M[b])
  assert np.array_equal(mat(star((b,)+w)),M[phi[b]]@Y)
  checks+=1
print({'seed':261109,'word_pairs':100,'two_sided_matrix_updates':checks,'exact_integer_dtype':'object','passed':True})
