#!/usr/bin/env python3
"""Exact finite closure certificate for all-word Hankel rank, not sampled rank."""
import sympy as s,json
from spine import machine

def basis_closure(start,matrices):
 B=[start];words=[()];j=0
 while j<len(B):
  x=B[j];w=words[j];j+=1
  for c,M in enumerate(matrices):
   y=M*x
   if s.Matrix.hstack(*B,y).rank()>len(B):B.append(y);words.append(w+(c,))
 return s.Matrix.hstack(*B),words

def certificate():
 init,a,E=machine((2,0));D=len(init);M=[]
 for es in E:
  A=s.zeros(D)
  for i,j in es:A[j,i]+=1
  M.append(A)
 B,p=basis_closure(s.Matrix(init),M);C,q=basis_closure(s.Matrix(a),[A.T for A in M]);G=C.T*B
 cs=list(G.rref()[1]);rs=list(G[:,cs].T.rref()[1]);K=G.extract(rs,cs)
 # Closure tested explicitly rather than inferred from observed stabilization.
 for A in M:
  if B.row_join(A*B).rank()!=B.cols or C.row_join(A.T*C).rank()!=C.cols:raise RuntimeError('closure')
 f=s.Matrix([a])*B[:,cs];start=K.inv()*C[:,rs].T*s.Matrix(init)
 reduced=[K.inv()*C[:,rs].T*A*B[:,cs] for A in M]
 # Equality in the full observable quotient, for all states in the reachable span.
 for A,R in zip(M,reduced):
  if C.T*(A*B[:,cs]-B[:,cs]*R)!=s.zeros(C.cols,len(cs)):raise RuntimeError('quotient')
 if C.T*(s.Matrix(init)-B[:,cs]*start)!=s.zeros(C.cols,1):raise RuntimeError('initial')
 return {'nonnegative_dimension':D,'reachable_dimension':B.cols,'observable_dimension':C.cols,
 'minimal_rational_dimension':len(cs),'minor_determinant':str(K.det()),
 'prefix_words':[p[j] for j in cs],'suffix_words':[list(reversed(q[i])) for i in rs],
 'initial':[str(x) for x in start],'output':[str(x) for x in f],
 'matrices':[[[str(R[i,j]) for j in range(R.cols)] for i in range(R.rows)] for R in reduced]}
if __name__=='__main__':print(json.dumps(certificate(),indent=2))
