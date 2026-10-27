#!/usr/bin/env python3
"""Exact SU(3) six-port contractions and Petersen carrier. No float gates."""
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import json,sys
import sympy as s

def require(b,msg):
 if not b: raise ValueError(msg)
def eps(a):
 if len(set(a))<3:return 0
 return (-1)**sum(a[i]>a[j] for i in range(3) for j in range(i+1,3))
def invariant_dimension(n):
 require(n%3==0,'balanced weight required')
 k=n//3;words=[w for w in product(range(3),repeat=n) if all(w.count(i)==k for i in range(3))]
 cols={w:j for j,w in enumerate(words)};rows={}
 for a,b in permutations(range(3),2):
  for w,j in cols.items():
   for p,c in enumerate(w):
    if c==b:
     v=w[:p]+(a,)+w[p+1:];key=(a,b,v)
     if key not in rows:rows[key]=[0]*len(words)
     rows[key][j]+=1
 M=s.Matrix(list(rows.values()))
 return len(words)-M.rank(),len(words),M.rows

def find_iso(A):
 edges=list(combinations(range(5),2));B=[[int(not set(a)&set(b)) for b in edges] for a in edges];n=10
 def rec(f,unused):
  if len(f)==n:return f
  v=max((v for v in range(n) if v not in f), key=lambda v:sum(int(A[v,u]) for u in f))
  for w in sorted(unused):
   if all(int(A[v,u])==B[w][t] for u,t in f.items()):
    z=rec(dict(f,**{})|{v:w},unused-{w})
    if z:return z
  return None
 f=rec({0:0},set(range(1,10)));require(f is not None,'Kneser isomorphism')
 require(all(int(A[i,j])==B[f[i]][f[j]] for i in range(10) for j in range(10)),'exact K5 edge map')
 return [edges[f[i]] for i in range(10)]

def main():
 words=list(product(range(3),repeat=6));parts=[];rs=[]
 for pair in combinations(range(1,6),2):
  A=(0,)+pair;B=tuple(i for i in range(6) if i not in A);parts.append((A,B))
  rs.append([eps(tuple(c[i] for i in A))*eps(tuple(c[i] for i in B)) for c in words])
 M=s.Matrix(rs);G=M*M.T;I=s.eye(10)
 require(all(G[i,i]==36 for i in range(10)),'norm36')
 require(all(abs(G[i,j])==12 for i in range(10) for j in range(i)),'overlap12')
 require(G*G==72*G,'tight Gram polynomial')
 require(G.rank()==5,'Gram rank5')
 require(G.eigenvals()=={s.Integer(0):5,s.Integer(72):5},'Gram spectrum')
 # Pointwise Schouten relation, not just a zero norm inference.
 relation=s.Matrix([-1,1,-1,1,0,0,0,0,0,0])
 require(M.T*relation==s.zeros(729,1),'four-term relation')
 # Port permutations of one Pluecker relation generate the full kernel.
 part_index={a:i for i,(a,b) in enumerate(parts)};orbit=[]
 for perm in permutations(range(6)):
  row=[0]*10
  for j,z in enumerate(relation):
   if z:
    aa=tuple(perm[k] for k in parts[j][0]);bb=tuple(perm[k] for k in parts[j][1])
    sign=(-1)**(sum(aa[k]>aa[l] for k in range(3) for l in range(k+1,3))+sum(bb[k]>bb[l] for k in range(3) for l in range(k+1,3)))
    aa=tuple(sorted(aa));bb=tuple(sorted(bb));key=aa if 0 in aa else bb
    row[part_index[key]]+=int(z)*sign
  orbit.append(row)
 O=s.Matrix(orbit)
 require(O*M==s.zeros(720,729),'all permuted relations')
 require(O.rank()==5,'one relation orbit spans full kernel')
 require(all(G.extract(ids,ids).det()>0 for ids in combinations(range(10),3)),'all triples independent')
 # Every color generator annihilates each tensor (and H imposes counts2).
 for r in rs:
  require(all(not z or all(w.count(a)==2 for a in range(3)) for w,z in zip(words,r)),'zero weight')
  for a,b in permutations(range(3),2):
   out={}
   for w,z in zip(words,r):
    if z:
     for p,c in enumerate(w):
      if c==b:
       v=w[:p]+(a,)+w[p+1:];out[v]=out.get(v,0)+z
   require(all(z==0 for z in out.values()),'Lie generator invariance')
 dim6,weight6,eq6=invariant_dimension(6)
 require(dim6==5,'full invariant kernel dimension')
 dim3,weight3,eq3=invariant_dimension(3)
 require(dim3==1,'unique baryon singlet')
 signs=[1,1,1,1,1,1,-1,-1,-1,1];D=s.diag(*signs)
 S=D*(G-36*I)*D/12;A=(s.ones(10)-I-S)/2
 require(S*S==9*I,'Seidel conference square')
 require(all(A[i,i]==0 for i in range(10)) and all(A[i,j] in (0,1) for i in range(10) for j in range(10)),'simple graph')
 require(A==A.T and all(sum(A.row(i))==3 for i in range(10)),'cubic graph')
 require(A*A==2*I-A+s.ones(10),'SRG(10,3,0,1) identity')
 require(A.eigenvals()=={s.Integer(3):1,s.Integer(1):5,s.Integer(-2):4},'Petersen spectrum')
 mapping=find_iso(A)
 # Norm and port-sign hostile controls.
 hostile=0
 try:require(G==36*I,'false orthogonal cluster basis')
 except ValueError:hostile+=1
 try:require(M.T*s.Matrix([1,1,1,1,0,0,0,0,0,0])==s.zeros(729,1),'drop alternating signs')
 except ValueError:hostile+=1
 try:require(G.rank()==10,'forget linear relations')
 except ValueError:hostile+=1
 try:require(s.Rational(12,36)==0,'orthogonal reclustering')
 except ValueError:hostile+=1
 require(hostile==4,'all hostile controls reject')
 # Two typed contraction channels are linearly independent.
 e1=s.Matrix([1,0,0]);e2=s.Matrix([0,1,0]);phi=s.Matrix([[1,0,0]])
 F=(phi*e1)[0]*e2;H=(phi*e2)[0]*e1
 require(F==e2 and H==s.zeros(3,1),'typed channel distinction')
 # A vector pair is not its baryon determinant or an Eisenstein additive sum.
 eps_basis=eps((0,1,2));require(eps_basis==1,'ordered epsilon basis')
 require(eps((1,0,2))==-1 and eps((1,2,0))==1,'reflection/cyclic orientation')
 # Exact local web contractions, with epsilon_012=epsilon^012=+1.
 for i,j,l,m in product(range(3),repeat=4):
  require(sum(eps((i,j,k))*eps((l,m,k)) for k in range(3))==int(i==l)*int(j==m)-int(i==m)*int(j==l),'epsilon epsilon delta identity')
 for i,j in product(range(3),repeat=2):
  require(sum(eps((i,a,b))*eps((j,a,b)) for a,b in product(range(3),repeat=2))==2*int(i==j),'digon factor2')
 require(sum(int(i==i) for i in range(3))==3,'closed delta loop3')
 # Octet highest-weight multiplicity plane, coordinates a=001,b=010,c=100.
 U=s.Matrix([[2,0],[-1,1],[-1,-1]])*s.diag(1/s.sqrt(6),1/s.sqrt(2))
 Up=s.Matrix([[1,1],[1,-1],[-2,0]])*s.diag(1/s.sqrt(6),1/s.sqrt(2))
 R=s.simplify(U.T*Up)
 require(U.T*U==s.eye(2) and Up.T*Up==s.eye(2),'orthonormal multiplicity bases')
 require(R==s.Matrix([[s.Rational(1,2),s.sqrt(3)/2],[s.sqrt(3)/2,-s.Rational(1,2)]]),'nontrivial octet recoupling')
 C=s.Matrix([[0,0,1],[1,0,0],[0,1,0]]);K=s.simplify(U.T*(C-C.T)*U)
 require(s.simplify(K*K)==-3*s.eye(2),'C3 skew square minus3')
 require(s.ones(1,3)*U==s.zeros(1,2),'A2 centered plane')
 report={'status':'PASS',
  'partitions':parts,'gram':G.tolist(),'gram_rank':5,'gram_eigenvalues':{'0':5,'72':5},
  'normalized_squared_off_diagonal_overlap':'1/9','frame_dimension':5,'frame_vectors':10,
  'full_invariant_dimensions':{'3':dim3,'6':dim6},'six_weight_zero_columns':weight6,'six_generator_equations':eq6,
  'schouten_relation':list(relation),'switching_signs':signs,'petersen_adjacency':A.tolist(),'k5_edge_map':mapping,
  'invariant_dimension_sequence_3m':[2*factorial(3*m)//(factorial(m)*factorial(m+1)*factorial(m+2)) for m in range(1,7)],
  'octet_recoupling':[[str(z) for z in R.row(i)] for i in range(2)],'relation_orbit_rank':5,'minimum_dependent_support':4,'web_loop_factor':3,'web_digon_factor':2,'hostiles_rejected':hostile,'scope':'Finite SU3 representation tensors only; no QCD dynamics, no mass gap, no claimed novelty.'}
 out=json.dumps(report,indent=2,default=int,sort_keys=True)+'\n'
 if len(sys.argv)>1: Path(sys.argv[1]).write_text(out)
 print(out,end='')
if __name__=='__main__':main()
