#!/usr/bin/env python3
"""Finite positive-representation orbit lattices, exact HNF controls."""
import sys,pathlib,json,math
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'wildcard'))
from check_automaton import compile_automaton

def build(k,arcs):
 v,acc,tr=compile_automaton(k,arcs);d=len(v)
 matrices=[]
 for rows in tr:
  m=S.zeros(d)
  for i,js in enumerate(rows):
   for j in js:m[j,i]+=1
  matrices.append(S.diag(m,1))
 B=S.Matrix(v+[1]);history=[]
 while True:
  new=hermite_normal_form(S.Matrix.hstack(B,*[m*B for m in matrices]))
  history.append(new.cols)
  if new==B:break
  B=new
 return B,matrices,S.Matrix([[int(a) for a in acc]+[0]]),history

def rows_close(B,M):
 while True:
  C=S.Matrix.vstack(B,*[B*m for m in M]);_,p=C.T.rref();N=C[list(p),:]
  if N.rows==B.rows:return N
  B=N

def main():
 out=[]
 for k,arcs in [(1,[[False]]),(2,[[False,True],[False,False]])]:
  B,M,L,hist=build(k,arcs);constant=S.zeros(1,B.rows);constant[0,-1]=1
  parity=[int(x) for x in (L-constant)*B]
  row={'k':k,'lattice_rank':B.cols,'rank_history':hist,'gcd_H_minus_constant':math.gcd(*parity),'pairs':[]}
  if k==2:
   for p in range(4):
    for q in range(4):
     if p==q:continue
     R=rows_close(S.Matrix.vstack(L*M[p]-2*L-constant,L*M[q]-2*L+constant),[M[p],M[q]])
     A=S.Matrix.vstack(R*B,constant*B);target=S.zeros(A.rows,1);target[-1]=1
     H=hermite_normal_form(A);J=hermite_normal_form(A.row_join(target))
     feasible=(H==J)
     # Does every rational solution of the cylinder constraints with const1 force H1?
     aug=S.Matrix.vstack(R*B,constant*B);y=S.zeros(aug.rows,1);y[-1]=1
     try:
      sol,params=aug.gauss_jordan_solve(y);hs=S.expand((L*B*sol)[0]); hs=str(hs)
     except ValueError:hs='inconsistent'
     coords=None
     try:
      zz,pp=H.gauss_jordan_solve(target);coords=[str(x) for x in zz]
     except ValueError:pass
     row['pairs'].append({'plus':p,'minus':q,'integer_affine_feasible':feasible,'rational_H':hs,'constraint_HNF_target_coordinates':coords})
  out.append(row)
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
