#!/usr/bin/env python3
"""Independent integer Smith check of the palette audit appendix."""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
import json

def need(x,m):
 if not x:raise RuntimeError(m)
records=[]
for q in range(2,9):
 A=sp.ones(q)-sp.eye(q)
 rows=[]
 for i in range(q):
  for j in range(q):
   z=[0]*(q*q)
   for k in range(q):z[i*q+k]=int(A[j,k])
   rows.append(z)
 B=sp.Matrix(rows)
 D=smith_normal_form(B,domain=sp.ZZ)
 pos=sorted(abs(int(D[i,i])) for i in range(q*q))
 need(pos==sorted([1]*(q*(q-1))+[q-1]*q),'positive Smith invariants')
 empty=[int(i==j) for i in range(q) for j in range(q)]
 D2=smith_normal_form(B.col_join(sp.Matrix([empty])),domain=sp.ZZ)
 whole=sorted(abs(int(D2[i,i])) for i in range(q*q))
 need(whole==sorted([1]*(q*q-q+1)+[q-1]*(q-1)),'identity Smith invariants')
 records.append({'q':q,'positive_invariants':pos,'including_identity_invariants':whole,'index_positive':int(abs(B.det())),'index_including_identity':(q-1)**(q-1)})
print(json.dumps({'status':'PASS','universe':'q2..8 exact integer matrices of singleton-generated positive words, with and without identity','records':records},indent=2))
