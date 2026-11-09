import sys,json,itertools,os
from pathlib import Path
import numpy as np
import sympy as S
p=Path(os.environ.get('MATH_REPO',str(Path(__file__).resolve().parents[4])))/'04-computation/experiments/ribbon_spectrum_20261104'
sys.path.insert(0,str(p))
from intrep import build
from ribbon import letters
basis,A,y0,h,*_=build(2,False)
L=letters(2); mats=[np.array(A[a],dtype=np.int64) for a in L];n=len(y0)
def singular(M,p=1000003):
 R=[[int(v)%p for v in row] for row in M]
 for c in range(n):
  q=next((r for r in range(c,n) if R[r][c]),None)
  if q is None:return True
  R[c],R[q]=R[q],R[c];inv=pow(R[c][c],p-2,p)
  for r in range(c+1,n):
   if R[r][c]:
    v=R[r][c]*inv%p
    for j in range(c+1,n):R[r][j]=(R[r][j]-v*R[c][j])%p
 return False
found=[];tested=0
front=[((),np.eye(n,dtype=np.int64))]
for length in range(1,7):
 front=[(w+(i,),mats[i]@M) for w,M in front for i in range(4)]
 for word,M in front:
  tested+=1
  for a in (2,4,8):
   B=M-a*np.eye(n,dtype=np.int64)
   if singular(B) and S.Matrix(B).det()==0:
    found.append({'word':word,'radix':a,'eigen1':S.Matrix(M-np.eye(n,dtype=np.int64)).det()==0})
 print(json.dumps({'length':length,'tested':tested,'found':found}),flush=True)
