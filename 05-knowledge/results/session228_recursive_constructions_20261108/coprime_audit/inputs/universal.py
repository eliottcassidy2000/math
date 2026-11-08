import sys
sys.path.append('/workspace/scratch/07d5532f2fce/session228-fiber-switch')
from transfer import *
import sympy as s
x,y,a,b=s.symbols('x y a b');weights={(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b}
F=s.zeros(6)
for j,(p,r) in enumerate(ST):
 for (pp,rr),c in weights.items():
  if not r and rr==2:continue
  F[ST.index((pp or p,int(bool(r or rr)))),j]+=c
S0,S1=map(s.Matrix,S)
co=[F.diff(v) for v in [x,y,a,b]]
K=s.Matrix.hstack(*co);K=s.Matrix.hstack(*K.columnspace())
L=s.Matrix.hstack(S0*K,S1*K);L=s.Matrix.hstack(*L.columnspace())
print('F',F,'Krank',K.rank(),'L',L,'FL',s.simplify(F*L),'invariant',L.row_join(S0*L).row_join(S1*L).rank()==L.rank())
def require(x,msg):
 if not x:raise RuntimeError(msg)
W=s.Matrix.hstack(s.eye(6)[:,1]-s.eye(6)[:,5],s.eye(6)[:,3]-s.eye(6)[:,5])
require(F*W==s.zeros(6,2),'F annihilation')
require(W.row_join(S0*W).row_join(S1*W).rank()==2,'W invariance')
for C in co:
 require(W.row_join(S0*C).row_join(S1*C).rank()==2,'S_b image F')
print('PASS: universal symbolic separated-defect law')
