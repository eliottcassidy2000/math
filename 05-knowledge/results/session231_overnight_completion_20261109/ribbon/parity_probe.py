import os,sys,json,random
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(os.environ.get('MATH_REPO',str(Path(__file__).resolve().parents[4])))/'04-computation/experiments/ribbon_spectrum_20261104'))
from intrep import build
_,A,_,_,*_=build(2,False);mats=[np.array(M,dtype=np.int64) for M in A.values()];n=14

def rank(M,p):
 R=M.copy()%p;r=0
 for c in range(n):
  nz=np.flatnonzero(R[r:,c])
  if not len(nz):continue
  q=r+int(nz[0]);R[[r,q]]=R[[q,r]]
  coeff=R[r+1:,c]*pow(int(R[r,c]),p-2,p)%p
  R[r+1:,c+1:]=(R[r+1:,c+1:]-coeff[:,None]*R[r,c+1:])%p
  r+=1
  if r==n:break
 return r

def stable(word,p):
 M=np.eye(n,dtype=np.int64)
 for a in word:M=mats[a]@M%p
 R=np.eye(n,dtype=np.int64)
 for _ in range(n):R=R@M%p
 return rank(R,p)
rng=random.Random(231)
fail=[]
for i in range(2000):
 word=[rng.randrange(4) for _ in range(rng.randint(1,40))]
 r2=stable(word,2);rq=stable(word,1000003)
 if r2!=rq:fail.append({'word':word,'rank2':r2,'rankp':rq});break
print(json.dumps({'samples':i+1,'failures':fail}))
