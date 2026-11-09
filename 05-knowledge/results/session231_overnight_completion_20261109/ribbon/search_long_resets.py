import numpy as np,json,os
from pathlib import Path
p=Path(os.environ.get('MATH_REPO',str(Path(__file__).resolve().parents[4])))/'04-computation/experiments/ribbon_spectrum_20261104/span2_intrep.txt'
a=list(map(int,p.read_text().split()));it=iter(a);n=next(it);q=next(it);y=[next(it) for _ in range(n)];h=[next(it) for _ in range(n)];mats=[];masks=[]
for _ in range(q):
 mask=next(it);masks.append(mask);mats.append(np.array([[next(it) for _ in range(n)] for _ in range(n)],dtype=np.int64))
p=1000003
def singular(M):
 R=M.copy()%p
 for c in range(n):
  nz=np.flatnonzero(R[c:,c])
  if not len(nz):return True
  j=c+int(nz[0]);R[[c,j]]=R[[j,c]]
  coeff=R[c+1:,c]*pow(int(R[c,c]),p-2,p)%p
  R[c+1:,c+1:] = (R[c+1:,c+1:] - coeff[:,None]*R[c,c+1:])%p
 return False
front=[((),np.eye(n,dtype=np.int64))];found=[];tested=0
for length in range(1,10):
 front=[(w+(i,),mats[i]@M) for w,M in front for i in range(1,q)]
 for w,M in front:
  tested+=1
  for radix in (2,4,8):
   if singular(M-radix*np.eye(n,dtype=np.int64)):found.append([[masks[i] for i in w],radix])
 print(json.dumps({'span':2,'alphabet':[1,2,3],'length':length,'tested':tested,'modular_singular_candidates':found}),flush=True)
