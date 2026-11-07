import numpy as np,itertools,time,json
from pathlib import Path

def load(w):
 x=list(map(int,Path(__file__).with_name(f'span{w}_intrep.txt').read_text().split()));r,l=x[:2];v=x[2:2+r];h=x[2+r:2+2*r];pos=2+2*r; A={}
 for _ in range(l):
  mask=x[pos];pos+=1;A[mask]=np.array(x[pos:pos+r*r],dtype=np.int64).reshape(r,r);pos+=r*r
 return A,np.array(v),np.array(h)
def singular(M,p=1009):
 a=M.copy()%p;n=len(a)
 for k in range(n):
  nz=np.nonzero(a[k:,k])[0]
  if not len(nz):return True
  j=k+nz[0];a[[k,j]]=a[[j,k]]
  v=pow(int(a[k,k]),-1,p); z=a[k+1:,k]*v%p
  a[k+1:,k:]=(a[k+1:,k:]-z[:,None]*a[k,k:])%p
 return False
if __name__=='__main__':
 import sys
 w=int(sys.argv[1]);d=int(sys.argv[2]);A,_,_=load(w);n=len(next(iter(A.values())));I=np.eye(n,dtype=np.int64);got=[];start=time.time()
 def visit(B,word):
  if word and singular(B-2*I):
   got.append(word);print('candidate',word,flush=True)
  if len(word)==d:return
  for c,C in A.items():visit(C@B%1009,word+[c])
 visit(I,[])
 print('done',w,d,len(got),time.time()-start)
 Path(__file__).with_name(f'screen{w}_{d}.json').write_text(json.dumps(got))
