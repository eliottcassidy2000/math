from itertools import product
from collections import defaultdict
import numpy as np
T=(1,0,2,1);B=(3,1,3,1)
Q=np.array(((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1)),dtype=np.int64)
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D

def state(w,s=T):
 for c in w:s=step(s,c)
 return s

def words(N):
 for n in range(N+1):
  for w in product('OZJ',repeat=n): yield ''.join(w)

def candidates(seed,total=5):
 out={}
 for na in range(total+1):
  for a in words(na):
   if len(a)!=na:continue
   for b in words(total-na):
    for c in 'OZ':
     for off in range(-3,3):
      s=np.array(state(a+c*(4+off)+b,seed),dtype=np.int64)
      t=np.array(state(a+c*(5+off)+b,seed),dtype=np.int64)
      A=t-s; C=2*s-t
      k=tuple(A)+tuple(C)
      if k not in out:out[k]=(a,c,off,b)
 return list(out.items())

if __name__=='__main__':
 V=candidates(T,4);W=candidates(B,4)
 print('candidates',len(V),len(W),flush=True)
 av=np.array([x[0][:4]for x in V]); bv=np.array([x[0][4:]for x in V]);aw=np.array([x[0][:4]for x in W]);bw=np.array([x[0][4:]for x in W]);matches=[]
 for i in range(len(V)):
  mask=(aw@(Q@av[i])==1024)&(aw@(Q@bv[i])+bw@(Q@av[i])==-64)&(bw@(Q@bv[i])==3)
  for j in np.where(mask)[0]: matches.append((V[i][1],W[j][1]))
 print('balanced matches',len(matches));print(matches[:100])
