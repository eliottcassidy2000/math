"""Exact one-apex modular synthesis; standard library only."""
import json
from collections import deque
from itertools import product

def need(ok,msg):
 if not ok: raise RuntimeError(msg)
def state(w):
 s=p=0; h=1
 for b in w:
  h+=(1-b)*(s+p);s,p=2*s+b,b
 return s,p,h

def cycles(w):
 return 1+2*sum((1 if j-i==1 else 2**(j-i-2)) for i in range(len(w)) for j in range(i+1,len(w)) if w[i]==1 and w[j]==0)
def dp(w):
 n=len(w)+1; out=[0]*n
 for i,b in enumerate(w,1):
  if b:out[0]|=1<<i
  else:out[i]|=1
 for i in range(1,n):
  for j in range(i+1,n):out[i]|=1<<j
 D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if D[mask][v]:
    avail=out[v]&~mask
    while avail:
     z=avail&-avail;avail-=z;u=z.bit_length()-1
     D[mask|z][u]+=D[mask][v]
 return sum(D[-1])
def step(x,h,label,q):
 if label=='1':return (2*x%q,h)
 if label=='01':return ((4*x-2)%q,(h+x)%q)
 if label=='001':return ((8*x-6)%q,(h+3*x-2)%q)
 raise ValueError(label)
def bfs(q):
 root=(2%q,1%q);seen={root:None};Q=deque([root])
 while Q:
  z=Q.popleft()
  for lab in ('1','01','001'):
   v=step(*z,lab,q)
   if v not in seen:seen[v]=(z,lab);Q.append(v)
 need(len(seen)==q*q,('transitivity',q,len(seen)))
 witnesses={}
 for r in range(q):
  v=(2%q,(r-2)%q);blocks=[]
  while seen[v] is not None:
   v,lab=seen[v];blocks.append(lab)
  w='1'+''.join(reversed(blocks))+'0';bits=tuple(map(int,w));h=state(bits)[2]
  need(h%q==r,('residue',q,r,w,h));need(h==cycles(bits),'cycles')
  need(len(w)+1<=3*q*q,'length')
  witnesses[r]=(w,h)
 return witnesses

def main():
 count=0
 for m in range(10):
  for w in product((0,1),repeat=m):
   h=state(w)[2]
   need(h==cycles(w)==dp(w),('DP',w,h));count+=1
 rows=[]
 for q in range(1,102,2):
  W=bfs(q);rows.append({'q':q,'states':q*q,'max_order':max(len(w)+1 for w,h in W.values())})
 # Check symbolic affine identities over every odd modulus 3..25.
 for q in range(3,26,2):
  inv4=pow(4,-1,q);inv8=pow(8,-1,q)
  C=lambda x,h:((x-2)%q,(h+x*inv4)%q)
  Ci=lambda x,h:((x+2)%q,(h-(x+2)*inv4)%q)
  D=lambda x,h:((x-6)%q,(h+3*x*inv8-2)%q)
  Di=lambda x,h:((x+6)%q,(h-3*x*inv8-inv4)%q)
  for x,h in product(range(q),repeat=2):
   z=Ci(x,h);z=Di(*z);z=C(*z);z=D(*z)
   need(z==(x,(h+3*inv4)%q),('commutator',q,x,h,z))
   z=Di(x,h)
   for _ in range(3):z=C(*z)
   need(z==(x,(h+(3*x+22)*inv8)%q),('E',q,x,h,z))
 print(json.dumps({'direct_dp_words':count,'modular_bfs':rows,'witnesses_q7':bfs(7),'status':'PASS'},indent=2))
if __name__=='__main__':main()
