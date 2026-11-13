from functools import cache
from itertools import product
B=(2,80,3,7,173,143,189,15); ALL=255
@cache
def end(S,v):
 if not (S>>v&1): return 0
 if S==1<<v:return 1
 return sum(end(S^(1<<v),j) for j in range(8) if S>>j&1 and B[j]>>v&1)%3
@cache
def start(S,v):
 if not (S>>v&1): return 0
 if S==1<<v:return 1
 return sum(start(S^(1<<v),j) for j in range(8) if S>>j&1 and B[v]>>j&1)%3
@cache
def F(S,r):
 if not r:return int(not S)
 a=S;z=0
 while a:
  z+=sum(end(a,v) for v in range(8))*F(S^a,r-1)
  a=(a-1)&S
 return z%3
@cache
def base(r,loc,v):
 if not loc:return F(ALL,r)
 if len({i for i,j in loc})>1:return 0
 if len(loc)==2:return F(ALL^(1<<v),r-1)
 f=start if loc[0][1]==0 else end
 a=ALL;z=0
 while a:
  z+=f(a,v)*F(ALL^a,r-1);a=(a-1)&ALL
 return z%3
@cache
def walks(r,lim,empt=False):
 W=[]
 def go(a,n,left):
  if not n:
   W.append(a);return
  if empt:go(a+((),),n-1,left)
  for s in range(3):
   l=list(left);w=()
   for k in range(sum(l)):
    j=(s+k)%3;l[j]-=1;w=w+(j,)
    if l[j]<0:break
    go(a+(w,),n-1,tuple(l))
 go((),r,lim)
 return W
def extract(ws):
 R=[0]*3; ep={}
 for i,w in enumerate(ws):
  for k,j in enumerate(w):
   if k==0:ep[i,0]=(j,R[j],0)
   if k==len(w)-1:ep[i,1]=(j,R[j],1)
   R[j]+=1
 return R,ep
@cache
def mid(r,loc,v):
 if r==0 or r>=3+int(bool(loc)):return 0
 z=0
 for ws in walks(r,(2+int(bool(loc)),2,2)):
  rs,ep=extract(ws);req=[ep[x] for x in loc]
  if any(i!=0 for i,j,k in req):continue
  z+=base(rs[0],tuple(sorted((j,k) for i,j,k in req)),v)*F(ALL,rs[1])*F(ALL,rs[2])
 return z%3
def count(colors,v):
 z=0
 for ws in walks(3,(3,3,3),True):
  rs,ep=extract(ws)
  if not all(rs):continue
  for order in ('xy','yx'):
   if order=='yx' and not ws[1]:continue
   required=set();forbid=set();legal=True
   for i,x in enumerate(order):
    for a,d in ((i,1),(i+1,0)):
     t=ep.get((a,d))
     if t is None:continue
     j,r,s=t
     c=colors[j]
     wants=(x=='x' and d==1) or (x=='y' and d==0)
     matches=(x=='x' and c=='R') or (x=='y' and c=='B')
     if wants:
      if not matches:legal=False;break
      required.add(t)
     elif matches:
      forbid.add(t)
    if not legal:break
   if not legal or required&forbid:continue
   forbid=tuple(forbid)
   for mask in range(1<<len(forbid)):
    req=required|{t for k,t in enumerate(forbid) if mask>>k&1}
    q=1
    for j in range(3):
     q*=mid(rs[j],tuple(sorted((r,s) for jj,r,s in req if jj==j)),v)
    z+=(-1)**mask.bit_count()*q
 return z%3
print([(v,count('RBB',v),count('RRB',v)) for v in range(8)])
