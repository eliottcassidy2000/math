#!/usr/bin/env python3
"""Exact twelve-value strong windows at 15*2^k, k>=8."""
import json,hashlib,itertools
from pathlib import Path
from verify import require,INIT,F,TRANS,row_step,row_word,state,dot,value,graph,dp
SEEDS=[('2222000000',-13,7),('2032',-11,1),('21300',-9,3),('323',-7,0),
 ('1001133110',-5,7),('2020000010',-3,7),('230',-1,0),('110000022',1,6),
 ('31120',3,3),('3000000012',5,7),('130210',7,4),('1000103210',9,7)]
def strong(r):
 n=len(r)
 for i in range(n):
  seen=1<<i
  while True:
   nxt=seen
   for j in range(n):
    if seen>>j&1:nxt|=r[j]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def literal(r):
 n=len(r);total=0
 def walk(mask,end):
  nonlocal total
  if mask==(1<<n)-1:total+=1;return
  out=r[end]&~mask
  while out:
   b=out&-out;out-=b;walk(mask|b,b.bit_length()-1)
 for i in range(n):walk(1<<i,i)
 return total

def rooted_profile(r,b):
 n=len(r);N=1<<n;ends=[[0]*n for _ in range(N)];starts=[[0]*n for _ in range(N)]
 for i in range(n):ends[1<<i][i]=1
 starts[1<<b][b]=1
 for mask in range(1,N):
  for i in range(n):
   if not(ends[mask][i] or starts[mask][i]):continue
   todo=r[i]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
    ends[mask|bit][j]+=ends[mask][i]
    starts[mask|bit][j]+=starts[mask][i]
 h=[sum(x)for x in ends];h[0]=1
 R=sum(sum(starts[mask])*h[(N-1)^mask] for mask in range(N)if mask>>b&1)
 return h[-1],ends[-1][b],R,h[(N-1)^(1<<b)]

def main():
 r1=row_step(F,1);r2=row_step(r1,1);r3=row_step(r2,1)
 require(all(c-3*b+2*a==0 for a,b,c in zip(r1,r2,r3)),'full-state tail annihilator')
 ar=[];br=[]
 for o,x,y,z in itertools.product(range(2),range(3),range(3),range(2)):
  if o==0:
   br.append(int(x!=2 and z==1))
   ar.append(int(x!=2 and z==0)-int(x!=2 and y==2 and z==1))
  else:
   br.append(int(y==1));ar.append(-int(x==2 and y==1))
 require(all(a==2*x-y for a,x,y in zip(ar,r1,r2)),'combinatorial constant row')
 require(all(2*b==y-x for b,x,y in zip(br,r1,r2)),'combinatorial exponential row')
 records=[]
 for s,a,q in SEEDS:
  u=tuple(map(int,s));v=state(u);h1=dot(r1,v);h2=dot(r2,v)
  require(2*h1-h2==a and h2-h1==30*2**q,'exact tail coefficients')
  h0,e,R,D=rooted_profile(graph(u),1)
  require(h0-e-R+2*D==a and R-D==15*2**q,'rooted two-block coefficients')
  require(u[0]!=0 and (2 in u or any(x==3 and y==0 for i,x in enumerate(u)for y in u[i+1:])),'all-tail strongness')
  for t in (1,2):require(value(u+(1,)*t)==dp(graph(u+(1,)*t)),'small physical DP')
  if len(u)<=6:require(h1==literal(graph(u+(1,))),'literal small seed')
  k=8;w=u+(1,)*(k-q);h=15*2**k+a
  require(value(w)==h==dp(graph(w)),'actual k8 witness')
  require(strong(graph(w)),'actual k8 strongness')
  require(len(w)+2-k in (4,5),'order offset')
  records.append({'prefix':s,'offset':a,'q':q,'slope':15*2**q,'H_tail1':h1,
   'H_tail2':h2,'rooted_profile':{'h':h0,'e':e,'R':R,'D':D},'order_minus_k':len(u)+2-q,'k8_word':''.join(map(str,w)),'k8_H':h})
 require([a for s,a,q in SEEDS]==list(range(-13,10,2)),'consecutive bank')
 result={'status':'PASS','k_min':8,'center':'15*2^k','odd_offsets':list(range(-13,10,2)),
 'seed_count':12,'rooted_profiles':12,'combinatorial_coefficient_coordinates':72,'physical_DP_controls':36,'literal_seed_controls':sum(len(s)<=6 for s,a,q in SEEDS),'records':records}
 raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
 Path(__file__).with_name('window15_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','strong_windows':'all odd H from 15*2^k-13 to 15*2^k+9, k>=8',
  'physical_DP_controls':36,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
