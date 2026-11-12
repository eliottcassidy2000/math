#!/usr/bin/env python3
"""Literal14-vertex full-central six-set from unit-padding of oneC5 core."""
from mixed import packet,attach,fuse,check
import json

def graph_fuse(A,B):
 n=len(A);m=len(B);T=[0]*(n+m-1)
 for G,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i,u in enumerate(ids):
   for j,v in enumerate(ids):
    if G[i]>>j&1:T[u]|=1<<v
 for i in range(1,n):
  for j in range(n,n+m-1):T[i]|=1<<j
 return tuple(T)
def unitgraph(p,n):
 D=tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
 return attach(D,((1<<n)-1)^p)
def main():
 D=(6,12,24,17,3);specs=[(6689,7,12,7,0),(6691,7,12,11,0),(6695,3,4,7,8),(6697,3,4,11,8),(6701,1,0,3,34),(6703,1,0,5,34)]
 out=[]
 for H,a,p,c,q in specs:
  L=unitgraph(p,a);X=attach(D,31^c);R=unitgraph(q,8-a)
  T=graph_fuse(graph_fuse(L,X),R);pr=packet(T);algebra=fuse(fuse(packet(L),packet(X)),packet(R))
  check(pr==algebra,'independent literal/fusion agreement')
  check(pr==(H,15,30720,24834,8259,0,15),'full-central six-set')
  out.append(dict(H=H,rooted_rows=T,packet=pr,spec=dict(prefix_length=a,incoming_prefix=p,incoming_core=c,incoming_suffix=q)))
 # Exact failed sliding-pump: raw+6 accompanied by signed endpoint defects.
 X=packet(attach(D,31^7));A=packet(unitgraph(1,4));B=packet(unitgraph(8,4))
 lhs=fuse(A,X);rhs=fuse(X,B);delta=[b-a for a,b in zip(lhs,rhs)]
 check(delta==[6,0,0,0,0,6,6],'signed sliding defects')
 L=packet(unitgraph(4,3));R=packet(unitgraph(0,1));pump=[]
 for k in range(1,9):
  bank=[]
  for j in range(k+1):
   P=L
   for _ in range(j):P=fuse(P,A)
   P=fuse(P,X)
   for _ in range(k-j):P=fuse(P,B)
   P=fuse(P,R);bank.append(P)
  check(len({x[1:] for x in bank})==1,'pump central after sandwich')
  check([x[0]-bank[0][0] for x in bank]==[0]+[-6]*k,'onlyfirst shift survives')
  pump.append(dict(k=k,normalized_offsets=[(x[0]-bank[0][0])//2 for x in bank]))
 print(json.dumps(dict(status='PASS',witnesses=out,normalized_digits=[0,1,3,4,6,7],raw_slide_delta=delta,pump_controls=pump,scope='Actual full-central six-set, not six consecutive odds and not a saturatedmod15bank'),sort_keys=True,indent=2))
if __name__=='__main__':main()
