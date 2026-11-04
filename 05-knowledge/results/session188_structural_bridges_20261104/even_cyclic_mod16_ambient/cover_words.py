from collections import defaultdict
from itertools import product
import json
from pathlib import Path
R=Path(__file__).resolve().parent

def tensor(a,c,reverse=False,maxh=5):
 Q=[4,1,2] if reverse else [2,4,1];adj=list(Q)+[16,0]
 for i in range(3):
  if a>>i&1:adj[i]|=8
  else:adj[3]|=1<<i
  if c>>i&1:adj[i]|=16
  else:adj[4]|=1<<i
 states={(0,0,0,0,-1,0):1};out=defaultdict(int)
 for length in range(18):
  nxt=defaultdict(int)
  for (i,j,k,m,last,b),v in states.items():
   cnt=(i,j,k)
   if m==3 and min(cnt)>0:out[b+1,i,j,k]=out[b+1,i,j,k]+v
   for z in range(5):
    if z<3:
     if cnt[z]==5:continue
     cc=list(cnt);cc[z]+=1;mm=m
    else:
     bit=1<<(z-3)
     if m&bit:continue
     cc=cnt;mm=m|bit
    key=tuple(cc)+(mm,z,b)
    if last<0 or adj[last]>>z&1:nxt[key]=nxt[key]+v
    if last>=0 and b<maxh-1:
     key=tuple(cc)+(mm,z,b+1);nxt[key]=nxt[key]+v
  states={s:v for s,v in nxt.items() if v}
 return dict(out)

