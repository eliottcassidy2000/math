#!/usr/bin/env python3
"""All labelled tournaments through order five; exact H=13 classification."""
import hashlib,itertools,json
from pathlib import Path

def need(x,m):
 if not x:raise ValueError(m)

def rows_from_mask(n,mask):
 rows=[0]*n;t=0
 for i in range(n):
  for j in range(i+1,n):
   if mask>>t&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
   t+=1
 return rows

def strong(rows):
 n=len(rows)
 for i in range(n):
  seen=1<<i
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def direct(rows):
 n=len(rows);paths=[p for p in itertools.permutations(range(n)) if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]
 cycles=[(0,)+p for p in itertools.permutations(range(1,n)) if n>=3 and all(rows[u]>>v&1 for u,v in zip((0,)+p,p+(0,)))];s=[0]*n;e=[0]*n
 for p in paths:s[p[0]]+=1;e[p[-1]]+=1
 return paths,cycles,s,e

def dp(rows):
 n=len(rows);N=1<<n;d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    old=M^(1<<v)
    if old:d[M][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
 return sum(d[-1]),d[-1]

def canonical(rows):
 n=len(rows)
 return min(tuple(sum(1<<j for j,v in enumerate(p) if rows[u]>>v&1) for u in p) for p in itertools.permutations(range(n)))

def pad(rows,a,b):
 n=len(rows);N=a+n+b;out=[0]*N
 for i in range(N):
  for j in range(i+1,N):
   if a<=i<a+n and a<=j<a+n:
    u,v=i-a,j-a
    if rows[u]>>v&1:out[i]|=1<<j
    else:out[j]|=1<<i
   else:out[i]|=1<<j
 return out

def main():
 inventory=[];classes={}
 for n in range(1,6):
  total=1<<(n*(n-1)//2);hits=0;maximum=0
  for mask in range(total):
   rows=rows_from_mask(n,mask);paths,cycles,s,e=direct(rows);h=len(paths);z,end=dp(rows)
   need(z==h and end==e,'independent recurrence');maximum=max(maximum,h)
   if h!=13:continue
   hits+=1;need(strong(rows) and len(cycles)==2 and min(s)>=2 and min(e)>=2,'H13 labelled classification')
   key=canonical(rows);classes[key]=classes.get(key,0)+1
  inventory.append({'n':n,'labelled_tournaments':total,'max_H':maximum,'H13_labelled':hits})
 expected=(2,12,9,17,7);need(classes=={expected:120},'unique H13 class')
 paths,cycles,s,e=direct(expected);need(cycles==[(0,1,2,3,4),(0,1,3,4,2)],'explicit cycles')
 need(5**5>13**3,'all-order strong cutoff')
 checks=[]
 for a in range(4):
  for b in range(4):
   rows=pad(expected,a,b);h,ends=dp(rows);rev=[((1<<len(rows))-1)^(1<<i)^r for i,r in enumerate(rows)];_,starts=dp(rev)
   need(h==13 and 1 not in starts and 1 not in ends,'padding endpoint control')
   need(starts==([13]+[0]*(len(rows)-1) if a else s+[0]*b),'exact start padding')
   need(ends==([0]*(len(rows)-1)+[13] if b else [0]*a+e),'exact end padding')
   checks.append({'a':a,'b':b,'s':starts,'e':ends})
 receipt={'status':'PASS','inventory':inventory,'core':{'rows':expected,'paths':paths,'cycles':cycles,'s':s,'e':e,'labelled_orbit_size':120},'padding_checks':checks,'cutoff':{'5^5':5**5,'13^3':13**3}}
 raw=(json.dumps(receipt,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('thirteen_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','labelled_tournaments':sum(x['labelled_tournaments'] for x in inventory),'H13_classes':1,'H13_cycles':2,'padding_checks':16,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
