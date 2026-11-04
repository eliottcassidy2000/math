#!/usr/bin/env python3
"""Self-contained exact verification of the one-to-three-cycle ternary fan."""
import hashlib,itertools,json
from pathlib import Path
CORES={'P':[14,20,8,2,13],'A':[14,52,40,2,45,9],'B':[14,52,40,34,45,1]}
EXPECTED={'P':(9,11,1),'A':(29,33,3),'B':(31,33,3)}
R,V=1,3

def need(b,m):
 if not b:raise ValueError(m)

def strong(rows):
 n=len(rows)
 for v in range(n):
  seen=1<<v
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def path_dp(rows,origin=None):
 n=len(rows);N=1<<n;d=[[0]*n for _ in range(N)]
 for v in range(n):
  if origin is None or origin==v:d[1<<v][v]=1
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    old=M^(1<<v)
    if old:d[M][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
 return d

def profile(rows,root):
 n=len(rows);N=1<<n;end=path_dp(rows);tot=[sum(x) for x in end];tot[0]=1
 converse=[((1<<n)-1)^(1<<v)^rows[v] for v in range(n)];starts=path_dp(converse)[-1]
 J=[]
 for v in range(n):
  d=path_dp(rows,v);J.append(sum(d[M][root]*tot[(N-1)^M] for M in range(1,N)))
 return {'H':tot[-1],'F':sum(end[M][root]*tot[(N-1)^M] for M in range(1,N)),'s':starts,'e':end[-1],'J':J}

def direct(rows,root):
 n=len(rows);N=1<<n;paths={}
 for M in range(N):paths[M]=[p for p in itertools.permutations(v for v in range(n) if M>>v&1) if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]
 s=[0]*n;e=[0]*n;J=[0]*n;F=0
 for p in paths[N-1]:s[p[0]]+=1;e[p[-1]]+=1
 for M,words in paths.items():
  for p in words:
   if p and p[-1]==root:F+=len(paths[(N-1)^M]);J[p[0]]+=len(paths[(N-1)^M])
 return {'H':len(paths[N-1]),'F':F,'s':s,'e':e,'J':J}

def cycles(rows):
 n=len(rows)
 return [(0,)+p for p in itertools.permutations(range(1,n)) if all(rows[u]>>v&1 for u,v in zip((0,)+p,p+(0,)))]

def wrap(core,k):
 n=len(core);w=n+k;rows=core+[0]*(k+1);rows[R]|=1<<w
 rows[w]=sum(1<<v for v in range(w) if v!=R)
 for i in range(k):rows[n+i]=((1<<n)-1)|sum(1<<(n+j) for j in range(i+1,k))
 return rows,w

def main():
 receipt={'cores':{},'wrappers':[]}
 for name,rows in CORES.items():
  p=profile(rows,R);need(p==direct(rows,R),'core permutation/DP comparison')
  need((p['H'],p['F'],p['e'][R])==EXPECTED[name] and strong(rows),'core target')
  receipt['cores'][name]={'rows':rows,**p,'cycles':cycles(rows)}
 for name,out in [('A',9),('B',1)]:
  need([r&31 for r in CORES[name][:5]]==CORES['P'] and CORES[name][5]==out,'literal one-birth extension')
 need([(i,j) for i in range(6) for j in range(i+1,6) if (CORES['A'][i]^CORES['B'][i])>>j&1]==[(3,5)],'single physical arc flip')
 for k in range(1,7):
  N=9+11*(1<<k)
  for name,core in CORES.items():
   rows,w=wrap(core,k);d=path_dp(rows);h=sum(d[-1]);rev=[((1<<len(rows))-1)^(1<<v)^r for v,r in enumerate(rows)];s=path_dp(rev)[-1];p=receipt['cores'][name]
   need(strong(rows) and h==p['H']+(1<<k)*p['F'],'physical lifted count/strongness')
   need(s[:len(core)]==p['J'] and s[w]==p['H'],'full old/outer start vector')
   need(s[len(core):w]==[(1<<(k-i-1))*p['F'] for i in range(k)],'clock start vector')
   need(h=={'P':N,'A':3*N+2,'B':3*N+4}[name],'ternary identity')
   need(s[V]=={'P':1,'A':3,'B':4}[name],'same physical start-budget port')
   entry={'name':name,'k':k,'rows':rows,'H':h,'s':s,'cycle_count_by_core_endpoint':p['e'][R]}
   if k<=2:
    cs=cycles(rows);need(len(cs)==p['e'][R],'independent lifted cycle enumeration');entry['cycles']=cs
   receipt['wrappers'].append(entry)
 raw=(json.dumps(receipt,indent=2,sort_keys=True)+'\n').encode();Path(__file__).with_name('cycle_fan_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','independent_core_profiles':3,'physical_wrappers':18,'direct_cycle_checks':6,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
