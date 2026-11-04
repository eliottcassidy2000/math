#!/usr/bin/env python3
"""Self-contained verification of both physical continuations of the three-cycle lower branch."""
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


LEFT={
 'parent':([14,52,40,2,45,9],29,33,3),
 'lower_birth':([78,116,104,2,109,73,8],51,57,3),
 'lower':([206,244,232,130,237,201,8,64],89,99,3),
 'upper_birth':([78,116,104,66,109,9,32],65,73,4),
 'upper_flip':([78,124,104,64,109,9,32],27,27,1),
 'upper':([206,252,104,64,109,137,32,92],91,99,2),
}
def named(core,k):
 rows,w=wrap(core,k);return rows,list(range(len(core)))+['z'+str(i) for i in range(k)]+['w']
def changes(a,b):
 ar,an=a;br,bn=b;index={v:i for i,v in enumerate(bn)};need(set(an)<=set(bn),'old name inclusion');out=[]
 for i,u in enumerate(an):
  for j in range(i+1,len(an)):
   if ((ar[i]>>j)&1)!=((br[index[u]]>>index[an[j]])&1):out.append((u,an[j]))
 return out

def main():
 receipt={'cores':{},'wrappers':[]}
 for name,(rows,h,f,e) in LEFT.items():
  p=profile(rows,1);need(p==direct(rows,1),'literal core comparison')
  need((p['H'],p['F'],p['e'][1])==(h,f,e) and strong(rows),'exact core target')
  receipt['cores'][name]={'rows':rows,**p}
 need(LEFT['lower_birth'][0][6]==8 and LEFT['lower'][0][7]==64,'lower sole exits')
 need(LEFT['upper_birth'][0][6]==32 and LEFT['upper'][0][7]==92,'upper birth outneighbors')
 for k in range(1,7):
  x=1<<k;wrapped={name:named(rows,k) for name,(rows,h,f,e) in LEFT.items()}
  for a,b,d in [('parent','lower_birth',[]),('lower_birth','lower',[]),('parent','upper_birth',[]),('upper_birth','upper_flip',[(1,3)]),('upper_flip','upper',[])]:need(changes(wrapped[a],wrapped[b])==d,'literal named step')
  for name,(rows,names) in wrapped.items():
   core,h,f,e=LEFT[name];actual=sum(path_dp(rows)[-1]);need(actual==h+x*f and strong(rows),'wrapper count/strongness');entry={'name':name,'k':k,'rows':rows,'names':names,'H':actual,'cycles_by_core_endpoint':e}
   if k==1:cs=cycles(rows);need(len(cs)==e,'independent cycle count');entry['cycles']=cs
   receipt['wrappers'].append(entry)
  L=29+33*x;need(89+99*x==3*L+2 and 91+99*x==3*L+4,'ternary targets')
  need(51+57*x<89+99*x and 65+73*x<91+99*x and 27+27*x<91+99*x,'target nonovershoot')
 raw=(json.dumps(receipt,indent=2,sort_keys=True)+'\n').encode();Path(__file__).with_name('left_continuation_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','independent_core_profiles':6,'physical_wrappers':36,'named_move_checks':30,'direct_cycle_checks':6,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
