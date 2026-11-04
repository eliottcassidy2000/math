#!/usr/bin/env python3
"""Self-contained verification of both physical continuations of the three-cycle upper branch."""
import hashlib,itertools,json
from pathlib import Path

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


B=[14,52,40,34,45,1]
D=[78,52,104,34,109,65,10]
C=[10,52,41,34,45,1]
E=[10,116,41,34,45,1,61]
F=[10,244,41,34,45,1,61,125]
STATES={'parent':(B,(1,)), 'bridge':(B,(1,4)), 'reroot':(B,(4,)), 'lower':(D,(4,)), 'flip':(C,(1,)), 'birth':(E,(1,)), 'upper':(F,(1,))}

def outer(core,S,k):
 n=len(core);w=n+k;rows=core+[0]*(k+1)
 for r in S:rows[r]|=1<<w
 rows[w]=sum(1<<v for v in range(w) if v not in S)
 for i in range(k):rows[n+i]=((1<<n)-1)|sum(1<<(n+j) for j in range(i+1,k))
 return rows,list(range(n))+['z'+str(i) for i in range(k)]+['w']

def changed(a,b):
 rows,names=a;rr,nn=b;idx={v:i for i,v in enumerate(nn)}
 need(set(names)<=set(nn),'old name inclusion')
 out=[]
 for i,u in enumerate(names):
  for j in range(i+1,len(names)):
   v=names[j]
   if ((rows[i]>>j)&1)!=((rr[idx[u]]>>idx[v])&1):out.append((u,v))
 return out

def cross(rows,q,u):
 n=len(rows);count=0
 for word in itertools.permutations(range(n)):
  for cut in range(1,n):
   a,b=word[:cut],word[cut:]
   if a[-1]==q and b[0]==u and all(rows[x]>>y&1 for x,y in zip(a,a[1:])) and all(rows[x]>>y&1 for x,y in zip(b,b[1:])):count+=1
 return count

def main():
 receipt={'cores':{},'wrappers':[],'cross_1_4':cross(B,1,4),'cross_4_1':cross(B,4,1)}
 need((receipt['cross_1_4'],receipt['cross_4_1'])==(19,3),'cross-cover constants')
 specs={'B1':(B,1,31,33,3),'B4':(B,4,31,34,3),'D4':(D,4,95,99,7),'C1':(C,1,25,27,2),'E1':(E,1,49,51,2),'F1':(F,1,97,99,2)}
 for name,(rows,r,h,f,e) in specs.items():
  p=profile(rows,r);need(p==direct(rows,r),'core permutation/DP check')
  need((p['H'],p['F'],p['e'][r])==(h,f,e) and strong(rows),'core target')
  receipt['cores'][name]={'rows':rows,'root':r,**p}
 need([x&63 for x in D[:6]]==B and D[6]==10,'lower birth with outneighbors1,3')
 need([x&63 for x in E[:6]]==C and E[6]==61,'first upper birth')
 need([x&127 for x in F[:7]]==E and F[7]==125,'second upper birth')
 for k in range(1,7):
  x=1<<k;expect={'parent':31+33*x,'bridge':9+67*x,'reroot':31+34*x,'lower':95+99*x,'flip':25+27*x,'birth':49+51*x,'upper':97+99*x};ec={'parent':3,'bridge':6,'reroot':3,'lower':7,'flip':2,'birth':2,'upper':2}
  wrapped={name:outer(core,S,k) for name,(core,S) in STATES.items()}
  for a,b,difference in [('parent','bridge',[(4,'w')]),('bridge','reroot',[(1,'w')]),('reroot','lower',[]),('parent','flip',[(0,2)]),('flip','birth',[]),('birth','upper',[])]:need(changed(wrapped[a],wrapped[b])==difference,'physical named move')
  for name,(rows,names) in wrapped.items():
   h=sum(path_dp(rows)[-1]);need(strong(rows) and h==expect[name],'physical count/strongness')
   core,S=STATES[name];p=profile(core,S[0]);cyc=sum(p['e'][r] for r in S);need(cyc==ec[name],'cycle endpoint law')
   entry={'name':name,'k':k,'rows':rows,'names':names,'H':h,'cycles_by_core_endpoints':cyc}
   if k==1:cs=cycles(rows);need(len(cs)==cyc,'independent cycle enumeration');entry['cycles']=cs
   receipt['wrappers'].append(entry)
  M=expect['parent'];need(expect['lower']==3*M+2 and expect['upper']==3*M+4,'exact ternary children')
  need(all(expect[z]<expect['lower'] for z in ['parent','bridge','reroot']),'lower no target overshoot')
  need(all(expect[z]<expect['upper'] for z in ['parent','flip','birth']),'upper no target overshoot')
 raw=(json.dumps(receipt,indent=2,sort_keys=True)+'\n').encode();Path(__file__).with_name('reroot_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','independent_core_profiles':6,'physical_wrappers':42,'named_move_checks':36,'direct_cycle_checks':7,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
