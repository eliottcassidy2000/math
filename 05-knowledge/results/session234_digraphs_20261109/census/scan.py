import itertools,json,collections,time
from pathlib import Path

def matrix(n,code):
 a=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n):
   s=code%3;code//=3
   if s in (0,2):a[i][j]=1
   if s in (1,2):a[j][i]=1
 return a

def count(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for j in range(n):
   if s>>j&1:
    p=s^(1<<j)
    if p:d[s][j]=sum(d[p][i] for i in range(n) if p>>i&1 and a[i][j])
 return sum(d[-1])
def direct(a):
 return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))
def canonical(a):
 n=len(a)
 return min(''.join(str(a[p[i]][p[j]]) for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))
def scc(a):
 n=len(a);r=[row[:] for row in a]
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k]&r[k][j]
 unseen=set(range(n));cs=[]
 while unseen:
  i=min(unseen);c={j for j in unseen if i==j or r[i][j] and r[j][i]};unseen-=c;cs.append(c)
 cs.sort(key=lambda c:-sum(a[i][j] for i in c for j in range(n) if j not in c))
 return [(len(c),count([[a[i][j] for j in sorted(c)] for i in sorted(c)])) for c in cs]
out={}
for n in range(1,6):
 hist=collections.Counter();examples={};classes={7:{},21:{}};mismatches=[]
 for code in range(3**(n*(n-1)//2)):
  a=matrix(n,code);h=count(a);hist[h]+=1
  # independent direct permutations across every labelled graph
  hd=direct(a)
  if h!=hd:raise AssertionError((n,code,h,hd))
  if h in classes:
   c=canonical(a)
   if c not in classes[h]:classes[h][c]={'matrix':a,'labelled_count':0,'scc':scc(a),'code':code}
   classes[h][c]['labelled_count']+=1
 out[n]={'labelled_total':sum(hist.values()),'histogram':dict(sorted(hist.items())),'target_classes':classes,'direct_verification':'all labelled graphs'}
 print(n,sum(hist.values()),'targets',hist[7],hist[21],'classes',len(classes[7]),len(classes[21]),flush=True)
 Path('results.json').write_text(json.dumps(out,indent=2))
