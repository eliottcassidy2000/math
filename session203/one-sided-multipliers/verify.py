"""Literal graph controls for PROOF.md; no threshold claim."""
import itertools,json,random
from pathlib import Path

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def matrix(a):
 n=len(a);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for S in range(1,1<<n):
  for v in range(n):
   if not D[S][v]:continue
   z=a[v]&~S
   while z:
    b=z&-z;z-=b;w=b.bit_length()-1
    D[S|b][w]+=D[S][v]
 return D

def profile(a):
 n=len(a);N=(1<<n)-1;E=matrix(a)
 rev=[sum(1<<j for j in range(n)if a[j]>>i&1)for i in range(n)]
 S=matrix(rev);H=[sum(x)for x in E];H[0]=1
 R=sum(H[X]*S[N^X][0]for X in range(N+1)if not X&1)
 L=sum(E[X][0]*H[N^X]for X in range(N+1)if X&1)
 return dict(h=H[N],e=E[N][0],s=S[N][0],R=R,L=L,d=H[N^1])

def Ecore(u):
 if u==1:return [2,0]
 m=(u-1).bit_length();a=[2,0]+[0]*m
 for j in range(m):
  v=j+2;a[v]|=1
  if (u-1)>>(m-1-j)&1:a[1]|=1<<v
  else:a[v]|=2
  for k in range(j+1,m):a[v]|=1<<(k+2)
 return a

def join(g,u):
 n=len(g);a=list(g)+[0]*(len(u)-1)
 for j in range(1,len(u)):
  v=n+j-1
  if u[0]>>j&1:a[0]|=1<<v
  else:a[v]|=1
  for i in range(1,n):a[i]|=1<<v
  for k in range(1,j):
   w=n+k-1
   if u[j]>>k&1:a[v]|=1<<w
   else:a[w]|=1<<v
 return a

inputs=[]
for n in range(1,4):
 for code in range(1<<(n*(n-1)//2)):
  u=[0]*n
  for k,(i,j)in enumerate(itertools.combinations(range(n),2)):
   if code>>k&1:u[i]|=1<<j
   else:u[j]|=1<<i
  inputs.append(u)
rng=random.Random(2033)
for n in (4,5):
 for _ in range(6):
  u=[0]*n
  for i in range(n):
   for j in range(i+1,n):
    if rng.randrange(2):u[i]|=1<<j
    else:u[j]|=1<<i
  inputs.append(u)

profiles=[];checks=0;examples=[]
for factor in (1,3,5,9,17,33):
 g=Ecore(factor);p=profile(g)
 need((p['h'],p['e'],p['R'],p['d'])==(2*factor-1,factor-1,2*factor,factor),'core profile')
 profiles.append(dict(factor=factor,rows=g,profile=p))
 for u in inputs:
  up=profile(u);h=sum(matrix(join(g,u))[-1]);base=sum(matrix(join([2,0],u))[-1])
  need(base==up['h']+up['L']-up['s'],'base context')
  need(h==factor*(up['h']+up['L'])-up['s'],'exact multiplier graph')
  need(h==factor*base+(factor-1)*up['s'],'relative context')
  checks+=1
  if factor in (3,9) and len(u)==3 and up['h']==3:examples.append(dict(factor=factor,input=u,input_profile=up,base=base,output=h))

out=dict(status='PASS_LITERAL_GRAPH_CONTROLS',rooted_inputs=len(inputs),graph_pairs=checks,
         profiles=profiles,examples=examples,
         scope='Finite controls for the all-order context formula and inherited zero-slack classification; no coverage theorem')
print(json.dumps(out,indent=2))
Path(__file__).with_name('ONE_SIDED_MULTIPLIER_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
