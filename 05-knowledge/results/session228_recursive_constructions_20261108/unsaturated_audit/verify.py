#!/usr/bin/env python3
"""Independent standard-library audit; imports no producer code; no assertions."""
import itertools,json,math
from pathlib import Path

def ck(v,s):
 if not v: raise RuntimeError(s)
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def mv(a,v):return [sum(x*y for x,y in zip(r,v)) for r in a]
def diff(a,b):return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def endpoints(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v);d[mask][v]=sum(d[rest][u] for u in range(n) if a[u]>>v&1) if rest else 1
 return d

def local(a,out):
 n=len(a);full=(1<<n)-1;e=endpoints(a);rev=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)];f=endpoints(rev)
 h=[sum(r) for r in e];h[0]=1;m=[[0]*6 for _ in range(6)]
 for p in range(full+1):
  q=full^p
  for col in range(6):
   oldlast,begun=divmod(col,2)
   qcount=h[q] if begun or not q else sum(f[q][j] for j in range(n) if out>>j&1)
   if not p:m[2*oldlast+int(bool(begun or q))][col]+=qcount
   else:
    for j in range(n):m[2*(2 if out>>j&1 else 1)+int(bool(begun or q))][col]+=e[p][j]*qcount
 alpha=[sum(f[q][j]*h[full^q] for q in range(1,full+1)) for j in range(n)]
 return m,h[-1],sum(h[p]*h[full^p] for p in range(full+1)),alpha

def graph(a,masks,bits):
 blocks=[([0],1)]
 for i,b in enumerate(bits):
  blocks.append((a,masks[b]))
  if i<len(bits)-1:blocks.append(([0],1))
 blocks.append(([0],0));n=1+sum(len(b) for b,m in blocks);adj=[0]*n;off=1
 for b,m in blocks:
  for u,row in enumerate(b):
   adj[off+u]=row<<off
   for v in range(off+len(b),n):adj[off+u]|=1<<v
   if m>>u&1:adj[0]|=1<<(off+u)
   else:adj[off+u]|=1
  off+=len(b)
 return adj

def literal(a):
 n=len(a);full=(1<<n)-1;e=endpoints(a);h=[sum(r) for r in e];h[0]=1
 # Root-started paths are counted with a separate forward DP.
 r=[[0]*n for _ in range(1<<n)];r[1][0]=1
 for s in range(1,full+1,2):
  for v,z in enumerate(r[s]):
   if not z:continue
   available=a[v]&~s&full
   while available:
    bit=available&-available;available-=bit;u=bit.bit_length()-1;r[s|bit][u]+=z
 U=sum(sum(r[s])*h[full^s] for s in range(1,full+1,2))
 return [h[full],e[full][0],U,h[full^1]]

def strong(a):
 for v in range(len(a)):
  reached=1<<v
  while True:
   more=reached
   for u in range(len(a)):
    if reached>>u&1:more|=a[u]
   if more==reached:break
   reached=more
  if reached!=(1<<len(a))-1:return False
 return True

def main():
 p=json.loads((Path(__file__).parent/'inputs/frozen.json').read_text());a=p['adjacency'];masks=p['root_out_masks']
 records=[local(a,m) for m in masks];M=[r[0] for r in records];S=[local([0],b)[0] for b in range(2)];F=diff(M[1],M[0]);v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
 ck(records[0][1:]==(19,256,[21,20,29,35,73,59]),'block invariants')
 ck(strong(a),'strong block');ck([sum(records[0][3][j] for j in range(6) if m>>j&1) for m in masks]==[94,94],'alpha collision')
 neg=local(p['negative_adjacency'],0);sums=[sum(neg[3][j] for j in range(4) if m>>j&1) for m in range(16)]
 ck(neg[1:]==(5,32,[5,6,9,7]) and len(set(sums))==16,'negative control')
 negative_blocks=0
 for code in range(64):
  adj=[0]*4
  for i,(u,w) in enumerate(itertools.combinations(range(4),2)):
   if code>>i&1:adj[u]|=1<<w
   else:adj[w]|=1<<u
  rec=local(adj,0)
  if rec[1:3]==(5,32):
   negative_blocks+=1;ck(len({sum(rec[3][j] for j in range(4) if m>>j&1) for m in range(16)})==16,'all order4 controls')
 ck(negative_blocks==24,'order4 labeled control count')
 def count(bits):
  x=mv(S[1],v)
  for i,b in enumerate(bits):
   x=mv(M[b],x)
   if i<len(bits)-1:x=mv(S[1],x)
  return sum(x*y for x,y in zip(ell,mv(S[0],x)))
 literalrows=[]
 for t in (1,2):
  for bits in itertools.product(range(2),repeat=t):
   adj=graph(a,masks,bits);r=literal(adj);ck(r==[count(bits),19**t,(700*512**t-207*19**t)//493,19**t],'literal full profile');ck(strong(adj),'literal strong');literalrows.append({'bits':bits,'profile':r})
  rows=[r['profile'][1:] for r in literalrows if len(r['bits'])==t];ck(len(set(map(tuple,rows)))==1,'literal sidecar')
 # Exact matrix certificates imply all-order recurrences, not finite extrapolation.
 A=mm(S[1],M[0]);seed=mv(S[1],v)
 def shift(c):return [[A[i][j]-c*int(i==j) for j in range(6)] for i in range(6)]
 P2=mm(shift(512),shift(19));P3=mm(P2,shift(19))
 ck(mv(P3,seed)==[0]*6,'baseline recurrence certificate')
 ck(not any(map(any,mm(F,P2))),'defect recurrence certificate')
 K=[[0,0,1,0,-1,0],[0,0,0,1,0,-1]]
 ck(all(mv(F,k)==[0]*6 for k in K),'F kills K')
 ck(all(all(z[0]==z[1]==0 and z[2]==-z[4] and z[3]==-z[5] for z in [mv(B,k)]) for B in M+S for k in K),'K invariant')
 def exact(a,b):
  ck(a%b==0,'exact divisibility');return a//b
 def base(t):return exact(19950*512**t+(13862*t-2289)*19**t,15979)
 def weight(t,i):return exact(-(8400*512**i*19**(t-1-i)+(474 if i==t-1 else 1460)*19**(t-1)),493)
 for t in range(1,4):ck(count((0,)*t)==base(t),'baseline initial values')
 for i in (0,1):
  z=mv(F,seed if i==0 else mv(A,seed))
  ck(z[2]==-2*19**i and z[3]==-(8400*512**i+474*19**i)//493,'defect initial values')
 # Universal category-basis audit, including arbitrary intervening blocks.
 def category(p,r):
  C=[[0]*6 for _ in range(6)]
  for old in range(6):
   u,b=divmod(old,2)
   if not b and r==2:continue
   C[2*(p or u)+int(bool(b or r))][old]=1
  return C
 W=[[0,1,0,0,0,-1],[0,0,0,1,0,-1]];V=W+[[0,0,1,0,-1,0]]
 def inv(z):return z[0]==0 and z[2]+z[4]==0 and z[1]+z[3]+z[5]==0
 cats={(p,r):category(p,r) for p,r in itertools.product(range(3),repeat=2) if p or r}
 ck(all(inv(mv(C,z)) for C in cats.values() for z in V),'all eight category maps preserve V')
 ck(all(inv(mv(C,z)) and mv(C,z)[2]==mv(C,z)[4]==0 for C in S for z in V),'singletons V to W')
 universal=[]
 for j in range(4):
  x,y,a,b=[int(i==j) for i in range(4)];coeff={(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b}
  C=[[sum(c*cats[k][r][s] for k,c in coeff.items()) for s in range(6)] for r in range(6)]
  ck(all(inv(list(z)) for z in zip(*C)),'universal defect image V');ck(all(mv(C,z)==[0]*6 for z in W),'universal defect kills W');universal.append(C)
 ck(any(any(map(any,mm(C,D))) for C in universal for D in universal),'empty separation counterexample')
 ck(8400*(512-20*19)-19*(20*1460-474)>0,'worst final superincreasing comparison')
 finite=[]
 for t in range(1,11):
  vals=[count(b) for b in itertools.product(range(2),repeat=t)];ck(vals==[base(t)+sum(b[i]*weight(t,i) for i in range(t)) for b in itertools.product(range(2),repeat=t)],'all binary formulas');ck(len(set(vals))==2**t,'finite distinctness');ck(math.gcd(*vals)==1,'finite gcd');finite.append({'t':t,'base':vals[0],'gcd':math.gcd(*vals),'distinct':len(set(vals)),'weights':[count(tuple(int(i==j) for i in range(t)))-vals[0] for j in range(t)]})
 print(json.dumps({'status':'PASS','M':M,'S':S,'F':F,'literal':literalrows,'finite':finite,'negative_subset_sums':sums,'negative_labeled_blocks':negative_blocks,'scope':'Finite literal checks plus exact all-order matrix certificates; see AUDIT.md. No prime or interval coverage claim.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
