#!/usr/bin/env python3
"""Exact finite controls for the three text notes. No assertions are disabled by -O."""
from fractions import Fraction as Q
from decimal import Decimal, localcontext
from itertools import product
from math import comb,gcd
from collections import deque
import hashlib,json
CHECKS=0
def check(v,msg):
 global CHECKS
 CHECKS+=1
 if not v: raise RuntimeError(msg)
MATS={'L':(2,-1,1,0),'M':(2,1,1,0),'R':(1,2,0,1)}
def mul(x,y):
 a,b,c,d=x;e,f,g,h=y
 return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
def val(mat,x):
 a,b,c,d=mat
 return (a*x+b)/(c*x+d)
def det(mat):
 a,b,c,d=mat;return a*d-b*c
def pnum(x):
 n,m=x.numerator,x.denominator
 return n**4-8*n**3*m+2*n*n*m*m+8*n*m**3+m**4
def cmp_mu(x):
 if x<=1:return -1
 if x>=2:return 1
 p=pnum(x);check(p!=0,'rational mu wall')
 return -1 if p>0 else 1
def decode(x):
 original=x;w=[]
 while x!=2:
  check(x>1 and x!=3,'Euclid decoder domain')
  old=x.numerator+x.denominator
  if x<2:w.append('L');x=1/(2-x)
  elif x<3:w.append('M');x=1/(x-2)
  else:w.append('R');x=x-2
  check(x.numerator+x.denominator<old,'Euclid descent')
  check((x.numerator+x.denominator)%2==1,'opposite parity')
 y=Q(2)
 for s in reversed(w):y=val(MATS[s],y)
 check(y==original,'finite suffix reconstruction')
 return ''.join(w)
def fieldmul(x,y):
 # basis 1,sqrt2,sqrt3,sqrt6, indexed by two exponent bits
 z=[0]*4
 for i,a in enumerate(x):
  for j,b in enumerate(y):
   overlap=i&j;factor=(2 if overlap&1 else 1)*(3 if overlap&2 else 1)
   z[i^j]+=a*b*factor
 return tuple(z)
def berggren_controls():
 mu=(2,-1,-1,1);z=(1,0,0,0)
 for coefficient in [-8,2,8,1]:
  z=fieldmul(z,mu);z=(z[0]+coefficient,z[1],z[2],z[3])
 check(z==(0,0,0,0),'quartic exact field identity')
 conjugates={(2,-s,-t,s*t) for s in [-1,1] for t in [-1,1]}
 check(len(conjugates)==4,'four distinct conjugates')
 states=[];matrix=MATS['L'];word='L'
 for _ in range(64):
  states.append((word,matrix))
  c2=det(matrix)*cmp_mu(val(matrix,Q(2)))
  if c2>0:s='L'
  else:s='M' if det(matrix)*cmp_mu(val(matrix,Q(3)))>0 else 'R'
  matrix=mul(matrix,MATS[s]);word+=s
 def side(i,q):
  mat=states[i][1]
  return det(mat)*cmp_mu(val(mat,q)) # sign(q-alpha_i)
 with localcontext() as ctx:
  ctx.prec=250
  muD=Decimal(2)-Decimal(3).sqrt()+Decimal(6).sqrt()-Decimal(2).sqrt()
  approx=[]
  for _,(a,b,c,d) in states:approx.append((Decimal(d)*muD-b)/(a-Decimal(c)*muD))
  order=sorted(range(len(states)),key=lambda i:approx[i])
  samples=[]
  for i,j in zip(order,order[1:]):
   midpoint=(approx[i]+approx[j])/2
   for exponent in range(1,300):
    denominator=2**exponent;numerator=int(midpoint*denominator)
    if numerator%2==0:numerator+=1
    q=Q(numerator,denominator)
    if q>1 and side(i,q)>0 and side(j,q)<0:
     samples.append(q);break
   else:raise RuntimeError('dyadic separator not found')
  for exponent in range(1,300):
   q=Q(2**exponent+1,2**exponent)
   if all(side(i,q)<0 for i in range(64)):samples.insert(0,q);break
  q=Q(2*int(max(approx))+5,2)
  check(all(side(i,q)>0 for i in range(64)),'upper separator')
  samples.append(q)
 suffixes=[decode(q) for q in samples]
 signatures=[]
 for i,(w,mat) in enumerate(states):
  sig=tuple(pnum(val(mat,q))>0 for q in samples)
  signatures.append(sig)
 for i in range(64):
  for j in range(i):check(signatures[i]!=signatures[j],'distinct finite-word residuals')
 # Complete chronological ternary tree to depth 9.
 nodes=[(1,2)];visited=0
 for depth in range(10):
  new=[]
  for m,n in nodes:
   visited+=1;a=n*n-m*m;b=2*m*n;c=n*n+m*m
   check(c*c-4*a*b==pnum(Q(n,m)), 'Pythagorean quartic identity')
   check(gcd(m,n)==1 and (m+n)%2==1,'primitive Euclid parameters')
   t=Q(n,m)
   for letter in ['L','R']:
    y=t
    for _ in range(4):y=val(MATS[letter],y)
    check(pnum(y)>0,'four-step universal cone entry')
   new.extend([(n,2*n-m),(n,2*n+m),(m,n+2*m)])
  nodes=new
 return {'boundary_prefixes':64,'distinguished_pairs':2016,'finite_suffixes':len(samples),
         'max_suffix_word_length':max(map(len,suffixes)),'ternary_nodes_checked':visited,
         'boundary_itinerary_64':states[-1][0],
         'max_matrix_entry_bits':max(abs(v).bit_length() for _,mat in states for v in mat),
         'signature_sha256':hashlib.sha256(repr(signatures).encode()).hexdigest(),
         'suffix_samples':[{'slope':str(q),'word':w} for q,w in zip(samples,suffixes)]}
M=12
FIXED={1:(5,6),2:(1,3),3:(8,11),4:(3,7),5:(0,5),7:(4,11),8:(1,9),10:(0,10),12:(0,12)}
FREE=(6,9,11)
def rank(edges):
 parent=list(range(M+1))
 def root(v):
  while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
  return v
 r=0
 for x,y in edges:
  x,y=root(x),root(y)
  if x!=y:parent[x]=y;r+=1
 return r

def release_controls():
 rows=[]
 for release in [None]+sorted(FIXED):
  fixed=[e for d,e in FIXED.items() if d!=release]
  colors=sorted(FREE+(() if release is None else (release,)))
  states=[]
  for es in product(*[[(i,i+d) for i in range(M+1-d)] for d in colors]):
   if rank(fixed+list(es))==M:states.append(es)
  adj=[[] for _ in states]
  for a in range(len(states)):
   for b in range(a):
    if sum(e!=f for e,f in zip(states[a],states[b]))<=2:adj[a].append(b);adj[b].append(a)
  remaining=set(range(len(states)));components=[]
  while remaining:
   todo=[remaining.pop()];comp=set(todo)
   while todo:
    for w in adj[todo.pop()]:
     if w in remaining:remaining.remove(w);comp.add(w);todo.append(w)
   components.append(comp)
  ranks={d:rank(fixed+[(i,i+d) for i in range(M+1-d)])-len(fixed) for d in colors}
  old=[i for i,s in enumerate(states) if release is None or s[colors.index(release)]==FIXED[release]]
  check(len(old)==5,'original fiber retained')
  if release is not None:
   for d in FREE:
    Ed=[(i,i+d) for i in range(M+1-d)]
    before=rank(list(FIXED.values())+Ed)-len(FIXED)
    after=ranks[d]
    reconnect=rank(fixed+Ed+[FIXED[release]])==rank(fixed+Ed)
    check(after-before==int(reconnect),'exact one-edge rank gain')
  row={'release':release,'states':len(states),'components':sorted(map(len,components)),
       'ranks':ranks,'staircase':all(q>=len(colors)-i for i,q in enumerate(sorted(ranks.values(),reverse=True)))}
  if release==1:
   path=[((5,6),(2,8),(0,9),(0,11)),((1,2),(0,6),(0,9),(0,11)),
         ((1,2),(0,6),(2,11),(1,12)),((5,6),(4,10),(2,11),(1,12))]
   for s in path:check(s in states,'explicit repair state')
   for a,b in zip(path,path[1:]):check(sum(x!=y for x,y in zip(a,b))==2,'repair width')
   src,dst=map(states.index,[path[0],path[-1]])
   distances={src:0};todo=deque([src])
   while todo:
    v=todo.popleft()
    for w in adj[v]:
     if w not in distances:distances[w]=distances[v]+1;todo.append(w)
   check(distances[dst]==3,'shortest three-move repair')
   row['shortest_repair_length']=3;row['explicit_path']=path
  rows.append(row)
 expected=[5,54,56,24,20,54,14,20,16,5]
 check([r['states'] for r in rows]==expected,'release state counts')
 check([r['release'] for r in rows if len(r['components'])==1]==[1,2,3,4,5,7,8,10],'exact successful releases')
 check([r['release'] for r in rows if r['staircase']]==[1,2,4,5,7,8],'rank-certified releases')
 return rows

def terras(x):return (3*x+1)//2 if x%2 else x//2
def branch(word):
 a,h,m=1,0,1
 for b in word:
  if b=='1':a*=3;h=3*h+m
  m*=2
 return a,h,m

def avoidance(word,nmax):
 L=len(word);dp=[0]*L;dp[0]=1;counts=[1]
 transitions={}
 for s in range(L):
  for b in '01':
   t=word[:s]+b
   if t.endswith(word):transitions[s,b]=None;continue
   transitions[s,b]=max(j for j in range(L) if t.endswith(word[:j]))
 for n in range(nmax):
  ndp=[0]*L
  for s,v in enumerate(dp):
   for b in '01':
    t=transitions[s,b]
    if t is not None:ndp[t]+=v
  dp=ndp;counts.append(sum(dp))
 return counts

def collatz_controls():
 for N in range(1,15):
  seen=set()
  for x in range(2**N):
   y=x;bits=0
   for i in range(N):bits|=(y%2)<<i;y=terras(y)
   seen.add(bits)
  check(len(seen)==2**N,'finite parity encoding bijection')
 K='110110110110110101011011101010';W=K+'10';a,h,m=branch(W)
 rho=(-h*pow(a,-1,m))%m
 check((a,h,m,rho)==(3486784401,15810005461,4294967296,2962948091),'explicit guarded KF branch')
 check(h<20*(m-a) and rho>20,'all-lift descent threshold')
 for lift in range(256):
  x=rho+lift*m;y=x;got=[]
  for _ in W:got.append(str(y%2));y=terras(y)
  check(''.join(got)==W and y==(a*x+h)//m and y<x,'actual guarded descent')
 check((a*rho+h)//m==2405410912,'least endpoint')
 borders=[j for j in range(1,len(W)+1) if W[:j]==W[-j:]]
 check(borders==[32],'unbordered KF')
 exact=avoidance(W,512);rec=[]
 for n in range(513):
  r=2**n if n<32 else 2*rec[n-1]-rec[n-32]
  rec.append(r)
  incl=sum((-1)**j*comb(n-31*j,j)*2**(n-32*j) for j in range(n//32+1))
  check(exact[n]==r==incl,'automaton recurrence inclusion-exclusion agreement')
 for q in range(1,257):
  x=2**q-1;y=x
  for j in range(q):
   check(y%2==1,'all-odd initial parity');y=terras(y)
   check(y==3**(j+1)*2**(q-j-1)-1 and y>x,'exact increasing Mersenne run')
 # Small-pattern exact disjoint block control; illustrates the general count without a 2^32 enumeration.
 small='101';L=len(small);k=4;count=0
 for x in range(2**(L*k)):
  y=x;bits=''
  for _ in range(L*k):bits+=str(y%2);y=terras(y)
  count+=all(bits[j*L:(j+1)*L]!=small for j in range(k))
 check(count==(2**L-1)**k,'disjoint cylinder block count')
 return {'word':W,'affine':[a,h,m],'residue':rho,'least_endpoint':(a*rho+h)//m,
         'borders':borders,'all_lifts_checked':256,'avoidance_lengths_checked':513,
         'avoidance_counts':{str(n):exact[n] for n in [32,33,63,64,96,128,256,512]},
         'parity_bijection_checked_through':14,'Mersenne_runs_checked_through':256}

def main():
 b=berggren_controls();g=release_controls();c=collatz_controls()
 print(json.dumps({'status':'PASS','checks':CHECKS,'berggren':b,'graceful':g,'collatz':c},indent=2,sort_keys=True))
if __name__=='__main__':main()
