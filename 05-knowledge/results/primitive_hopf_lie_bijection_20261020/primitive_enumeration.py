#!/usr/bin/env python3
"""Integer Burnside buckets + ordinary and signed primitive transforms.
Stdlib only. No graph atlas needed. Independent Lie-character cross-check.
"""
from __future__ import annotations
import argparse,csv,math,json,itertools
from collections import Counter,defaultdict

CHECKS=0
def check(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise AssertionError(msg)

def odd_partitions(n,least=1):
 if not n:yield ();return
 for k in range(least,n+1,2):
  for tail in odd_partitions(n-k,k):yield (k,)+tail

def zvalue(parts):
 z=1
 for d,m in Counter(parts).items():z*=d**m*math.factorial(m)
 return z

def pair_orbits(parts):
 return sum(d//2 for d in parts)+sum(math.gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])

def burnside(n,converse=False):
 """Scale by n!, bucket equal exponents, then perform one exact division."""
 fact=math.factorial(n);buckets=defaultdict(int)
 source=odd_partitions(n//2 if converse else n)
 for p in source:
  parts=tuple(2*d for d in p)+((1,) if n%2 else ()) if converse else p
  z=zvalue(parts);check(fact%z==0,'centralizer divides factorial')
  buckets[pair_orbits(parts)]+=fact//z
 numerator=sum(v<<k for k,v in buckets.items())
 check(numerator%fact==0,'Burnside integer')
 return numerator//fact

def ghost(a):
 N=len(a)-1;g=[0]*(N+1)
 for n in range(1,N+1):g[n]=n*a[n]-sum(g[k]*a[n-k] for k in range(1,n))
 return g

def inv_euler(a):
 g=ghost(a);N=len(a)-1;c=[0]*(N+1)
 for n in range(1,N+1):
  v=g[n]-sum(d*c[d] for d in range(1,n) if n%d==0)
  check(v%n==0,'Euler integer');c[n]=v//n
 return c

def mobius(n):
 r=1;p=2
 while p*p<=n:
  if n%p==0:
   n//=p;r=-r
   if n%p==0:return 0
   while n%p==0:n//=p
  p+=1
 return -r if n>1 else r

def powers(a):
 N=len(a)-1;ps=[[1]+[0]*N]
 for k in range(1,N+1):
  old=ps[-1];ps.append([sum(old[i]*a[n-i] for i in range(n+1)) for n in range(N+1)])
 return ps

def reconstruct_product(plus,minus):
 N=len(plus)-1;a=[1]+[0]*N
 # Direct finite Euler multiplication, independent from ghost inversion.
 for k in range(1,N+1):
  for count,sgn in ((plus[k],1),(minus[k],-1)):
   if not count:continue
   fac=[0]*(N+1);fac[0]=1;v=1
   for j in range(1,N//k+1):
    v=v*(count+j-1)//j;fac[j*k]=sgn**j*v
   a=[sum(a[i]*fac[n-i] for i in range(n+1)) for n in range(N+1)]
 return a

def oriented_converse_control(N=6):
 """Independent all-permutation Burnside for three pair states: absent, forward, back."""
 totals=[];fixed=[]
 for n in range(N+1):
  edges=list(itertools.combinations(range(n),2));ix={e:i for i,e in enumerate(edges)}
  ordinary_sum=converse_sum=0
  for p in itertools.permutations(range(n)):
   dest=[];flips=[]
   for u,v in edges:
    a,b=p[u],p[v];dest.append(ix[tuple(sorted((a,b)))]);flips.append(a>b)
   seen=set();ordinary_good=converse_good=0
   for j in range(len(edges)):
    if j in seen:continue
    k=j;length=parity=0
    while k not in seen:
     seen.add(k);length+=1;parity^=flips[k];k=dest[k]
    ordinary_good+=parity==0
    converse_good+=(parity^(length%2))==0
   ordinary_sum+=3**ordinary_good;converse_sum+=3**converse_good
  fact=math.factorial(n)
  check(ordinary_sum%fact==converse_sum%fact==0,'oriented Burnside integrality')
  totals.append(ordinary_sum//fact);fixed.append(converse_sum//fact)
 connected=inv_euler(totals);fake=inv_euler(fixed);true=[0]*(N+1)
 for n in range(1,N+1):
  true[n]=fake[n]-((connected[n//2]-true[n//2])//2 if n%2==0 else 0)
  check(true[n]>=0 and (connected[n]-true[n])%2==0,'oriented converse pairs')
 check(connected[3]==5 and true[3]==3,'five connected oriented triples; one converse pair')
 check(fake[6]==727 and true[6]==726,'paired-star historical correction independently recomputed')
 return totals,fixed,connected,true,fake

def run(N):
 T=[burnside(n) for n in range(N+1)];R=[burnside(n,True) for n in range(N+1)]
 c=inv_euler(T);h=ghost(R);d=[0]*(N+1);s=[0]*(N+1);q=[0]*(N+1)
 for n in range(1,N+1):
  s[n]=T[n]-sum(s[k]*T[n-k] for k in range(1,n))
  q[n]=R[n]-sum(s[j]*R[n-2*j] for j in range(1,n//2+1))
  v=h[n]-sum(k*(d[k] if (n//k)%2 else c[k]) for k in range(1,n) if n%k==0)
  check(v%n==0,'signed primitive integer');d[n]=v//n
 plus=[(x+y)//2 for x,y in zip(c,d)];minus=[(x-y)//2 for x,y in zip(c,d)]
 for n in range(N+1):
  check((c[n]+d[n])%2==0 and plus[n]>=0 and minus[n]>=0,'parity eigendimensions')
  check(0<=q[n]<=s[n] and (s[n]-q[n])%2==0,'strong converse orbits')
 check(reconstruct_product(plus,minus)==R,'signed Euler product reconstruction')
 check(reconstruct_product(c,[0]*(N+1))==T,'ordinary Euler product reconstruction')
 ps=powers(s);pq=powers(q);lengthrows=[]
 for n in range(1,N+1):
  dimsum=tracesum=0
  for m in range(1,n+1):
   dimnum=tracenum=0
   for k in range(1,math.gcd(n,m)+1):
    if n%k or m%k:continue
    mu=mobius(k)
    dimnum+=mu*ps[m//k][n//k]
    tracenum+=mu*(pq if k%2 else ps)[m//k][n//k]
   check(dimnum%m==tracenum%m==0,'Witt character integer')
   dim=dimnum//m;theta=tracenum//m;eta=(-1)**(m-1)*theta
   check(dim>=abs(eta) and (dim+eta)%2==0,'length eigendimensions')
   dimsum+=dim;tracesum+=eta
   if dim:lengthrows.append({'n':n,'length':m,'dimension':dim,'converse_trace':eta})
  check(dimsum==c[n] and tracesum==d[n],'Lie and graph transforms agree')
 print('n tournaments strong self_converse self_converse_strong connected_sign_even even_edges odd_edges primitive_trace')
 rows=[]
 for n in range(N+1):
  r=[n,T[n],s[n],R[n],q[n],c[n],plus[n],minus[n],d[n]];rows.append(r)
  if n<=16:print(*r)
 # The old InvEuler(self-converse oriented graphs) table contains paired atoms.
 ot,of,oc,oq,oe=oriented_converse_control()
 print('Independent oriented-graph all-permutation Burnside to n=6:')
 print('all=',ot,'self_converse=',of)
 print('connected=',oc,'connected_self_converse=',oq,'InvEuler_self_converse=',oe)
 print('Historical InvEuler(A005639) at n=6 is 727; connected self-converse oriented = 726 (one paired-star atom removed).')
 print('N=',N,'ALL CHECKS PASSED:',CHECKS)
 return rows,lengthrows

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-n',type=int,default=60);p.add_argument('--prefix',default='primitive_enumeration');a=p.parse_args()
 if not 0<=a.max_n<=100:raise SystemExit('Require 0 <= max_n <= 100; default resource boundary 60.')
 rows,lengthrows=run(a.max_n)
 with open(a.prefix+'.csv','w') as f:
  w=csv.writer(f);w.writerow(['n','tournaments','strong','self_converse','self_converse_strong','connected_sign_even','even_edges','odd_edges','primitive_trace']);w.writerows(rows)
 with open(a.prefix+'_length.json','w') as f:json.dump(lengthrows,f,indent=2);f.write('\n')
