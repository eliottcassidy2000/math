#!/usr/bin/env python3
"""Exact Legendre quotient and Frobenius tests, standard library only."""
from math import comb
import json
checks=0

def ck(b,m):
 global checks
 checks+=1
 if not b:raise RuntimeError(m)

def trim(a):
 while a and a[-1]==0:a.pop()
 return a

def add(a,b,p):
 return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])
def neg(a,p):return [(-v)%p for v in a]
def mul(a,b,p):
 c=[0]*(max(0,len(a)+len(b)-1))
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
 return trim(c)
def divrem(a,b,p):
 a=a[:];q=[0]*max(0,len(a)-len(b)+1);ib=pow(b[-1],-1,p)
 while len(a)>=len(b):
  i=len(a)-len(b);v=a[-1]*ib%p;q[i]=v
  for j,x in enumerate(b):a[i+j]=(a[i+j]-v*x)%p
  trim(a)
 return trim(q),a
def mod(a,h,p):return divrem(a,h,p)[1]
def pw(a,n,h,p):
 ans=[1]
 while n:
  if n&1:ans=mod(mul(ans,a,p),h,p)
  a=mod(mul(a,a,p),h,p);n//=2
 return ans
def gcd(a,b,p):
 while b:a,b=b,mod(a,b,p)
 return [(x*pow(a[-1],-1,p))%p for x in a]
def inv(a,h,p):
 r0,r1=h,a;s0,s1=[],[1]
 while r1:
  q,r=divrem(r0,r1,p);r0,r1=r1,r;s0,s1=s1,add(s0,neg(mul(q,s1,p),p),p)
 ck(len(r0)==1,'unit denominator')
 return mod([(v*pow(r0[0],-1,p))%p for v in s0],h,p)
def evaluate(a,b,h,p):
 out=[]
 for v in reversed(a):out=add(mod(mul(out,b,p),h,p),[v],p)
 return out

def minpoly(j,h,p):
 n=len(h)-1;basis={};v=[1]
 for k in range(n+1):
  vec=v+[0]*(n-len(v));tag=[0]*k+[1]
  for pivot,(old,t) in sorted(basis.items()):
   c=vec[pivot]
   if c:
    vec=[(x-c*y)%p for x,y in zip(vec,old)]
    tag=add(tag,[(-c*y)%p for y in t],p)
  if not any(vec):
   ck(tag[-1]==1,'monic minimal relation')
   ck(not evaluate(tag,j,h,p),'minimal relation verified')
   return tag
  pivot=next(i for i,x in enumerate(vec) if x);c=pow(vec[pivot],-1,p)
  basis[pivot]=([(x*c)%p for x in vec],[(x*c)%p for x in tag])
  v=mod(mul(v,j,p),h,p)
 raise RuntimeError('no relation')

def primes(B):return [n for n in range(5,B+1) if all(n%d for d in range(2,int(n**.5)+1))]
rows=[];ogg=[]
for p in primes(199):
 n=(p-1)//2;h=[comb(n,i)**2%p for i in range(n+1)]
 derivative=[i*h[i]%p for i in range(1,len(h))]
 ck(len(gcd(h,derivative,p))==1,('squarefree H',p))
 den=[0,0,1,-2%p,1];num=[256*x%p for x in mul(mul([1,p-1,1],[1,p-1,1],p),[1,p-1,1],p)]
 j=mod(mul(num,inv(den,h,p),p),h,p)
 ck(not add(pw([0,1],p*p,h,p),[0,p-1],p),('Legendre Fp2',p))
 # Exact S3 substitutions preserve H and its j-observation.
 for g in ([1,p-1],inv([0,1],h,p)):
  ck(not evaluate(h,g,h,p),('H stable under level permutation',p))
  ck(evaluate(j,g,h,p)==j,('j invariant under level permutation',p))
 sj=minpoly(j,h,p);X=[0,1]
 fixed=len(gcd(sj,add(pw(X,p,sj,p),[0,p-1],p),p))-1
 allfixed=not add(pw(j,p,h,p),neg(j,p),p)
 ck(allfixed==(fixed==len(sj)-1),('two quotient Frobenius tests',p))
 expected=p//12+{1:0,5:1,7:1,11:2}[p%12]
 ck(len(sj)-1==expected,('classical degree',p))
 ck((len(sj)-1-fixed)%2==0,('quadratic pair count',p))
 lf=len(gcd(h,add(pw(X,p,h,p),[0,p-1],p),p))-1
 if allfixed:ogg.append(p)
 rows.append({'p':p,'j_polynomial':sj,'j_count':len(sj)-1,'rational_j':fixed,'nonrational_j_pairs':(len(sj)-1-fixed)//2,'rational_level_parameters':lf,'all_j_rational':allfixed})
ck(ogg==[5,7,11,13,17,19,23,29,31,41,47,59,71],'Ogg primes within finite scope')
r5=next(r for r in rows if r['p']==5);ck(r5['rational_j']==1 and r5['rational_level_parameters']==0,'coarse-fixed not mark-fixed')
r37=next(r for r in rows if r['p']==37);ck(r37['nonrational_j_pairs']==1,'first non-Ogg pair')
p=71;n=35;h=[comb(n,i)**2%p for i in range(n+1)];fibres={}
for x in range(2,p):
 if sum(c*pow(x,i,p) for i,c in enumerate(h))%p==0:
  v=256*pow(1-x+x*x,3,p)*pow(x*x*(1-x)**2,-1,p)%p
  fibres.setdefault(v,[]).append(x)
ck({j:len(v) for j,v in fibres.items()}=={17:6,24:3,40:6,48:6},'71 named fibres')
ck(mul([29,1],[31,31,1],37)==r37['j_polynomial'],'37 exact factor and sign')
print(json.dumps({'status':'PASS' ,'checks':checks,'prime_bound':199,'odd_primes_at_least5':len(rows),'ogg_in_scope':ogg,'rows':rows},sort_keys=True))
