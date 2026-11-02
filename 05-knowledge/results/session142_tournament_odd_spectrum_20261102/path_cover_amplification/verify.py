#!/usr/bin/env python3
"""Exact finite controls for path-cover amplification; standard library only."""
import argparse,json,math
from functools import lru_cache
from itertools import permutations

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def tournament(n,code):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j) if code>>k&1 else (j,i);a[u]|=1<<v;k+=1
 return tuple(a)
def hp(a):
 n=len(a);N=1<<n;D=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    D[mask][v]=1 if not rest else sum(D[rest][u] for u in range(n) if rest>>u&1 and a[u]>>v&1)
  H[mask]=sum(D[mask])
 return H
def covers(H,n):
 @lru_cache(None)
 def rec(mask):
  if not mask:return (1,)+(0,)*n
  first=mask&-mask;rest=mask^first;result=[0]*(n+1);s=rest
  while True:
   b=s|first;tail=rec(mask^b)
   for k in range(1,n+1):result[k]+=H[b]*tail[k-1]
   if not s:break
   s=(s-1)&rest
  return tuple(result)
 return rec((1<<n)-1)
def subst(q,port,a):
 # Labels: module first, then quotient exterior in increasing old label.
 n=len(a);outside=[v for v in range(len(q)) if v!=port];b=list(a)+[0]*len(outside)
 for j,v in enumerate(outside):
  z=n+j
  for i in range(n):
   if q[port]>>v&1:b[i]|=1<<z
   else:b[z]|=1<<i
  for k,w in enumerate(outside):
   if q[v]>>w&1:b[z]|=1<<(n+k)
 return tuple(b)
def word_coeff(q,port,k):
 m=len(q);outside=tuple(v for v in range(m) if v!=port)
 @lru_cache(None)
 def rec(mask,remaining,last):
  if not mask and not remaining:return 1
  total=0
  if remaining and last!=port and (last==m or q[last]>>port&1):total+=rec(mask,remaining-1,port)
  for v in outside:
   if mask>>v&1 and (last==m or q[last]>>v&1):total+=rec(mask^(1<<v),remaining,v)
  return total
 return rec(sum(1<<v for v in outside),k,m)
def cyclic(q,p):
 # Reachability, deliberately separate from the triangle test.
 reached=0;front=q[p]
 while front:
  b=front&-front;front-=b;v=b.bit_length()-1
  if reached&b:continue
  reached|=b;front|=q[v]&~reached
 return bool(reached>>p&1)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-order',type=int,default=5);args=ap.parse_args();N=args.max_order
 need(1<=N<=5,'max-order must be1..5')
 rows=[];equalities=[];checks=0;literal_checks=0
 c3=(2,4,1)
 for n in range(1,N+1):
  count=0;mineven=None
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code);H=hp(a);h=H[-1];p=covers(H,n)
   need(p[1]==h and p[n]==1 and h%2==1,'basic profile/parity')
   if n<=4:
    literal=0;direct=[set() for _ in range(n+1)];fibres=[{} for _ in range(n+1)]
    for w in permutations(range(n)):
     fullpath=all(a[x]>>y&1 for x,y in zip(w,w[1:]));literal+=fullpath
     for cuts in range(1<<(n-1)):
      blocks=[];start=0
      for j in range(n-1):
       if cuts>>j&1:blocks.append(w[start:j+1]);start=j+1
      blocks.append(w[start:]);k=len(blocks)
      if all(all(a[x]>>y&1 for x,y in zip(b,b[1:])) for b in blocks):
       cover=tuple(sorted(blocks));direct[k].add(cover)
       if fullpath:fibres[k][cover]=fibres[k].get(cover,0)+1
    need(literal==h,'literal/DP mismatch')
    for k in range(1,n+1):
     need(len(direct[k])==p[k],'literal unordered-cover mismatch')
     need(sum(fibres[k].values())==h*math.comb(n-1,k-1),'literal incidence total')
     need(all(z<=math.factorial(k) for z in fibres[k].values()),'literal fibre bound')
    literal_checks+=1
   for k in range(1,n+1):
    need(math.factorial(k)*p[k]>=h*math.comb(n-1,k-1),'incidence bound');checks+=1
   p2=p[2] if n>=2 else 0
   actual=hp(subst(c3,0,a))[-1]
   need(actual==3*h+2*p2,'C3 kernel')
   need(actual>=(n+2)*h+(1 if n%2==0 else 0),'amplification')
   if actual==(n+2)*h+(1 if n%2==0 else 0) and len(equalities)<12:
    equalities.append({'n':n,'code':code,'out_masks':a,'H':h,'p2_unordered':p2,'C3_H':actual})
   count+=1
  rows.append({'order':n,'labelled_tournaments':count})
 context_count=active=inactive=kernelchecks=0
 modules=[tournament(1,0),tournament(2,0),tournament(3,0),tournament(3,2)]
 # code2 is cyclic under this lexicographic convention; verify rather than assume.
 need(hp(modules[-1])[-1]==3,'cyclic control seed')
 for m in range(1,5):
  for code in range(1<<(m*(m-1)//2)):
   q=tournament(m,code);hq=hp(q)[-1]
   for port in range(m):
    W=[0]+[word_coeff(q,port,k) for k in range(1,5)]
    need(W[1]==hq,'single-run coefficient')
    act=cyclic(q,port)
    triangle=any(q[port]>>u&1 and q[u]>>v&1 and q[v]>>port&1 for u in range(m) for v in range(m))
    need(act==triangle,'cycle/triangle equivalence')
    if act:need(W[2]>=1 and hq>=3,'active port coefficients');active+=1
    else:
     need(all(x==0 for x in W[2:]),'inactive return')
     need(hq==hp(q)[((1<<m)-1)^(1<<port)],'inactive deleted-port factor')
     inactive+=1
    for a in modules:
     n=len(a);H=hp(a);p=covers(H,n);actual=hp(subst(q,port,a))[-1]
     predicted=sum(math.factorial(k)*W[k]*p[k] for k in range(1,n+1))
     need(actual==predicted,'whole kernel mismatch')
     lower=H[-1]*sum(W[k]*math.comb(n-1,k-1) for k in range(1,n+1))
     need(actual>=lower,'binomial kernel bound');kernelchecks+=1
     if not act:need(actual==hq*H[-1],'inactive multiplicativity')
    context_count+=1
 print(json.dumps({'status':'PASS','scope':'Finite controls; universal theorem proved separately; unrestricted spectrum OPEN',
  'max_order':N,'labelled_universes':rows,'incidence_inequalities_checked':checks,'literal_permutation_and_unordered_cover_checks':literal_checks,
  'pointed_contexts_orders1to4':context_count,'active_contexts':active,'inactive_contexts':inactive,
  'literal_context_module_kernel_checks':kernelchecks,'sharp_controls':equalities,
  'hostile_scope':'S0/R1 new vertices are nonhomogeneous to A+x, so this theorem does not apply to those operations'},sort_keys=True,indent=2))
if __name__=='__main__':main()
