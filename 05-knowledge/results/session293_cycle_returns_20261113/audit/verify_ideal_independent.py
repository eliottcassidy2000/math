#!/usr/bin/env python3
"""Independent uniform-module ideal checks, explicit assertions active under -O."""
import json,math
from functools import lru_cache
from itertools import permutations,combinations
checks=0

def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def tour(n,mask):
 r=[0]*n
 for j,(a,b) in enumerate(combinations(range(n),2)):
  if mask>>j&1:r[a]|=1<<b
  else:r[b]|=1<<a
 return tuple(r)
def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
@lru_cache(None)
def h(r,S):
 if not S:return 1
 dp={(1<<i,i):1 for i in range(len(r)) if S>>i&1}
 for _ in range(1,S.bit_count()):
  nxt={}
  for (used,last),num in dp.items():
   for j in range(len(r)):
    if (S&~used&r[last])>>j&1:
     key=(used|(1<<j),j);nxt[key]=nxt.get(key,0)+num
  dp=nxt
 return sum(dp.values())
@lru_cache(None)
def fcover(r,k,S):
 if k==0:return int(S==0)
 if S.bit_count()<k:return 0
 ans=0;U=S
 while U:
  ans+=h(r,U)*fcover(r,k-1,S^U);U=(U-1)&S
 return ans

def coefficients(ext,out):
 # out = exterior vertices beaten by C; complement beats C.
 n=len(ext);a=[0]*(n+2)
 for perm in permutations(range(n)):
  gaps=[g for g in range(n+1) if (g==0 or perm[g-1] not in out) and (g==n or perm[g] in out)]
  for occ in range(1<<len(gaps)):
   chosen={g for j,g in enumerate(gaps) if occ>>j&1}
   valid=all(g in chosen or bool(ext[perm[g-1]]>>perm[g]&1) for g in range(1,n))
   if valid:a[len(chosen)]+=1
 return a

def insertion(ext,out,C):
 n=len(ext);r=list(ext)
 for j in range(n):
  if j not in out:r[j]|=((1<<len(C))-1)<<n
 r.extend((row<<n)|sum(1<<j for j in out) for row in C)
 return tuple(r)

probe=[]
for m in range(5):
 a=coefficients(tt(2*m),set(range(0,2*m,2)))
 for k in range(1,len(a)):
  need(a[k]%math.factorial(k-1)==0,'probe excursion divisibility')
 need(a[m+1]==math.factorial(m),'probe leading coefficient')
 need(not any(a[m+2:]),'probe degree')
 b=[a[k]//math.factorial(k-1) for k in range(1,m+2)]
 probe.append({'m':m,'a':a,'b':b})

contexts=0
for n in range(0,5):
 for bits in range(1<<(n*(n-1)//2)):
  ext=tour(n,bits)
  for color in range(1<<n):
   out={i for i in range(n) if color>>i&1}
   a=coefficients(ext,out)
   for k in range(1,len(a)):need(a[k]%math.factorial(k-1)==0,'all small context excursion divisibility')
   contexts+=1

examples=[]
for mask in (8,76):
 C=tour(5,mask);S=31
 F=[fcover(C,k,S) for k in range(1,6)]
 D=[math.factorial(k-1)*F[k-1] for k in range(1,6)]
 g=math.gcd(*D)
 aut=0
 for p in permutations(range(5)):
  aut+=all(bool(C[i]>>j&1)==bool(C[p[i]]>>p[j]&1) for i in range(5) for j in range(5))
 expect=(9,1,[9,66,186,240,120],3) if mask==8 else (15,5,[15,90,210,240,120],15)
 need((h(C,S),aut,F,g)==expect,'published module example')
 observed=[]
 for m in range(5):
  T=insertion(tt(2*m),set(range(0,2*m,2)),C)
  observed.append(h(T,(1<<len(T))-1))
  expected=sum(probe[m]['a'][k]*F[k-1] for k in range(1,min(5,len(probe[m]['a'])-1)+1))
  need(observed[-1]==expected,'direct physical triangular probe')
 need(math.gcd(*observed)==g,'probe gcd matches scaled cover ideal')
 # Explicit unitriangular inversion, no rational solve or assumed integer inverse.
 recovered=[]
 for m,pval in enumerate(observed):
  recovered.append(pval-sum(probe[m]['b'][j]*recovered[j] for j in range(m)))
 need(recovered==D,'integral triangular recovery')
 examples.append({'mask':mask,'H':h(C,S),'Aut':aut,'F':F,'D':D,'g':g,'probes':observed})
 # All three-vertex contexts, not only the alternating probes.
 for bits in range(8):
  ext=tour(3,bits)
  for color in range(8):
   out={i for i in range(3) if color>>i&1};a=coefficients(ext,out)
   T=insertion(ext,out,C);value=h(T,(1<<len(T))-1)
   need(value==sum(a[k]*F[k-1] for k in range(1,min(5,len(a)-1)+1)),'all n3 context direct cover identity')
   need(value%g==0,'all n3 context divisor')
# Explicit exploratory witness where the excursion scaling adds a divisor.
C=tour(8,180752028);full=255
F=[fcover(C,k,full) for k in range(1,9)]
D=[math.factorial(k-1)*F[k-1] for k in range(1,9)]
need(C==(56,117,129,22,100,76,141,59),'order8 witness decoding')
need(F==[513,6102,32508,99432,187200,216000,141120,40320],'order8 cover fingerprint')
need(math.gcd(*F)==9 and math.gcd(*D)==27,'strict excursion-factorial strengthening')
# No minimal-order assertion: this was found by exploratory sampling.
automorphisms8=sum(all(bool(C[i]>>j&1)==bool(C[p[i]]>>p[j]&1)
                       for i in range(8) for j in range(8))
                   for p in permutations(range(8)))
need(automorphisms8==1,'asymmetric order8 witness')
extra={'Aut':automorphisms8,'n':8,'mask':180752028,'rows':list(C),'F':F,'D':D,'ordered_cover_gcd':math.gcd(*F),'uniform_context_ideal':math.gcd(*D)}
print(json.dumps({'strict_strengthening':extra,'status':'PASS','checks':checks,'small_contexts':contexts,'probe_coefficients':probe,'examples':examples},indent=2,sort_keys=True))
