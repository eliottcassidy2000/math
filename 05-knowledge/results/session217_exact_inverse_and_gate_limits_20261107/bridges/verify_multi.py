#!/usr/bin/env python3
import itertools,json
# Import definitions only; no imported verification result is used.
s=open(__file__.replace('verify_multi.py','verify.py')).read();exec(s[:s.index('if hp([40')])
def ordinal_pad(r,m):
 n=len(r)
 return [x|(((1<<m)-1)<<n) for x in r]+[sum(1<<j for j in range(n+i+1,n+m)) for i in range(m)]
def substitute(q,blocks):
 labels=[(i,v) for i,b in enumerate(blocks) for v in range(len(b))]
 return [sum(1<<b for b,(j,w) in enumerate(labels) if (i==j and blocks[i][v]>>w&1) or(i!=j and q[i]>>j&1)) for i,v in labels]
count=safe_count=unsafe_count=0
seeds=([0],[2,4,1])
for n in range(2,5):
 masks=range(1<<(n*(n-1)//2)) if n<=3 else [i*4 for i in range(16)]
 for mask in masks:
  q=tournament(n,mask);hq=hp(q)
  for holes in itertools.combinations(range(n),2):
   safe=all(single_scc(q,i) for i in holes)
   for inputs in itertools.product(seeds,repeat=2):
    product=hp(inputs[0])*hp(inputs[1])
    for m in range(2):
     blocks=[[0] for _ in q]
     for i,u in zip(holes,inputs):blocks[i]=ordinal_pad(u,m)
     actual=hp(substitute(q,blocks));count+=1
     if safe:
      safe_count+=1
      if actual!=hq*product:raise RuntimeError('safe multi-hole factorization')
     else:
      unsafe_count+=1
      # Includes the disjoint one-run baseline and an unsafe two-run term.
      if actual<(hq+2**(m+1)-2)*product:raise RuntimeError('uniform unsafe lower bound')
print(json.dumps(dict(status='PASS',literal_composites=count,safe=safe_count,unsafe=unsafe_count,universe='all order2..3 quotients and order4 masks 0,4,..60; all two-hole choices; singleton/C3 input pairs; simultaneous ordinal padding m=0,1',checks=['independent literal substitution followed by subset DP','safe output H(Q) product h','unsafe uniform lower bound including one-run baseline']),indent=2))
