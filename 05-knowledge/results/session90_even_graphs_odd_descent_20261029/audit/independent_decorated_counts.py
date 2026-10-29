#!/usr/bin/env python3
"""Independent exhaustive audit of signed decorated graph/pair counts.
No imports from producer scripts. Checks use explicit exceptions, including -O.
"""
from itertools import permutations
from collections import Counter
import json

CHECKS=0
def check(p,msg):
 global CHECKS
 CHECKS+=1
 if not p: raise RuntimeError(msg)

def setup(n):
 es=[(i,j) for i in range(n) for j in range(i+1,n)]; ei={e:i for i,e in enumerate(es)}; m=len(es)
 acts=[]
 for p in permutations(range(n)):
  dest=[ei[tuple(sorted((p[i],p[j])))] for i,j in es]
  c=sum(1<<dest[k] for k,(i,j) in enumerate(es) if p[i]>p[j])
  ep=[sum(((x>>k)&1)<<dest[k] for k in range(m)) for x in range(1<<m)]
  vp=[sum(((d>>i)&1)<<p[i] for i in range(n)) for d in range(1<<n)]
  rev=sum(1<<k for k,(i,j) in enumerate(es) if p[i]>p[j])
  acts.append((p,ep,vp,c,rev))
 return es,acts

def audit(n):
 es,acts=setup(n);m=len(es);q=1<<m;N=1<<(2*m+n)
 def parts(z): return z&(q-1),(z>>m)&(q-1),z>>(2*m)
 def pack(x,r,d):return x|(r<<m)|(d<<(2*m))
 seenP=set();seenG=set(); pc=gc=good=base_even_good=0; p_hist=Counter();g_hist=Counter(); witness=None
 for z in range(N):
  x,r,d=parts(z)
  if z not in seenP:
   orb={pack(ep[x]^c,ep[r],vp[d]) for _,ep,vp,c,rev in acts}
   seenP.update(orb);pc+=1;p_hist[len(acts)//len(orb)]+=1
  if z not in seenG:
   orb={pack(ep[x],ep[r],vp[d]) for _,ep,vp,c,rev in acts}
   seenG.update(orb);gc+=1
   st=[a for a in acts if a[1][x]==x and a[1][r]==r and a[2][d]==d]
   even=all((x&a[4]).bit_count()%2==0 for a in st)
   if even:
    good+=1;g_hist[len(st)]+=1
    base_even=all((x&a[4]).bit_count()%2==0 for a in acts if a[1][x]==x)
    base_even_good+=base_even
    if not base_even and witness is None:witness=[x,r,d]
   check(len(orb)*len(st)==len(acts),'orbit stabilizer')
 # Per-permutation independent signed fixed traces versus tournament traces.
 for p,ep,vp,c,rev in acts:
  signed=sum((-1)**((x&rev).bit_count()%2) for x in range(q) if ep[x]==x)
  fixed_t=sum(ep[x]^c==x for x in range(q))
  check(signed==fixed_t,('signed trace',n,p))
  decor=sum(ep[r]==r for r in range(q))*sum(vp[d]==d for d in range(1<<n))
  check(signed*decor==fixed_t*decor,'decorated signed trace')
 check(pc==good,('class count',n,pc,good))
 check(len(seenP)==len(seenG)==N,'complete universes')
 return dict(n=n,labelled=N,pairs=pc,all_graph_triples=gc,relative_even=good,base_even_only=base_even_good,pair_stabilizers=dict(p_hist),graph_stabilizers=dict(g_hist),relative_not_base_witness=witness)

def main():
 rows=[audit(n) for n in range(5)]
 check([r['pairs'] for r in rows]==[1,2,8,88,2752],'known counts')
 print(json.dumps(dict(status='PASS',checks=CHECKS,rows=rows),indent=2))
if __name__=='__main__':main()
