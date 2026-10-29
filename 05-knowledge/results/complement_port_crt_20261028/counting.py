"""Classical Burnside/unit-tail kernel with ordered shared-port weights."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from math import gcd,factorial,comb
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)
def partitions(n,minimum=1):
 if n==0:yield ();return
 for a in range(minimum,n+1,2):
  for tail in partitions(n-a,a):yield (a,)+tail
def weight(parts):
 e=sum((a-1)//2 for a in parts)+sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
 z=1
 for a,m in Counter(parts).items():z*=a**m*factorial(m)
 return Q(2**e,z)
def direct(n,R):
 out=[Q(0)]*(R+1)
 for part in partitions(n):
  w=weight(part);j=part.count(1)
  for r in range(R+1):out[r]+=w*j**r
 return out
def fall(n,r):
 return factorial(n)//factorial(n-r) if r<=n else 0
def main(write=True,N=60,R=8):
 global CHECKS
 CHECKS=0;core={};rows=[]
 reference=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168,
            48542114686912,28401423719122304,31021002160355166848,
            63530415842308265100288,244912778438520759443245824,
            1783398846284777975419600287232,24605641171260376770598003978281472]
 core_terms=0
 for m in range(N+1):
  for part in partitions(m,3):
   key=(m,len(part));core[key]=core.get(key,Q(0))+weight(part);core_terms+=1
 for n in range(N+1):
  coeff=defaultdict(Q)
  for (m,t),w in core.items():
   if m<=n:
    j=n-m;coeff[j]+=w*Q(2**(j*(j-1)//2+j*t),factorial(j))
  counts=[sum(w*j**r for j,w in coeff.items()) for r in range(R+1)]
  need(all(c.denominator==1 for c in counts),'integer ordered-port orbit counts')
  if n<=26:need(counts==direct(n,R),'independent direct odd-partition formula')
  if n<len(reference):need(counts[0]==reference[n],'published A000568 prefix')
  distinct=[sum(w*fall(j,r) for j,w in coeff.items()) for r in range(n+1)]
  need(all(c.denominator==1 for c in distinct),'integer distinct-port counts')
  need(distinct[n]==2**(n*(n-1)//2),'all named ports recover labelled tournaments')
  if n:need(distinct[n-1]==distinct[n],'all but one port fixes the frame')
  # Stirling conversion is a separate recurrence over collision partitions.
  stirling=[[1]]
  for r in range(1,R+1):
   prev=stirling[-1];stirling.append([(prev[j-1] if j else 0)+(j*prev[j] if j<len(prev) else 0) for j in range(r+1)])
  for r in range(R+1):
   need(counts[r]==sum(stirling[r][j]*distinct[j] for j in range(min(r,n)+1)),'repeated ports are partitions of distinct roles')
  rows.append({'order':n,'ordered_port_counts':list(map(int,counts)),
               'distinct_port_counts':list(map(int,distinct)),
               'fixed_point_weight_coefficients':{str(j):str(w) for j,w in sorted(coeff.items())} if n in (11,12) else None})
 result={'status':'classical Burnside and inherited THM-514 extension, not a record or a class census',
         'maximum_order':N,'maximum_ordered_ports':R,'core_partition_terms':core_terms,
         'core_kernel_cells':len(core),'rows':rows,'checks':CHECKS}
 if write:(HERE/'counting.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'order':N,'core_terms':core_terms,'kernel_cells':len(core),
                   'n11_12':[(r['order'],r['ordered_port_counts'][:3]) for r in rows if r['order'] in (11,12)]}))
 return result
if __name__=='__main__':main()
