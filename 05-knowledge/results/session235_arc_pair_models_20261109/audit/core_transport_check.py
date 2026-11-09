#!/usr/bin/env python3
"""Independent literal checks of the frozen complete source-core transport."""
import json,itertools,collections,math
from fractions import Fraction
from pathlib import Path
P=Path(__file__).resolve().parent
j=json.loads((P.parent/'census/core_transport.json').read_text())
def h(a):
 return sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in itertools.permutations(range(len(a))))
def automorphisms(a):
 return {p for p in itertools.permutations(range(len(a))) if all(a[i][k]==a[p[i]][p[k]] for i in range(len(a)) for k in range(len(a)))}
summary=collections.Counter();coretypes=collections.Counter()
for c in j['cores']:
 a=[list(map(int,r)) for r in c['source_matrix']];b=[list(map(int,r)) for r in c['decoded_matrix']];n=len(a)
 assert b==[[int(i!=k and not a[k][i]) for k in range(n)] for i in range(n)]
 source,target=h(a),h(b);aa,ab=automorphisms(a),automorphisms(b)
 assert source==c['source_H'] and target==c['decoded_H'] and aa==ab and len(aa)==c['automorphism_order']
 summary[source,n,target,len(aa)]+=1;coretypes[source]+=1
assert dict(coretypes)=={3:2,7:3,21:32}
assert {r for (s,n,r,a) in summary if s==7}=={1}
assert {r for (s,n,r,a) in summary if s==21}=={1,3,5,7,15}
# Independent explicit source-word enumeration for each n checks GF coefficients.
checks=0
for n in range(1,13):
 for source in (7,21):
  cls=collections.Counter();lab=collections.Counter()
  for c in j['cores']:
   if c['source_H']==source and n>=c['n']:
    ways=n-c['n']+1;r=c['decoded_H'];cls[r]+=ways;lab[r]+=Fraction(ways,c['automorphism_order'])*math.factorial(n)
  if source==21 and n>=7:
   for a in [c for c in j['cores'] if c['source_H']==7]:
    for b in [c for c in j['cores'] if c['source_H']==3]:
     ways=2*sum(1 for x in range(n-6) for y in range(n-6-x))
     r=a['decoded_H']*b['decoded_H'];cls[r]+=ways;lab[r]+=Fraction(ways,a['automorphism_order']*b['automorphism_order'])*math.factorial(n)
  pred=collections.Counter();predlab=collections.Counter()
  for t in j['generating_function_terms']:
   if t['source_H']==source and n>=t['numerator_order']:
    ways=math.comb(n-t['numerator_order']+t['denominator_power']-1,t['denominator_power']-1)
    pred[t['decoded_H']]+=ways*t['isomorphism_coefficient'];predlab[t['decoded_H']]+=ways*Fraction(t['labelled_EGF_coefficient'])*math.factorial(n)
  assert cls==pred and lab==predlab,(n,source)
  checks+=1
report={'status':'PASS','source_cores':len(j['cores']),'all_core_paths_and_automorphisms_recomputed':True,'source_word_GF_checks_through_order':12,'checks':checks,'grouped_cores':[{'source_H':s,'n':n,'decoded_H':r,'automorphisms':a,'classes':c} for (s,n,r,a),c in sorted(summary.items())]}
(P/'core_transport_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
