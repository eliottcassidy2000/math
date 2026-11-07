from itertools import permutations, product
from fractions import Fraction
from collections import Counter
from pathlib import Path
import importlib.util,json
spec=importlib.util.spec_from_file_location('carrier',str(Path(__file__).with_name('carrier_support.py')))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def need(x,msg):
 if not x: raise RuntimeError(msg)

def gate(a):
 n=len(a)
 return [x|(1<<n) for x in a]+[1<<(n+1),(1<<n)-1]

def poly_from_hist(hist,n):
 ans=[Fraction(0) for _ in range(n+1)]
 for f,num in hist.items():
  p=[Fraction(num)]
  for j in range(n):
   np=[Fraction(0)]*(len(p)+1)
   for i,v in enumerate(p):np[i]+=v*(f-j);np[i+1]+=v
   p=[v/Fraction(j+1) for v in np]
  ans=[x+y for x,y in zip(ans,p)]
 need(all(x.denominator==1 for x in ans),'integral coefficients')
 return list(map(int,ans))

def sub_counts(a):
 n=len(a);out=[]
 for s in range(1<<n):
  vs=[v for v in range(n) if s>>v&1]
  aa=[sum((a[u]>>v&1)<<j for j,v in enumerate(vs)) for u in vs]
  out.append(c.hamilton(aa))
 return out

rows=[]
for name,a in [('A',[58,4,1,6,14,30]),('B',[38,8,2,5,15,30])]:
 n=len(a)
 need(all(not (a[u]>>u&1) for u in range(n)),'loops')
 need(all((a[u]>>v&1)+(a[v]>>u&1)==1 for u in range(n) for v in range(u+1,n)),'tournament')
 hist=Counter(sum(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(n)))
 coefficients=poly_from_hist(hist,n)
 cyc=c.cycles_from_mask(n,c.encode(a),c.potential_cycles(n))
 packing,edges=c.small_packing(cyc)
 betti=c.graph_betti(len(cyc),edges)
 sub=sub_counts(a);evaluations=[]
 for q in range(5):
  total=0
  for alloc in product(range(q),repeat=n):
   bins=[0]*q
   for i,b in enumerate(alloc):bins[b]|=1<<i
   value=1
   for s in bins:value*=sub[s]
   total+=value
  need(total==sum(v*q**j for j,v in enumerate(coefficients)),f'bin sum {name} q{q}')
  evaluations.append(total)
 b=a[:];nested=[]
 for r in range(3):nested.append(c.hamilton(b));b=gate(b)
 row={'name':name,'outgoing_row_bitmasks':a,'mask':c.encode(a),'forward_histogram':[hist[f] for f in range(n)],'P_coefficients_low_to_high':coefficients,'direct_bin_values_q0_to4':evaluations,'packing_coefficients':packing,'quotient_betti_Q_and_F2':betti,'nu':len(packing)-1,'nested_H_r0_to2':nested,'odd_cycles':[list(w) for _,w in cyc]}
 rows.append(row)
need(rows[0]['forward_histogram']==rows[1]['forward_histogram']==[17,105,238,238,105,17],'hist claim')
need(rows[0]['P_coefficients_low_to_high']==rows[1]['P_coefficients_low_to_high']==[0,0,8,0,8,0,1],'polynomial claim')
need(rows[0]['packing_coefficients']==[1,8],'A cycles')
need(rows[1]['packing_coefficients']==[1,6,1],'B cycles')
need(rows[0]['quotient_betti_Q_and_F2']==[8,0],'A topology')
need(rows[1]['quotient_betti_Q_and_F2']==[5,1],'B topology')
need(rows[0]['nested_H_r0_to2']==rows[1]['nested_H_r0_to2']==[17,241,4001],'nested counts')
result={'status':'PASS','scope':'Two specified order-six tournaments only; original frozen projective packet untouched','rows':rows,'all_nested_equality':'Follows from exact equal P and inherited gate operator; only r=0,1,2 were also counted directly.','minimality':'Below order6 no two nontrivial odd cycles can be disjoint, so I(z)=1+((P(1)-1)/2)z and the projective carrier is that many isolated points.'}
text=json.dumps(result,indent=2)+'\n';Path(__file__).with_name('RESULT.json').write_text(text);print(text)
