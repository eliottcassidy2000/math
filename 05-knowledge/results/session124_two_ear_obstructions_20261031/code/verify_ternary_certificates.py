"""Independent standard-library verifier: no producer import or symbolic algebra."""
from fractions import Fraction
import json,copy,itertools
from pathlib import Path
STATES=list(itertools.product(range(2),range(3),range(3),range(2)))
def advance(counts,c):
 result={}
 def add(s,n):result[s]=result.get(s,0)+n
 for (order,p,q,r),n in counts.items():
  first,second=(0,1) if order==0 else (1,0)
  fa=bool(c&(1<<first));fb=bool(c&(1<<second))
  add((order,2 if fa else 1,q,r),n)
  if q!=0 or fa:add((order,p,2 if fb else 1,r),n)
  if r!=0 or fb:add((order,p,q,1),n)
 return result
def apply(counts,word):
 for c in word:counts=advance(counts,c)
 return counts
def output(counts):
 return sum(n for (o,p,q,r),n in counts.items() if p!=2 and(q==1 or(q==0 and o==0)))
def check(cert):
 B=cert['radix'];bank=cert['bank'];checks=0
 for j,s in enumerate(STATES):
  total=Fraction(0)
  for term in cert['terms']:
   v=apply({s:1},term['word'])
   defect=sum(output(apply(v,w)) for w in bank)-B*B*output(v)
   total+=Fraction(term['coefficient'])*defect
  if total!=cert['target'][j]:return False
  expected=int(s[0]==0 and s[2]==0 and s[3]==0)
  if cert['target'][j]!=expected:return False
  checks+=1
 return checks
certs=json.loads(Path(__file__).with_name('ternary_certificates.json').read_text())
checks=0;mutations=0
for c in certs:
 if sum(Fraction(t['coefficient']) for t in c['terms'])!=0:raise RuntimeError('nonzero coefficient sum')
 result=check(c)
 if not result:raise RuntimeError('bad exact certificate')
 checks+=result
 bad=copy.deepcopy(c);bad['terms'][0]['coefficient']=str(Fraction(bad['terms'][0]['coefficient'])+1)
 if check(bad):raise RuntimeError('failed coefficient mutation')
 mutations+=1
 bad=copy.deepcopy(c);bad['target'][0]+=1
 if check(bad):raise RuntimeError('failed target mutation')
 mutations+=1
seed={(0,0,0,0):1,(1,0,0,0):1};seed=apply(seed,(1,))
children=[output(advance(seed,c)) for c in range(4)]
if children!=[3,1,5,1]:raise RuntimeError('convention hostile')
print(json.dumps(dict(exact_coordinate_checks=checks,rejected_mutations=mutations,seed1_children=children,independent_stdlib=True)))
