"""Independent old-coordinate evaluation of the complete short menu and graphs."""
from pathlib import Path
from itertools import product,permutations
import sys,json
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'dependencies'/'dependencies'))
import bilinear,c3

def move(v,k):
 D,E,U,d=v
 if k=='O':return D,0,2*U,d
 if k=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

def word(v,w):
 for letter in w:v=move(v,letter)
 return v

v,t=(5,1,5,1),(1,0,2,1)
W=[w for n in range(3)for w in product('OZJ',repeat=n)]
values=[]
for a in W:
 for b in W:
  if len(a)+len(b)<=2:values.append(bilinear.value('N',word(v,a),word(t,b)))
if len(values)!=34:raise RuntimeError('N count')
values.append(bilinear.value('M',v,t))
for k in ['J0','J1','J2','J3']:
 for a,b in product(range(3),repeat=2):values.append(c3.value(k,word(v,'J'*a),word(t,'J'*b)))
r=json.loads((P/'VALIDATION.json').read_text())
if len(values)!=71 or sorted(set(values))!=r['bad_pair']['complete_candidate_values'] or 81 in values or 83 in values:
 raise RuntimeError('branch mismatch')
graphs=[]
for q in r['literal_rebase_graphs']:
 T=q['rows'];n=len(T)
 H=sum(all(T[a]>>b&1 for a,b in zip(p,p[1:]))for p in permutations(range(n)))
 C=sum(all(T[a]>>b&1 for a,b in zip((0,)+p,p+(0,)))for p in permutations(range(1,n)))
 if(H,C)!=(q['H'],q['C']):raise RuntimeError(('graph',q,H,C))
 graphs.append(dict(label=q['label'],H=H,C=C,order=n))
out=dict(status='PASS',candidate_branches=len(values),distinct_values=len(set(values)),values=sorted(set(values)),literal_graphs=graphs)
(P/'INDEPENDENT_BRANCH_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
