"""Independent path-permutation bank and single-pass return census."""
from itertools import permutations
from pathlib import Path
import json
BASE=Path('/workspace/shared/session184/spectrum_four')
def require(x,m):
 if not x:raise RuntimeError(m)
def paths(g,V):
 return [p for p in permutations(V) if all(g[a]>>b&1 for a,b in zip(p,p[1:]))]
def signature(g,W):
 n=len(g);full=list(range(n));hp=paths(g,full);P=F=G=J=0
 for mask in range(1<<n):
  a=paths(g,[i for i in full if mask>>i&1]);b=paths(g,[i for i in full if not mask>>i&1])
  for u in a:
   for v in b:
    P+=1;ef=bool(u) and u[-1] in W;sg=bool(v) and v[0] in W
    F+=ef;G+=sg;J+=ef and sg
 return (len(hp),P,F,G,J,sum(p[-1] in W for p in hp),sum(p[0] in W for p in hp))
def combine(a,b):
 if a is None:return b
 h,p,f,g,j,e,s=a;H,P,F,G,J,E,S=b
 return h*H,p*P,H*f+p*F,P*g+h*G,H*j+h*J+e*S+g*F,h*E,H*s
bank=[]
for typ,g,marks in [('S',(0,),[set(),{0}]),('C',(2,4,1),[set(range(a)) for a in range(4)]),('Q',(6,12,8,1),[{i for i in range(4) if mask>>i&1} for mask in range(16)])]:
 for a,W in enumerate(marks):bank.append(((typ,a),signature(g,W)))
r=json.loads((BASE/'receipt.json').read_text());targets={}
for tr in r['transitions']:
 for scan in tr['complete_scans']:targets.setdefault((scan['F'],scan['c']),[])
maxF=max(f for f,c in targets);wantedF={f for f,c in targets};stack=[((),None)];visited=0;terminal=0
while stack:
 w,z=stack.pop();visited+=1
 for letter,b in bank:
  if not w and not letter[1]:continue
  zz=combine(z,b)
  if zz[2]>maxF:continue
  ww=w+(letter,)
  if letter[1] and zz[2] in wantedF:
   terminal+=1;key=zz[2],zz[0]-zz[4]
   if key in targets:targets[key].append((ww,zz))
  if zz[1]<=maxF:stack.append((ww,zz))
for tr in r['transitions']:
 actual=[]
 for scan in tr['complete_scans']:
  for word,z in targets[(scan['F'],scan['c'])]:actual.append((scan['unmarked_Q_prefixes'],word,z))
 expected=[(h['unmarked_Q_prefixes'],tuple(tuple(b) for b in h['word']),tuple(h['signature'])) for h in tr['hits']]
 require(sorted(actual)==sorted(expected),('return mismatch',tr['source_address'],tr['digit']))
exception=[tr for tr in r['transitions'] if tr['target_constant']==-35]
require(len(exception)==1 and {x['unmarked_Q_prefixes'] for x in exception[0]['complete_scans']}=={0,1},'C=-35 coverage')
print(json.dumps({'status':'PASS','independent_component_signatures':len(bank),'distinct_target_pairs':len(targets),'maximum_F':maxF,'one_pass_prefixes':visited,'tested_terminal_candidates':terminal,'transitions_rechecked':len(r['transitions']),'exception_unmarked_Q_prefix_counts':[0,1]},sort_keys=True,indent=2))
