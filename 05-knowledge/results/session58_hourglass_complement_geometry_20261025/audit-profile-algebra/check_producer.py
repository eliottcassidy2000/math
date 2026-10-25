from pathlib import Path
import json,hashlib
from fractions import Fraction as F
D=Path(__file__).resolve().parent.parent/'observer-profile-algebra';A=Path(__file__).parent
p=json.loads((D/'direct_certificate.json').read_text());a=json.loads((A/'independent_results.json').read_text())['census'][-1]
def req(v,m):
 if not v:raise RuntimeError(m)
req(p['representatives']==a['reps'],'representatives');req(p['columns']==a['cols'],'columns');req(p['matrix']==a['matrix'],'literal full matrix')
pos={int(k):F(v) for k,v in p['mixture_positive'].items()};neg={int(k):F(v) for k,v in p['mixture_negative'].items()}
req(sum(pos.values())==sum(neg.values())==1,'probability normalization');req(not set(pos)&set(neg),'disjoint supports')
for j in range(len(a['cols'])):
 req(sum((pos.get(r,0)-neg.get(r,0))*a['matrix'][i][j] for i,r in enumerate(a['reps']))==0,'independent literal mixture equality')
print(json.dumps({'entire_literal_matrix_matches':True,'mixtures_verified_against_independent_matrix':True,'columns':len(a['cols']),'rows':len(a['reps'])},indent=2))
