"""Independent C1 sublevel inventory using only three one-vertex source moves."""
from pathlib import Path
import json
from monotone import inventory
P=Path(__file__).resolve().parent
ONE=(1,1,1,1)
TWO=(1,0,2,1)
def readout(v):D,E,U,d=v;return D-E+U

def three_letter(cap,retain=False):
 seen={TWO};queue=[TWO];values=set()
 for D,E,U,d in queue:
  h=D-E+U;values.add(h)
  for w in [(D,0,2*U,d),(D-E+U,d,2*U-d,d),(D-E+U,D-E,D+U,D)]:
   if readout(w)<=cap and w not in seen:seen.add(w);queue.append(w)
 if retain:return values,seen|{ONE},len(seen)+1
 return values,None,len(seen)+1

v,s,n=three_letter(5001,True)
ref,rs,rn=inventory(5001,True)
if s!=rs or v!=set(ref):raise RuntimeError('independent capped source-set mismatch')
del s,rs,ref
values,_,count=three_letter(50001)
expected=json.loads((P/'VALIDATION.json').read_text())
missing=[h for h in range(3,50002,2)if h not in values]
if count!=expected['source_states'] or missing!=expected['C1_missing'] or len(values)!=expected['C1_values']:
 raise RuntimeError((count,missing,len(values)))
r=dict(status='PASS',method='three one-vertex O/Z/J moves, no indexed gadget constants',comparison_bound=5001,exact_compared_states=n,bound=50001,source_states=count,C1_values=len(values),C1_missing=missing)
(P/'INDEPENDENT_THREE_LETTER_AUDIT.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,sort_keys=True))
