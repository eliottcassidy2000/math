"""Exact finite-kernel leading contraction; no tournament samples used."""
import json
from fractions import Fraction
from collections import Counter, defaultdict
from pathlib import Path
src=Path(__file__).resolve().parents[1]/'session30-decorated-engine'/'symbolic-kernels.json'
cores=json.loads(src.read_text());total=defaultdict(Fraction);rows=[]
for ci,c in enumerate(cores):
 if c['weight']!=2:continue
 edges=[tuple(e) for e in c['edges']];deg=Counter()
 for a,b in edges:deg[a]+=1;deg[b]+=1
 leaves={i for i,v in deg.items() if v==1}
 leafat=Counter();mult=Counter()
 for a,b in edges:
  if a in leaves:leafat[b]+=1
  elif b in leaves:leafat[a]+=1
  elif a!=b:mult[(a,b)]+=1
 # Any odd simple S-edge leaves a vanishing colored signed density.
 # Otherwise odd remaining d weights average to zero.
 vanish=any(k%2 for k in mult.values()) or any(k%2 for k in leafat.values())
 poly=defaultdict(Fraction)
 if not vanish:
  for t in c['terms']:
   coeff=Fraction(t['numerator'],c['denominator']);power=0
   for (a,b),j in zip(edges,t['indices']):
    if a in leaves or b in leaves:
     # Row sum R+ = +zu d, R- = -zu d; column changes sign.
     # Since matrices' orientations follow edges, a leaf at a needs column sum.
     if j>=2:coeff=0;break
     coeff*= (1 if j==0 else -1)*( -1 if a in leaves else 1)
     power+=1
    elif a==b:
     coeff*=(-1,-1,1,-1)[j];power+=2
    else:
     if j>=2:coeff=0;break
     coeff*=(1,-1)[j];power+=1
   poly[power]+=coeff
 poly={k:v for k,v in poly.items() if v}
 for k,v in poly.items():total[k]+=v
 row={'core_index':ci,'V':c['V'],'edges':edges,'excess':c['excess'],'leaves':c['leaves'],'leading':{str(k):str(v) for k,v in sorted(poly.items())}}
 rows.append(row)
 print(row)
print('TOTAL',dict(total))
Path(__file__).with_name('a2-leading-certificate.json').write_text(json.dumps({'source':str(src),'method':'exact integer four-basis kernel contraction after leaf elimination and flat-spectrum leading reduction','cores':rows,'total':{str(k):str(v) for k,v in sorted(total.items())}},indent=2)+'\n')
assert dict(total)=={6:Fraction(26),4:Fraction(7)},dict(total)
