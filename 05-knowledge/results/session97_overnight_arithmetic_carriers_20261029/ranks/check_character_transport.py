"""Named local-character transport, with exact denominator/unit guards."""
import json
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent

def must(x):
 if not x:raise ValueError('transport check failed')
def mod(x,p):
 x=Q(x);must(x.denominator%p!=0);return x.numerator*pow(x.denominator,-1,p)%p
trans={d['curve']:d for d in json.loads((P/'specialization_certificates.json').read_text())}
out=[]
for cert in json.loads((P/'character_certificates.json').read_text()):
 i=cert['curve_id'];c=json.loads((P/f'curve{i}.json').read_text());t=trans[i];u,r=Q(t['u']),Q(t['r']);checked=0;excluded=[];missing_translation=0
 for row in cert['independent_rows']:
  p=row['p']
  if u.denominator%p==0 or u.numerator%p==0 or r.denominator%p==0:excluded.append(p);continue
  root=mod(u*u*row['root']+r,p);bits=0
  for k,(x,y) in enumerate(c['points']):
   x=Q(x);newx=u*u*x+r;d=(mod(newx,p)-root)%p;must(d!=0)
   bit=(pow(d,(p-1)//2,p)==p-1);bits|=bit<<k
   # Deliberately forget the root translation; count visibly wrong signs.
   wrong=(mod(newx,p)-mod(u*u*row['root'],p))%p
   if wrong and (pow(wrong,(p-1)//2,p)==p-1)!=bit:missing_translation+=1
  must(bits==row['bits']);checked+=1
 must(checked==len(c['points']))
 out.append({'curve':i,'rows_preserved':checked,'excluded_primes':excluded,'wrong_bits_if_root_translation_forgotten':missing_translation})
print(json.dumps(out,indent=2))
