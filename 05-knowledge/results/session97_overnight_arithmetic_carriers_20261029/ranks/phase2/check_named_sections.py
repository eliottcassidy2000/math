"""Three explicitly named generic sections; specialization proves independence."""
import sys,json
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from check_characters import check,require
from named_family import family
P=Path(__file__).resolve().parents[1]
trans=json.loads((P/'specialization_certificates.json').read_text());out=[]
for z in trans:
 i=z['curve'];t=Q(z['T']);a,pts=family(i,t);u,r,s,w=[Q(z[k])for k in ['u','r','s','w']]
 c=json.loads((P/f'curve{i}.json').read_text());minimal=[]
 for x,y in pts:
  require(y*y+a[0]*x*y+a[2]*y==x**3+a[1]*x*x+a[3]*x+a[4])
  X=(x-r)/u**2;Y=(y-s*u*u*X-w)/u**3;minimal.append([str(X),str(Y)])
 c['points']=minimal;c['rank_lower_bound']=3;cert=check(c);require(cert['character_rank']==3)
 out.append({'curve':i,'minimal_points':minimal,'certificate':cert})
print(json.dumps(out,indent=2))
