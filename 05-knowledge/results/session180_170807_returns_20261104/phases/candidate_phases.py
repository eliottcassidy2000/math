"""Recompute inherited field codecs, never infer them from source phase.
Codec copied mathematically from pinned return6325_completion/transport.py.
Input counts independently supplied by ports/joint.cpp.
"""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).parent

def need(x):
 if not x:raise RuntimeError('codec check failed')
def mul(a,b):
 c=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for i in (4,3):c[i-1]-=c[i];c[i-3]-=c[i]
 return tuple(x%7 for x in c[:3])
def power(a,k):
 out=(1,0,0)
 while k:
  if k&1:out=mul(out,a)
  a=mul(a,a);k//=2
 return out
def chart(g):
 h,x,y=g;z=tuple(v%7 for v in g);a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
 need(power(z,343)==z);need(((b+c)%7,(a-4*c)%7,b)==z)
 return {'observer_mod7':z,'F343_vector_codec':[z,power(z,7),power(z,49)],'F343_codec_degree':next(i for i in (1,3) if power(z,7**i)==z),'F49_charge':[a,b,c],'F49_frobenius':[(a+4*b)%7,(-b)%7,c]}
def main():
 # Validate local implementation against every historical saved state.
 old=json.loads((HERE.parent/'ports/input_audit.json').read_text())
 for s in old['states']:
  need(json.loads(json.dumps(chart(s['golden'])))==s['charts'])
 candidate=json.loads((HERE.parent/'returns/candidate.json').read_text())
 names=candidate.get('names',candidate.get('target_names'))
 lines=(HERE.parent/'ports/candidate.out').read_text().splitlines();h=int(lines[0].split()[1]);records=[]
 for line in lines[1:]:
  z=line.split()
  if z[0]!='PORT':continue
  idx,e1,s1,delh,e2,s2=map(int,z[1:]);need((s2-e1)%2==0);g=[h,(s2-e1)//2,(s2+e1)//2]
  records.append({'name':names[idx],'index':idx,'H':h,'E1':e1,'S2':s2,'golden':g,'charts':chart(g)})
 result={'status':'EXACT_CODECS_FROM_INDEPENDENT_ENDPOINT_COUNTS','modulus':[1,0,1,1],'source_ref':'4ad839b20dd280cdcfbadf86cb62b2be3c9615ad','historical_codec_controls':len(old['states']),'states':records,'scope':'Current observer charts; not inherited spectral-clock transport, automorphisms, or permission for phase resets.'}
 result['input_sha256']={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE.parent/'ports/input_audit.json',HERE.parent/'ports/candidate.out',HERE.parent/'returns/candidate.json')}
 (HERE/'candidate_phases.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
