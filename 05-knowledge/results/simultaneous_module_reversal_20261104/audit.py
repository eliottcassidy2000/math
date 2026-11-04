"""Independent cover-form quartics, literal interactions, and exact cube audit."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session188_structural_bridges_20261104/response/verify.py'
spec=importlib.util.spec_from_file_location('incoming_quartics',SOURCE)
cover=importlib.util.module_from_spec(spec);spec.loader.exec_module(cover)
need=cover.need


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'audit.json');args=p.parse_args()
 data=json.loads((HERE/'full_response.json').read_text());cert=json.loads((HERE/'certificate.json').read_text())
 states=data['states'];quartics=0;converse_pairs=0
 for s in states:
  n=len(s['rows']);ms=s['masks'];form=cover.quartic_normal_form(s['rows'])
  for i,a in enumerate(ms):
   for j,c in enumerate(ms):
    if a.bit_count()==c.bit_count()==2:
     need(s['coefficients'][i*len(ms)+j]==form[a|(c<<n)],'independent named quartic');quartics+=1
  reverse=next(t for t in states if t['family']==s['family'] and t['rows']==cover.opposite(s['rows']))
  difference=[a-b for a,b in zip(s['coefficients'],reverse['coefficients'])]
  need(all(not z for k,z in enumerate(difference) if ms[k//len(ms)].bit_count()+ms[k%len(ms)].bit_count()==4),'whole-converse quartic cancellation')
  converse_pairs+=1
 family=next(f for f in cert['families'] if f['name']=='cyclic-three-pairs')
 cube={states[i]['epsilon']:states[i] for i in family['states'] if states[i]['top']==0}
 # The mixed flip derivative retains a genuine a0 a2 c1 c4 coefficient.
 A,B=5,18;literal=[];total=0
 a=A
 while True:
  c=B
  while True:
   vals=[cover.hp(cover.extend(cube[e]['rows'],a,c)) for e in (0,1,2,3)]
   sign=(-1)**(4-a.bit_count()-c.bit_count());mixed=vals[0]-vals[1]-vals[2]+vals[3]
   total+=sign*mixed;literal.append(dict(a=a,c=c,counts=vals,sign=sign,mixed=mixed))
   if not c:break
   c=(c-1)&B
  if not a:break
  a=(a-1)&A
 need(total==-2,'literal mixed quartic')
 ms=cube[0]['masks'];index=ms.index(A)*len(ms)+ms.index(B)
 need(family['two_flip_interaction_coefficients'][index]==total,'mixed polynomial matches literal')
 top=next(states[i] for i in family['states'] if states[i]['top']==1 and states[i]['epsilon']==0)
 q=(1,5);values=[cover.hp(cover.extend(s['rows'],*q)) for s in (cube[0],top)]
 need(values==[251,249],'nested full-response separator')
 need(family['quotient_only_uniform']['nonzero_mod4']==0,'uniform modular blindness')
 need(family['quotient_only_uniform']['nonzero_exact']==24,'uniform exact separator count')
 receipt=dict(status='PASS independent quartic/interaction audit; scope remains disjoint all-order, nested finite',
              source_sha256={str(SOURCE.relative_to(HERE.parent)):sha256(SOURCE.read_bytes()).hexdigest(),
                             'audit.py':sha256(Path(__file__).read_bytes()).hexdigest()},
              tables_sha256=sha256((HERE/'full_response.json').read_bytes()).hexdigest(),
              full_quartic_coordinates=quartics,whole_converse_pairs=converse_pairs//2,
              mixed_quartic=dict(a=A,c=B,coefficient=total,literal_counts=literal),
              nested_separator=dict(a=q[0],c=q[1],counts=values,rows=[cube[0]['rows'],top['rows']]))
 args.output.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in receipt.items() if k not in ('source_sha256','mixed_quartic','nested_separator')},indent=2))


if __name__=='__main__':main()
