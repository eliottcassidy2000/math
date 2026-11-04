"""Lossless complete Boolean coefficient banks for all complete physical query cubes."""
from pathlib import Path
import json,hashlib
D=Path(__file__).resolve().parent;P=D.parent/'even_transfer'
if not P.exists():P=D.parent/'session189-even-transfer'
def need(c,s):
 if not c:raise RuntimeError(s)
def mobius(v,bits):
 out=v.copy()
 for bit in range(bits):
  for mask in range(1<<bits):
   if mask>>bit&1:out[mask]-=out[mask^(1<<bit)]
 return out

def main():
 data=json.loads((P/'physical_records.json').read_text());groups={}
 for r in data['records']:
  if r['ambient'] or len(r['Q'])!=2:continue
  key=('even',tuple(r['Q']),tuple(r['child_type_indices']));groups.setdefault(key,[]).append(r)
 groups[('outside-domain',tuple([2,4,1]),())]=data['complete_outside_domain_bank'];banks=[];coordinates=0;quartics=cubics=0
 for key,rs in groups.items():
  m=len(rs[0]['Q']);need(len(rs)==1<<(2*m),'complete physical block query cube');bycode={r['block_query'][0]|(r['block_query'][1]<<m):r for r in rs};need(set(bycode)==set(range(1<<(2*m))),'every physical query unique')
  sides=[]
  for side in range(2):
   channels=[]
   for h in range(5):
    vals=[bycode[x]['sides'][side]['ordered_covers'][h] for x in range(1<<(2*m))];coef=mobius(vals,2*m)
    need(all(c==0 for mask,c in enumerate(coef) if (mask&((1<<m)-1)).bit_count()>2 or (mask>>m).bit_count()>2),'full bidegree with retained zero coordinates')
    for x,value in enumerate(vals):
     S=x;back=0
     while True:
      back+=coef[S]
      if S==0:break
      S=(S-1)&x
     need(back==value,'every full Boolean query recovered from complete named bank')
    channels.append(dict(h=h+1,values=vals,coefficients=coef));coordinates+=len(coef);cubics+=sum(mask.bit_count()==3 for mask in range(len(coef)));quartics+=sum(mask.bit_count()==4 for mask in range(len(coef)))
   sides.append(channels)
  banks.append(dict(scope=key[0],quotient_rows=list(key[1]),child_type_indices=list(key[2]),child_sets=rs[0]['sides'][0]['child_sets'],physical_bit_expansion=[sum(1<<v for v in S) for S in rs[0]['sides'][0]['child_sets']],variable_order=['a_'+str(i) for i in range(m)]+['c_'+str(i) for i in range(m)],all_cover_channels=sides,physical_queries=[bycode[x]['sides'][0]['physical_query'] for x in range(1<<(2*m))]))
 p=D/'complete_named_banks.json';txt=json.dumps(banks,separators=(',',':'))+'\n'
 if p.exists():need(p.read_text()==txt,'complete named banks stable')
 else:p.write_text(txt)
 return dict(status='PASS',complete_physical_query_families=len(banks),ordered_cover_channels=10*len(banks),all_named_coefficient_coordinates=coordinates,all_cubic_coordinates_including_zero=cubics,all_quartic_coordinates_including_zero=quartics,scope='Complete block-uniform polynomial banks, with fixed physical child sets and query masks; no within-child nonuniform observer is asserted.',source_sha256=hashlib.sha256((D/'named_audit.py').read_bytes()).hexdigest(),bank_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
