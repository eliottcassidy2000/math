"""Freeze checks and complete named-bank inversion, plus new literal query cross-checks."""
from pathlib import Path
import json,hashlib,collections
D=Path(__file__).resolve().parent

def need(c,s):
 if not c:raise RuntimeError(s)
def locate(a,b):
 p=D.parent/a
 return p if p.exists() else D.parent/b

def main():
 O=locate('odd_transitive','session189-odd-transitive');A=locate('collision_next','session189-collision-next');sources={};allbanks=[]
 for label,P in [('odd',O),('anchor',A)]:
  release=json.loads((P/'RELEASE_FILES.json').read_text());pins={}
  for name,r in release.items():
   data=(P/name).read_bytes();pin=hashlib.sha256(data).hexdigest();need(len(data)==r['bytes'] and pin==r['sha256'],'each frozen allowlisted source byte');pins[name]=pin
  sources[label]=pins
  parts=sorted(P.glob('bank.part*.json')) if label=='odd' else []
  bank=sum([json.loads(p.read_text()) for p in parts],[]) if parts else json.loads((P/'bank.json').read_text());allbanks.extend((label,b) for b in bank)
 degree=collections.Counter();byrows={};coords=0
 for label,b in allbanks:
  R=b['rows'];n=len(R);m=[0]+[1<<i for i in range(n)]+[x for x in range(1<<n) if x.bit_count()==2];need(m==b['masks'],'complete canonical physical named basis');L=len(m);ix={s:i for i,s in enumerate(m)};need(len(b['coefficients'])==len(b['values'])==L*L,'every named coefficient and basis query retained')
  for ai,a in enumerate(m):
   for ci,c in enumerate(m):
    coeff=0;aa=a
    while True:
     cc=c
     while True:
      coeff+=(-1)**(a.bit_count()+c.bit_count()-aa.bit_count()-cc.bit_count())*b['values'][ix[aa]*L+ix[cc]]
      if cc==0:break
      cc=(cc-1)&c
     if aa==0:break
     aa=(aa-1)&a
    need(coeff==b['coefficients'][ai*L+ci],'independent full Mobius inversion including zero cubic and quartic slots');degree[a.bit_count()+c.bit_count()]+=1;coords+=1
  byrows.setdefault(tuple(R),[]).append((label,b))
 querydata=json.loads((D/'literal_queries.json').read_text());matches=[]
 for j,q in enumerate(querydata):
  for label,b in byrows.get(tuple(q['old_rows']),[]):
   m=b['masks'];L=len(m);a,c=q['physical_query'];v=sum(b['coefficients'][ai*L+ci] for ai,s in enumerate(m) if s&a==s for ci,t in enumerate(m) if t&c==t);need(v==q['complete_q'][0],'new literal scalar query agrees with complete existing named bank');matches.append(dict(query_id=j,packet=label,graph_id=b['id'],value=v))
 need(matches,'independent actual-query comparisons are nonempty')
 outputs={'source_manifest_recount.json':sources,'named_bank_query_matches.json':matches}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'frozen source and query match data stable')
  else:p.write_text(txt)
 return dict(status='PASS',full_named_graphs=len(allbanks),complete_named_coordinates=coords,degree_counts=dict(sorted(degree.items())),new_literal_query_bank_matches=len(matches),scope='Every existing named coefficient is independently inverted and source-pinned; full independent recount of all existing basis values is supplied by the separate primary audit. New literal queries here are recounted from scratch.',source_sha256=hashlib.sha256((D/'bank_integrity.py').read_bytes()).hexdigest(),output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
