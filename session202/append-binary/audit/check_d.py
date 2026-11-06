"""Independent D-interface controls. Imports only this audit's check.py, no author code."""
import contextlib,io,json,runpy
from pathlib import Path
ROOT=Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):x=runpy.run_path(str(ROOT/'check.py'))
need=x['need'];F=x['F'];apply=x['apply'];graph=x['graph'];count=x['count'];seen=x['seen']
limit=3333
plus={D for D,E,U,d in seen if 3<=D<=limit and U-D-E==1}
minus={D for D,E,U,d in seen if 3<=D<=limit and U-D-E==-1}
need(plus=={n for n in x['plus'] if 3<=n<=limit},'D plus spectrum')
need(minus=={n for n in x['minus'] if 3<=n<=limit},'D minus spectrum')
need(plus&minus=={5,9,17},'D intersection')
for s in seen:
 D,E,U,d=s
 need(U>=D and 2*E<=D-1,'physical D inequalities')
 need(2*F(s)>=3*D+1,'D growth')
literal=[]
for n,wm,wp in [(5,'Z','O'),(9,'ZZ','ZO'),(17,'ZZZ','ZJO')]:
 for w,sign in [(wm,-1),(wp,1)]:
  parent=w+'Z';child=parent+'J'
  need(apply(parent)[0]==n and apply(child)[0]==2*n+sign,'D algebra')
  pp,pc=count(graph(parent));cp,cc=count(graph(child))
  need(pp==n and cp==2*n+sign and cc==1,'D literal graphs')
  literal.append({'parent':n,'sign':sign,'parent_word':parent,'parent_paths':pp,'parent_cycles':pc,'child_word':child,'child_paths':cp,'child_cycles':cc})
need(apply('Z')[0]==3 and F(apply('Z'))==5,'D3 boundary')
print(json.dumps({'status':'PASS','physical_profiles':len(seen),'cover_cap':5000,'complete_D_parent_cap':limit,'plus_count':len(plus),'minus_count':len(minus),'both_child_parents':sorted(plus&minus),'literal_controls':literal,'D_parent3':'minus available; plus unavailable','scope':'Finite controls support separate all-order proof; no eventual-coverage claim.'},indent=2,sort_keys=True))
