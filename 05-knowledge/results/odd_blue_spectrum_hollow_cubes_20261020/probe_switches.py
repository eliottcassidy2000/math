import sys,json,collections,argparse
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument('--census',type=Path,default=Path('session14-metagraph-census'))
args=ap.parse_args()
sys.path.insert(0,str(args.census))
from census import schema,tournament,reflection

def hp_one(a):
 n=len(a);full=(1<<n)-1
 def go(path,used):
  if used==full:return path
  opts=(a[path[-1]] if path else full)&(full^used)
  while opts:
   b=opts&-opts; opts-=b
   p=go(path+[b.bit_length()-1],used|b)
   if p:return p
 return go([],0)

def apply(a,kind):
 a=a[:];n=len(a)
 pairs=[]
 if kind&1:pairs += [(i,n//2) for i in range(n) if i!=n//2]
 if kind&2:pairs += [(i,n-1-i) for i in range(n//2)]
 for u,v in pairs:a[u]^=1<<v;a[v]^=1<<u
 return a

out=[]
for n in [3,5,7]:
 d=json.load(open(args.census/f'census_n{n}.json')); tt,ss=schema(n)
 br={}
 for t,u in enumerate(d['pi']):
  if reflection(t,ss)==t:br.setdefault(u,t)
 for kind in [1,2,3]:
  rows=[]
  for u,t in br.items():
   a=apply(tournament(t,n,tt),kind);p=hp_one(a)
   z=sum((1-((a[p[n-x]]>>p[n-y])&1))<<j for j,(x,y) in enumerate(tt)); v=d['pi'][z]
   r,s=d['nodes'][u],d['nodes'][v]
   rows.append(dict(u=u,v=v,t=t,source_blue=r['blue'],target_blue=s['blue'],source_c3=r['c3'],target_c3=s['c3']))
  report=dict(n=n,kind=kind,blue_preserved=all(r['source_blue']==r['target_blue'] for r in rows),parity_flipped=all((r['source_c3']-r['target_c3'])%2 for r in rows),is_involution=all(rows[[q['u'] for q in rows].index(r['v'])]['v']==r['u'] for r in rows),failed=[r for r in rows if r['source_blue']!=r['target_blue']][:5],rows=rows)
  out.append(report);print({k:v for k,v in report.items() if k!='rows'})
open(Path(__file__).resolve().parent/'switches.json','w').write(json.dumps(out,indent=2)+'\n')
