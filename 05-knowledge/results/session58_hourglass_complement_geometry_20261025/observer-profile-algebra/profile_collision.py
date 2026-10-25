import itertools,json
from pathlib import Path
from analyze import req

def fwht(a):
 h=1
 while h<len(a):
  for i in range(0,len(a),h*2):
   for j in range(i,i+h):
    x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
  h*=2
 return a

def main():
 n=6;pairs=list(itertools.combinations(range(n),2));ix={e:k for k,e in enumerate(pairs)};N=1<<len(pairs)
 data=json.loads(Path(__file__).with_name('walsh_certificate.json').read_text())[-1];basis=[];ret=[]
 for row in data['rows']:
  x=row['mask'];bits=[k for k in range(len(pairs)) if x>>k&1];coef={}
  for p in itertools.permutations(range(n)):
   y=sum(1<<ix[tuple(sorted((p[pairs[k][0]],p[pairs[k][1]])))] for k in bits)
   sign=(-1)**sum(p[pairs[k][0]]>p[pairs[k][1]] for k in bits)
   req(y not in coef or coef[y]==sign,'even graph sign consistency');coef[y]=sign
  a=[0]*N
  for y,s in coef.items():a[y]=s
  basis.append(fwht(a));ret.append(row['signed_order_sum']!=0)
 classes={};profiles={}
 for x in range(N):
  full=tuple(a[x] for a in basis);visible=tuple(a[x] for a,keep in zip(basis,ret) if keep)
  classes.setdefault(full,x);profiles.setdefault(visible,{})[full]=classes[full]
 collisions=[list(d.values()) for d in profiles.values() if len(d)>1]
 out={'n':n,'full_classes':len(classes),'distinct_profiles':len(profiles),'collisions':collisions,'omitted_graph_masks':[r['mask'] for r in data['rows'] if not r['signed_order_sum']]}
 req(len(classes)==34,'complete signed basis class count')
 Path(__file__).with_name('collision_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
