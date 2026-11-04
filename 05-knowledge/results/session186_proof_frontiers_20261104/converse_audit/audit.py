"""Independent literal permutation audit, no imports from the author's checker."""
import itertools,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent/'converse'
def need(x,m):
 if not x:raise RuntimeError(m)
def poly(a):
 n=len(a);c={}
 for word in itertools.permutations(range(n+2)):
  pos=neg=0
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not a[u]>>v&1:break
   elif u>=n and v>=n:
    if u>v:break
   else:
    old=u if u<n else v;port=v if v>=n else u;b=1<<(old+n*(port-n))
    if u<n:pos|=b
    else:neg|=b
  else:
   sub=neg
   while True:
    m=pos|sub;c[m]=c.get(m,0)+(-1 if sub.bit_count()%2 else 1)
    if not sub:break
    sub=(sub-1)&neg
 return c
def covers(a):
 n=len(a);C=[[0]*n for _ in a]
 for p in itertools.permutations(range(n)):
  bad=[i for i in range(n-1) if not a[p[i]]>>p[i+1]&1]
  if len(bad)>1:continue
  splits=bad if bad else range(n-1)
  for i in splits:C[p[i]][p[i+1]]+=1
 return C
def run():
 data=json.loads((ROOT/'receipt.json').read_text());result=[]
 need(data==json.loads((ROOT/'receipt.optimized.json').read_text()),'normal optimized equality')
 for row in data['cases']:
  a=row['rows'];n=len(a);pa=poly(a)
  need(covers(a)==row['full_named_C'],'literal source C')
  need(covers(row['module_reversed_rows'])==row['reversed_full_named_C'],'literal reversed C')
  counts=[]
  for key,graphkey in [('full_named_module_reversal_polynomial','module_reversed_rows'),('full_named_global_converse_cubic','global_converse_rows')]:
   pb=poly(row[graphkey]);diff={k:pa.get(k,0)-pb.get(k,0) for k in set(pa)|set(pb)};diff={k:v for k,v in diff.items() if v}
   reported={x|(y<<n):v for x,y,v in row[key]['nonzero']}
   need(diff==reported,'ALL named coefficients '+key)
   need(max(map(int.bit_count,diff))==row[key]['degree'],'degree')
   counts.append(len(diff))
  result.append({'state':row['state'],'literal_C_match':True,'all_named_ambient_coefficients':counts[0],'all_named_global_coefficients':counts[1]})
 return {'status':'PASS: independent literal permutation count agrees with every supplied coefficient and C coordinate','cases':result}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
