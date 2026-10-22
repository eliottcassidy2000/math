import json,itertools,collections,pathlib
from experiment import frombits,bits,brute
ROOT=pathlib.Path(__file__).resolve().parent
rows=[json.loads(x) for x in open(ROOT/'data.jsonl')]
result={}
for n in [2,4,6]:
 pool={r['bits']:r for r in rows if r['name'].startswith(f'exhaustive-labelled-{n}-')};classes=[]
 while pool:
  b=min(pool);r=pool[b];S=frombits(n,b);orbit={bits(S[list(p)][:,list(p)]) for p in itertools.permutations(range(n))};assert orbit <= set(pool)
  assert brute(S)==int(r['H'])
  classes.append(dict(bits=b,H=r['H'],labelled_orbit=len(orbit),automorphisms=__import__('math').factorial(n)//len(orbit),y=r['y'],y_tree=r['y_tree']))
  for v in orbit:del pool[v]
 result[n]=classes
json.dump(result,open(ROOT/'small_classes.json','w'),indent=2)
print(json.dumps(result,indent=2))
