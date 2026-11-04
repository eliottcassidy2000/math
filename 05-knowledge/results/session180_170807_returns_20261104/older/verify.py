import itertools,json,pathlib,hashlib
P=pathlib.Path(__file__).parent
atlas=json.loads((P/'prefix_atlas.json').read_text()); pairs={(x['h'],x['p2']) for x in atlas['pairs']}
def literal(n,code):
 rows=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1: rows[i]|=1<<j
  else: rows[j]|=1<<i
 h=p2=0
 for pi in itertools.permutations(range(n)):
  bad=[k for k in range(n-1) if not rows[pi[k]]>>pi[k+1]&1]
  if not bad: h+=1;p2+=n+1
  elif len(bad)==1:p2+=1
 return h,p2
for x in atlas['pairs']:
 got=literal(x['n'],x['code'])
 if got!=(x['h'],x['p2']):raise RuntimeError((x,got))
# Independent permutation enumeration of the complete labelled universe n<=5.
small=set()
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
  h,p2=literal(n,code)
  if p2<=198:small.add((h,p2))
if small!={(x['h'],x['p2']) for x in atlas['pairs'] if x['n']<=5}:raise RuntimeError('small atlas mismatch')
rows=[json.loads(x) for x in (P.parent/'returns/hits.jsonl').read_text().splitlines()]
candidates=[x for x in rows if x.get('lattice')]
survivors=[x for x in candidates if (x['hA'],x['P2A']) in pairs]
out={'status':'PASS','atlas_rows':atlas['labelled_rows'],'witnesses_permutation_checked':len(atlas['pairs']),'all_labelled_orders_permutation_crosschecked':[1,2,3,4,5],'candidate_count':len(candidates),'candidate_max_p2':max(x['P2A'] for x in candidates),'realizable_candidates':survivors,'cxx_sha256':hashlib.sha256((P/'prefix_atlas.cpp').read_bytes()).hexdigest()}
(P/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
