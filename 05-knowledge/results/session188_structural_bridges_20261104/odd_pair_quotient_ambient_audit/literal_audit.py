"""Independent predecessor-DP replay of every retained ambient root/cube query."""
from pathlib import Path
import importlib.util,json,hashlib
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'odd_pair_quotient'
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
data=json.loads((SRC/'ambient.normal.json').read_text());v.need((SRC/'ambient.normal.json').read_bytes()==(SRC/'ambient.optimized.json').read_bytes(),'producer modes agree')
def root(Q,r,A):
 outside=[u for u in range(len(Q)) if u!=r];n=len(A);out=list(A)+[0]*len(outside)
 for i,u in enumerate(outside):
  out[n+i]=sum(1<<(n+j) for j,w in enumerate(outside) if Q[u]>>w&1)
  if Q[u]>>r&1:out[n+i]|=(1<<n)-1
  else:
   for w in range(n):out[w]|=1<<(n+i)
 return out
def coefficient(values,a,c):
 z=0;x=a
 while True:
  y=c
  while True:
   z+=(-1)**((a^x).bit_count()+(c^y).bit_count())*values[x,y]
   if not y:break
   y=(y-1)&c
  if not x:break
  x=(x-1)&a
 return z
counts=0;traces=[]
for r in data['contexts']:
 T,U=r['source'],r['target'];outside=sum(1<<u for u in range(6,len(T)) if T[u]&1);values={}
 for coarse,query,expected in zip(r['coarse_queries'],r['physical_queries'],r['source_target_values']):
  a,c=coarse;physical=[outside|sum(3<<(2*i) for i in range(3) if a>>i&1),outside|sum(3<<(2*i) for i in range(3) if c>>i&1)];v.need(physical==query,'physical root-incidence expansion');actual=[v.hp(v.extend(rows,*query)) for rows in (T,U)];v.need(actual==expected,'independent ambient root counts');values[a,c]=actual[0]-actual[1];counts+=2
 Q=[2,4,1];L=[sum(coefficient(values,1<<i,1<<j) for i in range(3) if i!=j)+sum(coefficient(values,(1<<i)|(1<<j),1<<j) for i in range(3) if Q[j]>>i&1) for j in range(3)]
 factor=v.hp(root(r['outer_quotient'],r['root'],[2,0]));v.need(factor==r['unit'] and factor%2 and L==r['L'] and sum(L)%8==4,'independent constant unit and cubic trace');traces.append(L)
cube=data['ambient_cube'];signatures=[]
for state in cube['states']:
 values=[v.hp(v.extend(state['physical_rows'],*query)) for query in cube['queries']];v.need(values==state['query_values'],'independent complete retained cube battery');counts+=len(values)
 for d in cube['decoders']:
  if d['theta']==state['theta']:v.need(values[d['query_index']]==d['values'][(state['bits']>>d['child'])&1],'pair decoder unaffected by other bits')
 signatures.append(tuple(values))
v.need(len(set(signatures))==16 and len(cube['queries'])==28 and cube['adaptive_bound']==25,'cube distinctness and query upper bounds')
print(json.dumps({'status':'PASS','literal_counts_replayed':counts,'contexts':len(data['contexts']),'root_traces':traces,'cube_states':len(signatures),'fixed_queries':28,'adaptive_bound':25,'method':'No source verifier import; independent predecessor DP and ordinary finite-difference extraction.','source_sha256':{p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ('AMBIENT.md','verify_ambient.py','ambient.normal.json','ambient.optimized.json')}},sort_keys=True,indent=2))
