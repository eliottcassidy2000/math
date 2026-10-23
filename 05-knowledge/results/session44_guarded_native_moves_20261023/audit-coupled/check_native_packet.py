"""Independent exact audit of a selected two-slot native packet.
No subject checker is imported. Run with Python or Python -O, optionally --repo PATH.
"""
from pathlib import Path
from itertools import permutations, combinations, product
from array import array
import argparse,csv,gzip,hashlib,json,sys
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path);args=p.parse_args()
here=Path(__file__).resolve()
candidates=([args.repo] if args.repo else [])+list(here.parents)+[here.parent.parent/'session35-publication/clean-repo']
repo=next((p for p in candidates if p and (p/'05-knowledge/results/metagraph_order8_20261023/owner_map.u16.gz').exists()),None)
if repo is None: raise RuntimeError('Supply --repo pointing to the pinned original atlas')
r=repo/'05-knowledge/results';f=r/'metagraph_order8_20261023'
checks=0

def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)

nodes=list(csv.DictReader((f/'merged_nodes.tsv').open(),delimiter='\t'))
owner=array('H');owner.frombytes(gzip.decompress((f/'owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(row[k]):int(row['merged_id']) for row in nodes for k in ('class_a','class_b')}
x=json.loads((r/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']['higher_layers'][-1]['lifts'][0]
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>i+1];full=(1<<21)-1

def matrix(c):
 A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(pairs):A[i][j]=(c>>(27-bit))&1;A[j][i]=1-A[i][j]
 return A

def kmatrix(A,P):
 keep=set(zip(P,P[1:]));return [[0 if u==v else A[u][v] if (u,v) in keep or (v,u) in keep else A[v][u] for v in range(8)] for u in range(8)]

def isomorphism(A,B):
 da=list(map(sum,A));db=list(map(sum,B));order=sorted(range(8),key=lambda u:db.count(da[u]));assign={}
 def visit(i):
  if i==8:return tuple(assign[u] for u in range(8))
  u=order[i]
  for v in range(8):
   if v in assign.values() or da[u]!=db[v]:continue
   if any(A[u][w]!=B[v][z] for w,z in assign.items()):continue
   assign[u]=v;z=visit(i+1)
   if z is not None:return z
   del assign[u]
  return None
 return visit(0)

source=1238;A=matrix(int(nodes[source]['canon_a']));records=[]
for P in permutations(range(8)):
 if not all(A[u][v] for u,v in zip(P,P[1:])):continue
 word=sum(A[P[i]][P[j]]<<b for b,(i,j) in enumerate(chords))
 need(merged[owner[word]]==source,'source owner')
 target=merged[owner[full^word]]
 records.append(dict(path=P,target=target,phase=(x[source]-x[target])%8,word=word))
need(len(records)==409,'complete source universe')
slots=[]
for start,end,target,phase in [(3,2,2397,5),(2,3,2790,4)]:
 rows=[a for a in records if a['path'][0]==start and a['path'][-1]==end and a['target']==target]
 need(len(rows)==2,'complete endpoint plus target fibre')
 need(all(a['phase']==phase for a in rows),'original marked phase')
 for a in rows:
  P=a['path'];K=kmatrix(A,P)
  need(kmatrix(K,P)==A,'literal complement involution')
  need(all(K[u][v] for u,v in zip(P,P[1:])),'same path in complement')
  row=nodes[target];witness=None
  for key in ('canon_a','canon_b'):
   if row.get(key) in (None,''):continue
   C=matrix(int(row[key]));witness=isomorphism(K,C)
   if witness is not None:break
  need(witness is not None,'independent target isomorphism, not only atlas lookup')
  a['target_isomorphism']=witness
  rev=P[::-1]
  reverse_word=sum(A[rev[j]][rev[i]]<<b for b,(i,j) in enumerate(chords))
  need(merged[owner[reverse_word]]==source and merged[owner[full^reverse_word]]==target,'literal converse with reversed clocks')
  need(P[2:4]==(6,7) and P[5:7]==(5,4),'two fixed indexed shared arcs')
 slots.append(rows)
states=[(a,b) for a,b in product(*slots) if a['path'][1:7]==b['path'][1:7]]
need(len(states)==2,'complete shared-frame state set')
need(len({a['path'] for slot in slots for a in slot})==4,'distinct actual path resources')
need(len({a['word'] for slot in slots for a in slot})==4,'distinct normalized presentations')
need(all(states[0][i]['path']!=states[1][i]['path'] for i in range(2)),'no one-slot edge')
changed=[[j for j in range(8) if states[0][i]['path'][j]!=states[1][i]['path'][j]] for i in range(2)]
need(changed==[[1,4],[1,4]],'four slot-clock edits, two shared variables')
need(all(states[0][i]['target']==states[1][i]['target'] for i in range(2)),'all per-slot harmonic characters preserved by literal equal targets')
inputs=[f/'merged_nodes.tsv',f/'owner_map.u16.gz',r/'metagraph_order8_electrical_20261023/receipt.json']
print(json.dumps({'status':'PASS','checks':checks,'source':source,'full_source_paths':len(records),'slots':slots,'compatible_states':[[a['path'] for a in st] for st in states],'changed_slot_clocks':changed,'minimum_slot_budget':2,'input_sha256':{str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'scope':'Selected observer packet, not a rewiring of the complete atlas; extension requires unused actual path resources or an explicit global permutation.'},indent=2))
