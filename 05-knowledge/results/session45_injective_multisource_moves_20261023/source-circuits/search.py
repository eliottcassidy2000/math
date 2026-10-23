import argparse,csv,gzip,json,sys,time
from pathlib import Path
from array import array
from itertools import combinations
from collections import defaultdict,Counter
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description="Regenerate the pinned native source-circuit discovery census.")
parser.add_argument('--repo',type=Path,help="Repository containing the original order-eight atlas; defaults to an enclosing repository.")
parser.add_argument('--out-dir',type=Path,help="Exact output directory; defaults to the original packet output location.")
args=parser.parse_args()
required=('05-knowledge/results/metagraph_order8_20261023/merged_nodes.tsv',
          '05-knowledge/results/metagraph_order8_20261023/owner_map.u16.gz',
          '05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json')
def has_inputs(root):return all((root/name).is_file() for name in required)
ROOT=args.repo.resolve() if args.repo else next((root for root in (HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents) if has_inputs(root)),None)
if ROOT is None or not has_inputs(ROOT):parser.error('Supply --repo with the original order-eight atlas and electrical receipt.')
H=args.out_dir.resolve() if args.out_dir else HERE
H.mkdir(parents=True,exist_ok=True)
base=ROOT/'05-knowledge/results/metagraph_order8_20261023';nodes=list(csv.DictReader((base/'merged_nodes.tsv').open(),delimiter='\t'))
D=json.load(open(ROOT/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json'))['order8'];X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4']
labels=[tuple(y[i]%2 for y in Y[:4])+(X[i]%128,) for i in range(len(nodes))];mods=[2,2,2,2,128]
owner=array('H');owner.frombytes(gzip.decompress((base/'owner_map.u16.gz').read_bytes()));
if sys.byteorder!='little':owner.byteswap()
merged={int(row[k]):int(row['merged_id']) for row in nodes for k in ('class_a','class_b')};chords=[(i,j) for i in range(8) for j in range(i+2,8)];pairs=list(combinations(range(8),2));full=(1<<21)-1
revmap=[chords.index((7-j,7-i)) for i,j in chords]
def resource(q):
 r=sum((q>>k&1)<<revmap[k] for k in range(21));return min(q,full^q,r,full^r)
records=[];start=time.time();sources=0
for s in range(len(nodes)):
 if nodes[s]['class_a']==nodes[s]['class_b']:continue
 canon=int(nodes[s]['canon_a']);A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(pairs):A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
 paths=[]
 def dfs(word,left):
  if not left:paths.append(word);return
  for v in left:
   if not word or A[word[-1]][v]:dfs(word+(v,),left-{v})
 dfs((),set(range(8)))
 if len(paths)*2!=int(nodes[s]['fibre']):continue
 sources+=1;fibres=defaultdict(list)
 for p in paths:
  q=sum(A[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(chords));t=merged[owner[full^q]]
  if t==s:continue
  r={'path':p,'target':t,'mask':q,'resource':resource(q)}
  for k in range(1,6):fibres[(p[0],p[-1],k,p[k],p[k+1],labels[t])].append(r)
 transitions=defaultdict(dict)
 for pins,rs in fibres.items():
  for r1,r2 in combinations(rs,2):
   t,u=r1['target'],r2['target']
   if t==u:continue
   if t>u:t,u=u,t;r1,r2=r2,r1
   cost=sum(a!=b for a,b in zip(r1['path'],r2['path']));key=(t,u);rkey=(r1['resource'],r2['resource'])
   cand={'source':s,'targets':[t,u],'cost':cost,'pins':pins[:5],'shared_clocks':[i for i in range(8) if r1['path'][i]==r2['path'][i]],'before':r1,'after':r2,'source_label':labels[s],'target_label':labels[t]}
   if rkey not in transitions[key] or pins[:5]<transitions[key][rkey]['pins']:transitions[key][rkey]=cand
 for targets,menu in transitions.items():
  menu=sorted(menu.values(),key=lambda r:(r['cost'],r['pins'],r['before']['path'],r['after']['path']));bestpair=None
  for i,a in enumerate(menu):
   for b in menu[i+1:]:
    if a['before']['resource']==b['before']['resource'] or a['after']['resource']==b['after']['resource']:continue
    pair=[a,b]
    if bestpair is None or sum(r['cost'] for r in pair)<sum(r['cost'] for r in bestpair):bestpair=pair
  records.append({'source':s,'targets':targets,'label':labels[s],'moves':menu,'best':menu[0],'bestpair':bestpair})
 if s%500==0:print(s,'menus',len(records),'sec',round(time.time()-start,2),flush=True)
print('ALL',sources,'eligible source-targetpairs',len(records),'capacity2',sum(r['bestpair'] is not None for r in records),'seconds',round(time.time()-start,2),flush=True)
(H/'all_menus.json').write_text(json.dumps(records,indent=2))
bytargets=defaultdict(list)
for r in records:bytargets[tuple(r['targets'])].append(r)
three=[];two=[]
for targets,rs in bytargets.items():
 for a,c in combinations(rs,2):
  if set(targets)&{a['source'],c['source']}:continue
  diff=tuple((u-v)%m for u,v,m in zip(a['label'],c['label'],mods))
  if any(diff) and all(2*x%m==0 for x,m in zip(diff,mods)) and a['bestpair'] and c['bestpair']:
   slots=a['bestpair']+c['bestpair'];two.append(dict(kind='two-source-double-rectangle',targets=targets,sources=[a['source'],c['source']],signs=[1,1,-1,-1],slots=slots,total_cost=sum(r['cost'] for r in slots)))
  for b in rs:
   if b['source'] in (a['source'],c['source']) or b['source'] in targets or not b['bestpair']:continue
   if all((u+w-2*v)%m==0 for u,w,v,m in zip(a['label'],c['label'],b['label'],mods)) and a['label']!=b['label']:
    slots=[a['best'],c['best']]+b['bestpair'];three.append(dict(kind='three-source-circuit',targets=targets,sources=[a['source'],b['source'],c['source']],signs=[1,1,-1,-1],slots=slots,total_cost=sum(r['cost'] for r in slots)))
for name,rs in [('three',three),('two',two)]:
 rs.sort(key=lambda r:(r['total_cost'],r['targets'],r['sources']));(H/(name+'_candidates.json')).write_text(json.dumps(rs,indent=2));print(name,'count',len(rs),flush=True)
 for r in rs[:10]:print(r['targets'],r['sources'],r['total_cost'],flush=True)
 if rs:(H/(name+'_candidate.json')).write_text(json.dumps(rs[0],indent=2))
