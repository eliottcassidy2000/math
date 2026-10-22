from pathlib import Path
from itertools import combinations,permutations
from collections import Counter,defaultdict
import json
BASE=Path(__file__).resolve().parent.parent/'session29-profiles'
DATA=json.loads((BASE/'seamed-results.json').read_text())
P={n:list(combinations(range(n),2)) for n in (4,5,10)}
def require(b,msg):
 if not b:raise RuntimeError(msg)
def adj(m,n):
 a=[[0]*n for _ in range(n)]
 for z,(i,j) in enumerate(P[n]):a[i][j]=(m>>z)&1;a[j][i]=1-a[i][j]
 return a
def mask(a):return sum(a[i][j]<<z for z,(i,j) in enumerate(P[len(a)]))
def canon_root(a,root):
 n=len(a);ans=None;w=None;autos=0
 for p0 in permutations([i for i in range(n) if i!=root]):
  p=tuple(p0)+(root,);m=mask([[a[i][j] for j in p] for i in p])
  if ans is None or m<ans:ans=m;w=p;autos=1
  elif m==ans:autos+=1
 inv=[w.index(i) for i in range(n)]
 return ans,inv,autos
def cp(t):return t ^ (1023^sum(1<<z for z,(i,j) in enumerate(P[5]) if j==i+1))
def op(a):return [[a[j][i] for j in range(len(a))] for i in range(len(a))]
def canonical_target(q,r):
 x=adj(q,5);y=op(adj(r,5));a=[[0]*10 for _ in range(10)]
 for i in range(10):
  for j in range(10):
   if i//5==j//5:a[i][j]=(x if i<5 else y)[i%5][j%5]
   elif i>=5:a[i][j]=1
 a[4][9]=1;a[9][4]=0
 return mask(a)
def cuts(a):
 n=len(a);ans=[]
 for s in range(1,(1<<n)-1):
  es=[(i,j) for i in range(n) if s>>i&1 for j in range(n) if not(s>>j&1) and a[i][j]]
  if len(es)==1:ans.append({'vertices':[i for i in range(n) if s>>i&1],'arc':list(es[0])})
 return ans
def strong(a):
 n=len(a)
 for rev in (False,True):
  seen={0};todo=[0]
  while todo:
   i=todo.pop()
   for j in range(n):
    if j not in seen and (a[j][i] if rev else a[i][j]):seen.add(j);todo.append(j)
  if len(seen)<n:return False
 return True
rooted={};entry={};classes=defaultdict(lambda:{'A':[],'B':[]})
for seed,tiles in DATA['tile_alphabet'].items():
 typ='A' if int(seed)==40 else 'B'
 for t in tiles:
  x=adj(cp(t),5);q,w,aut=canon_root(x,4);r,v,aut2=canon_root(op(x),0)
  rooted[t]=(q,w,aut);entry[t]=(r,v,aut2);classes[q][typ].append(t)
require(set(x[0] for x in entry.values())==set(classes),'same reflected root alphabet')
ordinary=defaultdict(list);merged=defaultdict(list);outrows=[];bridgehist=Counter()
for z,row in enumerate(DATA['rows']):
 t,u=row['tiles'];q,wx,ax=rooted[t];r,wy,ay=entry[u];a=adj(row['target'],10)
 ordinary[(q,r)].append(z);key=tuple(sorted((q,r)));merged[key].append(z)
 conv=q>r;w=([x+5 for x in wx]+wy) if conv else (wx+[y+5 for y in wy]);cm=canonical_target(*key);ca=adj(cm,10)
 require(all((a[j][i] if conv else a[i][j])==ca[w[i]][w[j]] for i in range(10) for j in range(10)),'explicit iso/converse witness')
 cc=cuts(a);balanced=[c for c in cc if len(c['vertices'])==5]
 require(len(balanced)==1 and balanced[0]['vertices']==list(range(5)) and balanced[0]['arc']==[4,5],'unique balanced seam')
 require(strong(a),'strong target')
 bridgehist[len(set(tuple(c['arc']) for c in cc))]+=1
 source='AA' if row['word']==[40,40] else ('BB' if row['word']==[41,41] else 'AB/BA')
 outrows.append({'row':z,'source':source,'tiles':[t,u],'target':row['target'],'H':row['H'],'blue':row['blue'],'ordinary_key':[q,r],'merged_key':list(key),'converse_used':conv,'old_to_canonical_vertex':w,'target_aut_order':ax*ay,'one_arc_cuts':cc})
table=[]
for key,ids in sorted(merged.items()):
 q,r=key;a,b=len(classes[q]['A']),len(classes[q]['B']);c,d=len(classes[r]['A']),len(classes[r]['B'])
 expect={'AA':a*c,'AB/BA':a*d+b*c,'BB':b*d}
 if q!=r:expect={k:2*v for k,v in expect.items()}
 counts={s:{'blue':0,'black':0} for s in expect}
 hs=set()
 for z in ids:
  row=outrows[z];counts[row['source']]['blue' if row['blue'] else 'black']+=1;hs.add(row['H'])
 require(all(sum(counts[s].values())==expect[s] for s in expect),'symmetric-square incidence')
 require(len(hs)==1,'isomorphism H invariance')
 if q==r:
  require(counts['AA']['blue']==a and counts['BB']['blue']==b and counts['AB/BA']['blue']==0,'diagonal blue formula')
 else:require(all(x['blue']==0 for x in counts.values()),'offdiagonal black')
 table.append({'root_pair':list(key),'self_converse':q==r,'H':hs.pop(),'canonical_target_mask':canonical_target(*key),'incidences':counts,'total':len(ids)})
# Independent exhaustive hostile test of the cut lemma at orders4 and5.
small=0
for n in (4,5):
 for m in range(1<<len(P[n])):
  a=adj(m,n)
  if not strong(a):continue
  small+=1;c=Counter(len(x['vertices']) for x in cuts(a))
  require(all(c[k]<=1 for k in range(2,n-1)),'all small strong tournament cut uniqueness')
res={'root_classes':{str(q):dict(v) for q,v in sorted(classes.items())},'root_class_count':len(classes),'ordinary_target_classes':len(ordinary),'merged_target_classes':len(merged),'self_converse_target_classes':sum(x['self_converse'] for x in table),'total_edges':len(outrows),'blue_edges':sum(x['blue'] for x in outrows),'strong_bridge_count_histogram':dict(sorted(bridgehist.items())),'exhaustive_small_strong_tournaments':small,'incidence_table':table,'witness_rows':outrows}
Path(__file__).with_name('wiring-results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k not in ('incidence_table','witness_rows')},indent=2))
