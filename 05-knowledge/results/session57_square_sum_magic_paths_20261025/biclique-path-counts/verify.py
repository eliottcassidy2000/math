import itertools as I, collections as C, json
V=range(4); S=[(0,0),(1,0),(0,1),(1,1)]
edges=list(I.product(V,V)); labels=[(S[j][0]-S[i][0],S[j][1]-S[i][1]) for i,j in edges]; labels_unique=sorted(set(labels)); li={v:i for i,v in enumerate(labels_unique)}
def emask(seq):
 return sum(1<<(a*4+b) for a,b in seq)
paths=[]; cycles={}
# Each unoriented Hamilton path has exactly one orientation beginning in left class.
for a in I.permutations(V):
 for b in I.permutations(V):
  es=[(a[k],b[k]) for k in V]+[(a[k+1],b[k]) for k in range(3)]
  m=emask(es); paths.append({'vertices':sum(([f'L{a[k]}',f'R{b[k]}'] for k in V),[]),'mask':m,'label_mask':sum(1<<z for z in set(li[labels[e]] for e in range(16) if m>>e&1))})
  cm=m | 1<<(a[0]*4+b[3]); cycles[cm]=1
assert len({p['mask'] for p in paths})==576
cycs=list(cycles)
counts={}
for name,ms in [('paths',[p['mask'] for p in paths]),('cycles',cycs)]:
 ec=[sum(m>>e&1 for m in ms) for e in range(16)]
 pc=C.defaultdict(set)
 for e,f in I.combinations(range(16),2):
  typ='adjacent' if edges[e][0]==edges[f][0] or edges[e][1]==edges[f][1] else 'disjoint'
  pc[typ].add(sum((m>>e&1) and(m>>f&1) for m in ms))
 counts[name]={'number':len(ms),'edge_incidence':sorted(set(ec)),'edge_pair_incidence':{k:sorted(v) for k,v in pc.items()}}
labelcovers=[]
for i,j in I.combinations(range(576),2):
 if paths[i]['label_mask']|paths[j]['label_mask']==511: labelcovers.append((i,j))
# Integer bitsets index paths containing a given edge. Count edge-cover triples exactly once.
inc=[sum(1<<i for i,p in enumerate(paths) if p['mask']>>e&1) for e in range(16)]
allp=(1<<576)-1; triples=0; witness=None
for i,j in I.combinations(range(576),2):
 missing=65535^(paths[i]['mask']|paths[j]['mask']); candidates=allp ^ ((1<<(j+1))-1)
 while missing and candidates:
  bit=missing & -missing; candidates &= inc[bit.bit_length()-1]; missing^=bit
 triples+=candidates.bit_count()
 if candidates and witness is None: witness=(i,j,(candidates&-candidates).bit_length()-1)
result={'offsets':S,'labels_by_edge':labels,'label_multiplicities':{str(k):v for k,v in C.Counter(labels).items()},'counts':counts,'path_label_count_distribution':dict(sorted(C.Counter(p['label_mask'].bit_count() for p in paths).items())),'min_label_cover':2,'number_label_cover_pairs':len(labelcovers),'label_cover_witness':[paths[i] for i in labelcovers[0]],'min_edge_cover':3,'number_edge_cover_triples':triples,'edge_cover_witness':[paths[i] for i in witness], 'cycle_cover_pairs':sum(a|b==65535 for a,b in I.combinations(cycs,2))}
assert labelcovers and witness
assert all(p['label_mask'].bit_count()<9 for p in paths)
assert all((paths[i]['mask']|paths[j]['mask']).bit_count()<16 for i,j in I.combinations(range(576),2))
open(__import__('pathlib').Path(__file__).with_name('results.json'),'w').write(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
assert counts['paths']=={'number':576,'edge_incidence':[252],'edge_pair_incidence':{'adjacent':[72],'disjoint':[120]}}
assert counts['cycles']=={'number':72,'edge_incidence':[36],'edge_pair_incidence':{'adjacent':[12],'disjoint':[20]}}
assert len(labelcovers)==20624 and triples==675264
union=0; intersection=65535
for p in paths: union|=p['mask']; intersection&=p['mask']
assert union==65535 and intersection==0
fixed=[]
for a in I.permutations(V):
 for b in I.permutations(V):
  for swap in [False,True]:
   target=[(b[j]*4+a[i]) if swap else (a[i]*4+b[j]) for i,j in edges]
   if all(labels[e]==labels[target[e]] for e in range(16)): fixed.append((a,b,swap))
assert fixed==[((0,1,2,3),(0,1,2,3),False),((3,2,1,0),(3,2,1,0),True)]
# A 3x3 magic arrangement of all nine formal coefficient pairs.
magic=[[(1,0),(-1,-1),(0,1)],[(-1,1),(0,0),(1,-1)],[(0,-1),(1,1),(-1,0)]]
assert set(sum(magic,[]))==set(labels_unique)
lines=magic+list(zip(*magic))+[[magic[i][i] for i in V if i<3],[magic[i][2-i] for i in V if i<3]]
assert all(tuple(map(sum,zip(*line)))==(0,0) for line in lines)
print('All exact count, cover, fixed-label symmetry, and magic-line assertions passed.')
uniondist=C.Counter((paths[i]['mask']|paths[j]['mask']).bit_count() for i,j in labelcovers)
minunion=min(uniondist)
best=next((i,j) for i,j in labelcovers if (paths[i]['mask']|paths[j]['mask']).bit_count()==minunion)
extra={'label_cover_pair_union_edge_distribution':dict(sorted(uniondist.items())),'minimum_union_edges':minunion,'maximum_shared_edges':14-minunion,'minimal_union_witness':[paths[i] for i in best]}
open(__import__('pathlib').Path(__file__).with_name('minimal_union.json'),'w').write(json.dumps(extra,indent=2))
print(json.dumps(extra,indent=2))
assert dict(sorted(uniondist.items()))=={9:204,10:1852,11:5600,12:7720,13:3760,14:1488}
assert minunion==9
best_union=paths[best[0]]['mask']|paths[best[1]]['mask']
assert len(set(labels[e] for e in range(16) if best_union>>e&1))==9
print('Minimal nine-edge union certificate verified.')
