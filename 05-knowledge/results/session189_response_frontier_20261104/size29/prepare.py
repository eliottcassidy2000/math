#!/usr/bin/env python3
"""Retain exact size29 root queries, normalized residues, and physical realizations."""
import collections,csv,gzip,itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent;B=P.parent/'session189-c3-frontier';S=P.parent/'session189-size28'
PAIRS=list(itertools.combinations(range(9),2))
FIVE_A=[4,1,2,7,15];FIVE_B=[16,1,3,7,14];FOUR=[2,4,1,7]
def load(p):return json.loads(p.read_text())
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def configs(i):
 if i<10:return 'five_6',[i-1]
 if i<46:return 'two_four',list(PAIRS[i-10])
 return 'five_2',[i-46]
def cells_for(typ,pos):
 cells=[list(range(3*i,3*i+3)) for i in range(9)]
 if typ=='two_four':cells[pos[0]].append(27);cells[pos[1]].append(28)
 else:cells[pos[0]].extend([27,28])
 return cells
def inflate(q,typ,pos):
 cells=cells_for(typ,pos);rows=[0]*29
 for i,cell in enumerate(cells):
  inner=FOUR if typ=='two_four' and i in pos else FIVE_A if typ=='five_6' and i in pos else FIVE_B if typ=='five_2' and i in pos else [2,4,1]
  cross=sum(1<<v for j,other in enumerate(cells) if q[i]>>j&1 for v in other)
  for a,v in enumerate(cell):rows[v]=cross|sum(1<<cell[b] for b in range(len(cell)) if inner[a]>>b&1)
 return rows
def expand(mask,typ,pos):return sum(1<<v for i,cell in enumerate(cells_for(typ,pos)) if mask>>i&1 for v in cell)
def aug(r,a,c):
 n=len(r);f=(1<<n)-1
 return [z|((a>>i&1)<<n)|((c>>i&1)<<(n+1)) for i,z in enumerate(r)]+[(f^a)|(1<<(n+1)),f^c]
def derive():
 cs=load(B/'candidates.json')['candidates'];rows=tab(P/'root_residues.tsv');extra=tab(P/'root_selected_residues.tsv');qs=[(q['a'],q['c']) for q in load(S/'battery.json')]+[tuple(q) for q in load(P/'root_selected_extras.json')]
 by=collections.defaultdict(list)
 for r in rows:by[int(r['node'])].append(r)
 ex=collections.defaultdict(list)
 for r in extra:ex[int(r['node'])].append(r)
 witnesses=[];miss=[]
 for c in cs:
  node=c['node'];rs=by[node]
  for i in range(1,55):
   v=i if i<46 else 0
   hit=next(((k,r) for k,r in enumerate(rs) if r[f'S{v}']!=r[f'T{v}']),None)
   if hit is None:
    miss.append(dict(node=node,configuration=i));hit=next(((14+k,r) for k,r in enumerate(ex[node]) if r[f'S{v}']!=r[f'T{v}']),None)
   if hit is None:raise RuntimeError(('unseparated',node,i))
   k,r=hit;typ,pos=configs(i);a,c=qs[k];factor=1 if typ=='five_2' else 3
   witnesses.append(dict(node=node,configuration=i,type=typ,positions=pos,query_index=k,a=a,c=c,normalized_source=int(r[f'S{v}']),normalized_target=int(r[f'T{v}']),physical_factor=factor,physical_a=expand(a,typ,pos),physical_c=expand(c,typ,pos)))
 return qs,witnesses,miss

def physical(w,q):
 typ,pos=w['type'],w['positions'];op=[511^z^(1<<i) for i,z in enumerate(q)];fac=w['physical_factor']
 return dict(w,source_residue=fac*w['normalized_source']%8,target_residue=fac*w['normalized_target']%8,augmented_source_rows=aug(inflate(q,typ,pos),w['physical_a'],w['physical_c']),augmented_target_rows=aug(inflate(op,typ,pos),w['physical_a'],w['physical_c']))
def main():
 cs=load(B/'candidates.json')['candidates'];byid={r['node']:r for r in cs};qs,witnesses,miss=derive()
 (P/'root_battery.json').write_text(json.dumps([dict(index=i,a=a,c=c) for i,(a,c) in enumerate(qs)],indent=2)+'\n')
 with (P/'root_witnesses.tsv').open('w') as f:
  keys=list(witnesses[0]);f.write('\t'.join(keys)+'\n')
  for w in witnesses:f.write('\t'.join(json.dumps(w[k],separators=(',',':')) if isinstance(w[k],list) else str(w[k]) for k in keys)+'\n')
 with (P/'root_physical_witnesses.jsonl.gz').open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as out:
   for w in witnesses:out.write((json.dumps(physical(w,byid[w['node']]['rows']),separators=(',',':'))+'\n').encode())
 summary=dict(quotients=len(cs),exceptional_five_placements=9*len(cs),two_four_placements=36*len(cs),inherited_five_placements=9*len(cs),initial_queries=14,initial_misses=len(miss),extra_queries=qs[14:],physical_witness_pairs=len(witnesses),physical_graphs=2*len(witnesses),endpoint_union_queries=len(set(qs+[(1<<i,0) for i in range(9)])),histogram=dict(collections.Counter(w['query_index'] for w in witnesses)))
 (P/'root_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
