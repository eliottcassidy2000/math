#!/usr/bin/env python3
"""Deterministic physical records and complete small named-response controls."""
import collections,csv,gzip,io,itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
B=P.parent/'session189-c3-frontier'
TYPES={'cycle_source':([2,4,1,7],[3,18,36,24]),'cycle_sink':([10,12,9,0],[3,18,36,24]),'strong_four':([6,4,8,3],[5,22,36,24])}

def load(p):return json.loads(p.read_text())
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def expand(m,j):return sum(7<<(3*i) for i in range(9) if m>>i&1)+((1<<27) if m>>j&1 else 0)
def inflate(q,j,typ):
 cells=[list(range(3*i,3*i+3))+([27] if i==j else []) for i in range(9)]
 r=[0]*28
 for i,c in enumerate(cells):
  inner=TYPES[typ][0] if i==j else [2,4,1]
  cross=sum(1<<v for k,d in enumerate(cells) if q[i]>>k&1 for v in d)
  for x,v in enumerate(c):r[v]=cross+sum(1<<c[y] for y in range(len(c)) if inner[x]>>y&1)
 return r

def aug(r,a,c):
 n=len(r);f=(1<<n)-1
 return [z+((a>>i&1)<<n)+((c>>i&1)<<(n+1)) for i,z in enumerate(r)]+[(f^a)|(1<<(n+1)),f^c]

def hamilton(r):
 n=len(r);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i,v in enumerate(dp[s]):
   todo=r[i]&~s
   while v and todo:
    b=todo&-todo;todo-=b;j=b.bit_length()-1;dp[s|b][j]+=v
 return sum(dp[-1])

def named_control(typ):
 r=TYPES[typ][0];n=4;op=[15^z^(1<<i) for i,z in enumerate(r)];records=[]
 for t in (r,op):
  values=[hamilton(aug(t,m&15,m>>4)) for m in range(256)];coeff=values[:]
  for j in range(8):
   for m in range(256):
    if m>>j&1:coeff[m]-=coeff[m^(1<<j)]
  records.append(dict(rows=t,values=values,coefficients=coeff))
 return dict(type=typ,variable_order=['a0','a1','a2','a3','c0','c1','c2','c3'],source=records[0],converse=records[1])

def main():
 candidates=load(B/'candidates.json')['candidates'];byid={r['node']:r for r in candidates};initial=tab(P/'residues.tsv');extra={(int(r['node']),int(r['position'])):r for r in tab(P/'selected_residues.tsv')};qs=[(r['a'],r['c']) for r in load(B/'battery.json')]+[(r['a'],r['c']) for r in load(P/'selected_extras.json')]
 (P/'battery.json').write_text(json.dumps([dict(index=i,a=a,c=c) for i,(a,c) in enumerate(qs)],indent=2)+'\n')
 witnesses=[];hist=collections.Counter()
 for r in initial:
  node,j=int(r['node']),int(r['position']);hits=[(k,int(r[f'S{k}']),int(r[f'T{k}'])) for k in range(12) if r[f'S{k}']!=r[f'T{k}']]
  if not hits:
   x=extra[node,j];hits=[(12+k,int(x[f'S{k}']),int(x[f'T{k}'])) for k in range(2) if x[f'S{k}']!=x[f'T{k}']]
  if not hits:raise RuntimeError('unseparated placement')
  k,s,t=hits[0];a,c=qs[k];hist[k]+=1;witnesses.append(dict(node=node,position=j,query_index=k,a=a,c=c,physical_a=expand(a,j),physical_c=expand(c,j),normalized_source=s,normalized_target=t))
 (P/'witnesses.tsv').write_text('node\tposition\tquery_index\ta\tc\tphysical_a\tphysical_c\tnormalized_source\tnormalized_target\n'+'\n'.join('\t'.join(map(str,w.values())) for w in witnesses)+'\n')
 with (P/'physical_witnesses.jsonl.gz').open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as out:
   for w in witnesses:
    q=byid[w['node']]['rows'];op=[511^z^(1<<i) for i,z in enumerate(q)]
    for typ,(_,profile) in TYPES.items():
     factor=profile[0];record=dict(w,type=typ,profile=profile,source_residue=w['normalized_source']*factor%8,target_residue=w['normalized_target']*factor%8,augmented_source_rows=aug(inflate(q,w['position'],typ),w['physical_a'],w['physical_c']),augmented_target_rows=aug(inflate(op,w['position'],typ),w['physical_a'],w['physical_c']))
     out.write((json.dumps(record,separators=(',',':'))+'\n').encode())
 (P/'named_controls.json').write_text(json.dumps([named_control(t) for t in TYPES],indent=2)+'\n')
 print(json.dumps(dict(placements=len(witnesses),physical_witnesses=3*len(witnesses),physical_graphs=6*len(witnesses),first_query_histogram=[hist[k] for k in range(14)]),indent=2))
if __name__=='__main__':main()
