from pathlib import Path
import json,itertools,hashlib
P=Path(__file__).resolve().parent;G=P.parent/'gadgets'
def need(x):
 if not x:raise RuntimeError('audit check failed')
def counts(rows,root):
 h=e=0
 for q in itertools.permutations(range(len(rows))):
  bad=[i for i in range(len(q)-1) if not(rows[q[i]]>>q[i+1]&1)]
  h+=not bad;e+=not bad or bad==[q.index(root)]
 return h,e
checks=0
for n in range(1,4):
 edges=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(edges)):
  b=[0]*n
  for k,(u,v) in enumerate(edges):
   if bits>>k&1:b[u]|=1<<v
   else:b[v]|=1<<u
  for W in range(1,1<<n):
   answers=[]
   for triangle in [False,True]:
    k=3 if triangle else 2;root=n+k;rows=[r for r in b]+[0]*(k+1)
    for i in range(n):
     if W>>i&1:rows[i]|=1<<root
     else:rows[root]|=1<<i
    for i in range(k):rows[n+i]|=(1<<n)-1;rows[root]|=1<<(n+i)
    for i in range(k):
     for j in range(i+1,k):
      if triangle and (i,j)==(0,2):rows[n+j]|=1<<(n+i)
      else:rows[n+i]|=1<<(n+j)
    answers.append(counts(rows,root))
   need(answers[1]==tuple(3*x for x in answers[0]));checks+=1
src=json.loads((G/'source.json').read_text());tri=json.loads((G/'triple.json').read_text());names=tri['names'];rows=tri['rows'];n=len(rows)
need(names[:-1]==src['names']);need(all(r&((1<<(n-1))-1)==s for r,s in zip(rows[:-1],src['rows'])))
a,b,x=map(names.index,[50,51,55]);root=names.index(39)
need(rows[a]>>b&1 and rows[b]>>x&1 and rows[x]>>a&1)
for v in range(n):
 if v not in [a,b,x]:need(len({(rows[v]>>z)&1 for z in [a,b,x]})==1)
raw=(G/'triple_deltas.txt').read_text().splitlines();need(raw[0]=='H 1537275');d=[list(map(int,s.split())) for s in raw[1:]]
need(len(d)==153);need({(u,v) for u,v,*_ in d}=={(u,v) for u in range(n) for v in range(n) if rows[u]>>v&1});need(len({(u,v) for u,v,*_ in d})==153)
need(all(delta==back-front for u,v,delta,front,back in d));need(all(delta not in [2,4] for u,v,delta,*_ in d))
zeros=[(names[u],names[v]) for u,v,delta,*_ in d if not delta];small=min((abs(delta),names[u],names[v],delta) for u,v,delta,*_ in d if delta)
need(zeros==[(6,8)]);need(small==(28062,22,2,-28062))
out={'small_literal_tripling_checks_including_nonstrong_suffixes':checks,'source_induced_names_and_edges_preserved':True,'homogeneous_C3_verified':True,'complete_distinct_edge_responses':len(d),'zero_responses':zeros,'minimum_nonzero_absolute_signed_response':small,'input_sha256':{f:hashlib.sha256((G/f).read_bytes()).hexdigest() for f in ['source.json','triple.json','triple_deltas.txt']}}
(P/'independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
