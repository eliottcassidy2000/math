from pathlib import Path
from collections import defaultdict
import json,math,sys,hashlib
root=Path(__file__).parent
if '--regenerate' in sys.argv:
 orders={};seeds=defaultdict(set);totalstates=0;contentone=0;sources=[]
 for n in range(3,10):
  p=root/f'SATURATED{n}.jsonl';rows=list(map(json.loads,p.read_text().splitlines()));sources.append(dict(order=n,rows=len(rows),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  for row in rows:
   h,q,a,g=row[:4];totalstates+=1
   if math.gcd(h,a,g)>1:continue
   contentone+=1
   if h not in orders:
    order=1
    while pow(2,order,h)!=1:order+=1
    orders[h]=order
   if orders[h]<12:seeds[h,q].add((a,g,a%2))
 data=dict(saturated_signature_rows=totalstates,content_one_rows=contentone,sources=sources,groups=[dict(h=h,q=q,order2=orders[h],seeds=sorted(vs)) for (h,q),vs in sorted(seeds.items())])
 (root/'MINIMALITY_SEEDS.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
else:
 data=json.loads((root/'MINIMALITY_SEEDS.json').read_text());totalstates=data['saturated_signature_rows'];contentone=data['content_one_rows'];orders={r['h']:r['order2'] for r in data['groups']};seeds={(r['h'],r['q']):set(map(tuple,r['seeds'])) for r in data['groups']}
results=[];alltested=0
for (h,q),starts in sorted(seeds.items()):
 lengths=list(range(orders[h],12,orders[h]));cache={};returns=0;gains=set();paths=0
 def edges(node):
  if node in cache:return cache[node]
  a,g,p=node;x=g%2;out=[]
  for y in (0,1):
   aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2
   if -h<aa<q and -h<gg<q:out.append(((aa,gg,y),a*y+h*p*y-x*gg-h*x*(gg%2),x,y))
  cache[node]=out;return out
 for start in sorted(starts):
  for length in lengths:
   stack=[(start,0,0,[])]
   while stack:
    v,depth,gain,path=stack.pop();paths+=1
    if depth==length:
     if v==start:
      returns+=1;gains.add(gain)
      if abs(gain)==2:raise RuntimeError(('short unit',h,q,start,path))
     continue
    for w,wt,x,y in edges(v):stack.append((w,depth+1,gain+wt,path+[(x,y)]))
 results.append(dict(h=h,q=q,order2=orders[h],physical_AG_seeds=len(starts),lengths=lengths,partial_walks=paths,closed_walks=returns,gains=sorted(gains)))
 alltested+=paths
print(json.dumps(dict(status='PASS no unit loop shorter than12 in complete core-order<=9 bank',saturated_signature_rows=totalstates,content_one_rows=contentone,partial_walks=alltested,groups=results),indent=2))
