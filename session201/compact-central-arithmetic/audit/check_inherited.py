"""Direct DFS recheck of the inherited 27 adjacency-list repairs."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
if not(ROOT/'MENUS.json').is_file():ROOT=ROOT/'session201-compact-central'
receipts=[]
for r in json.loads((ROOT/'dependencies/REPAIR_GRAPHS.json').read_text()):
 rows=r['rows'];n=len(rows);full=(1<<n)-1;paths=cycles=nodes=0
 if any(x<0 or x>full or x>>i&1 for i,x in enumerate(rows)):raise RuntimeError('invalid row')
 if any(((rows[i]>>j)&1)+((rows[j]>>i)&1)!=1 for i in range(n)for j in range(i)):raise RuntimeError('not tournament')
 def visit(first,last,left):
  global paths,cycles,nodes
  nodes+=1
  if not left:
   paths+=1
   if first==0 and rows[last]&1:cycles+=1
   return
  nxt=rows[last]&left
  while nxt:
   bit=nxt&-nxt;nxt-=bit;visit(first,bit.bit_length()-1,left^bit)
 for i in range(n):visit(i,i,full^(1<<i))
 if(paths,cycles)!=(r['H'],r['C']):raise RuntimeError('count mismatch')
 if n!=r['order']or n>r['H'].bit_length()-1+6:raise RuntimeError('size bound')
 receipts.append(dict(target=paths,cycles=cycles,vertices=n,partial_paths=nodes))
if len(receipts)!=27:raise RuntimeError('count of inherited repairs')
raw=json.loads((ROOT/'RAW-REPAIRS.json').read_text());repair=json.loads((ROOT/'REPAIRS.json').read_text())
for cert,(n,i,u,v,k)in zip(repair['repairs'],raw['repairs']):
 if(n,u,v,k)!=(cert['target'],cert['left_index'],cert['right_index'],len(cert['center'])):raise RuntimeError('raw profile metadata')
extra=sorted(n for n,i,u,v,k in raw['repairs']if n<=17383 and k>6)
if extra!=[1811,1877,1897,6151,6721,17383]:raise RuntimeError('six length-six repairs')
print(json.dumps(dict(status='PASS',inherited_graphs=receipts,length_six_additional_gaps=extra,
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
