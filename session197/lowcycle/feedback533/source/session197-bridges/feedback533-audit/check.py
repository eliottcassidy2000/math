from pathlib import Path
from itertools import combinations
import json,hashlib
HERE=Path(__file__).resolve().parent

def need(v,s):
 if not v:raise RuntimeError(s)
def walk_ie(a):
 n=len(a);total=0
 for S in range(1,1<<n):
  v=[(S>>u)&1 for u in range(n)]
  for _ in range(n-1):
   w=[0]*n
   for u,x in enumerate(v):
    if not x:continue
    U=a[u]&S
    while U:
     b=U&-U;U-=b;w[b.bit_length()-1]+=x
   v=w
  total+=(-1)**(n-S.bit_count())*sum(v)
 return total

def cycles(a):
 n=len(a);full=(1<<n)-1;count=0
 def rec(last,S):
  nonlocal count
  if S==full:count+=bool(a[last]&1);return
  U=a[last]&~S
  while U:
   b=U&-U;U-=b;rec(b.bit_length()-1,S|b)
 rec(0,1);return count

def fvs(a):
 n=len(a);best=0;witness=[]
 for S in range(1<<n):
  vs=[v for v in range(n) if S>>v&1];d=sorted((a[v]&S).bit_count() for v in vs)
  if d==list(range(len(vs))) and len(vs)>best:best=len(vs);witness=[v for v in range(n) if not S>>v&1]
 return n-best,witness

def iso(a,b):
 n=len(a);da=[x.bit_count() for x in a];db=[x.bit_count() for x in b];p=[-1]*n;found=[]
 def rec(k,used):
  if k==n:found.append(p[:]);return
  candidates=[]
  for u in range(n):
   if p[u]>=0:continue
   opts=[v for v in range(n) if not used>>v&1 and da[u]==db[v] and all(p[x]<0 or bool(a[u]>>x&1)==bool(b[v]>>p[x]&1) for x in range(n))]
   candidates.append((len(opts),u,opts))
  _,u,opts=min(candidates)
  for v in opts:p[u]=v;rec(k+1,used|1<<v);p[u]=-1
 rec(0,0);return found
cert=json.loads(Path('session197-lowcycle/feedback533_certificate.json').read_text());ind=[list(map(int,s.split(',')))[2:] for s in (HERE/'independent-hits.csv').read_text().splitlines()];src=[x['rows'] for x in cert['cores']]
need(len(ind)==len(src)==6,'six records');maps=[];seen=set();controls=[]
for a in ind:
 matches=[(j,iso(a,b)) for j,b in enumerate(src)];matches=[(j,m) for j,m in matches if m]
 need(len(matches)==1,'unique source isomorphism class');j,m=matches[0];need(j not in seen,'class repeated');seen.add(j);maps.append({'source_index':j,'mapping':m[0],'isomorphisms':len(m)})
 need(walk_ie(a)==533 and cycles(a)==2,'independent actual counts');f,deleted=fvs(a);need(f==3,'exhaustive feedback number');controls.append({'rows':a,'H':533,'C':2,'feedback':f,'deletion_witness':deleted})
# Check all listed certificates, including every pair deletion.
triangles=0
for item in cert['cores']:
 a=item['rows'];pairs=set()
 for w in item['pair_deletion_cycle_witnesses']:
  removed=tuple(w['deleted']);t=w['cycle'];need(len(t)==3 and len(set(t))==3 and not set(t)&set(removed),'triangle labels');need(all(a[t[i]]>>t[(i+1)%3]&1 for i in range(3)),'triangle edges');pairs.add(removed);triangles+=1
 need(pairs==set(combinations(range(11),2)),'every pair exactly represented')
 w=item['three_deletion_transitive_witness'];vs=w['transitive_order'];need(set(vs)==set(range(11))-set(w['deleted']) and len(vs)==8,'upper labels');need(all(a[vs[i]]>>vs[j]&1 for i in range(8) for j in range(i+1,8)),'upper order')
word=cert['changing_two_port_witness'];birth=word['birth_order'];moves=word['moves'];a=src[0];ports=birth[:2];need(a[ports[0]]>>ports[1]&1,'initial forward edge');arcs=1
for k,(mask,drop) in enumerate(moves,2):
 v=birth[k]
 for u in birth[:k]:
  wins=u in ports and bool(mask>>(ports.index(u))&1);need(bool(a[v]>>u&1)==wins,'physical birth orientation');arcs+=1
 ports=[u for j,u in enumerate(ports+[v]) if j!=drop]
need(ports==[1,0] and arcs==55,'complete physical word');need(533**3<5**12,'strict cutoff')
out={'status':'PASS','independent_bank_sha256':hashlib.sha256(Path('session197-bridges/c3-independent-records.csv').read_bytes()).hexdigest(),'six_complete_isomorphism_matches':maps,'all_subset_feedback_checks':6*2**11,'pair_triangle_certificates':triangles,'physical_word_arcs':arcs,'final_ports':ports,'controls':controls}
print(json.dumps(out,sort_keys=True,indent=2))
