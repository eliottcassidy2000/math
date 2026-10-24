#!/usr/bin/env python3
"""Independent sparse bridge/coordinate audit using prior native domain records."""
import argparse,hashlib,json
from pathlib import Path
from collections import defaultdict,Counter,deque
from itertools import combinations
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));p.add_argument('--output',type=Path);args=p.parse_args()
base=args.repo/'05-knowledge/results';incoming=base/'sixslot570_20261024';rows=sorted(json.loads((base/'native_packet_decomposition_20261023/native/domains.json').read_text())['sources'],key=lambda r:r['source']);data=json.loads((incoming/'bridges.json').read_text());checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
x=[r['code'][0]for r in rows];y=[sum(z<<j for j,z in enumerate(r['code'][1:]))for r in rows];w=[sum(a!=b for a,b in zip(r['old_before'],r['old_after']))for r in rows]
def f(xx,yy):return ((xx&1)^((yy>>3)&1))|((yy&7)<<1)
def inv(xx,pp):return (pp>>1)|(((pp&1)^(xx&1))<<3)
for xx in range(128):
 for yy in range(16):
  need(inv(xx,f(xx,yy))==yy,'invertible original/incoming parity-coordinate transform')
  for dx,dy in [(1,0),(0,1),(0,2),(0,4),(0,8)]:need(f((xx+dx)%128,yy^dy)==f(xx,yy)^f(dx,dy),'transform respects all additive generators')
lab=[tuple(map(int,z.split()))for z in (incoming/'labels.txt').read_text().splitlines()]
need(lab==[(f(xx,yy),xx)for xx,yy in zip(x,y)],'all actual labels agree under group automorphism')
need(data['sources']==[r['source']for r in rows],'source order agrees')
three=sum(1<<i for i,v in enumerate(w)if v==3);four=sum(1<<i for i,v in enumerate(w)if v==4)
def work(m):return 2*m.bit_count()+(m&three).bit_count()+2*(m&four).bit_count()
def q(m):
 xx=yy=0
 while m:
  bit=m&-m;i=bit.bit_length()-1;m^=bit;xx+=x[i];yy^=y[i]
 return xx%128,yy
class DSU:
 def __init__(self,n):self.p=list(range(n))
 def find(self,i):
  while self.p[i]!=i:self.p[i]=self.p[self.p[i]];i=self.p[i]
  return i
 def join(self,a,b):
  a,b=self.find(a),self.find(b)
  if a==b:return False
  self.p[a]=b;return True
banks={};where={};comps={};groups={};report=[]
for k in (4,5):
 g=defaultdict(list)
 for ii in combinations(range(28),k):
  m=sum(1<<i for i in ii);g[q(m)].append(m)
 cc={};ix={};counts=Counter();excess=0
 for key,ss in g.items():
  d=DSU(len(ss))
  for i,a in enumerate(ss):
   for j,b in enumerate(ss[:i]):
    if(a^b).bit_count()<=6:d.join(i,j)
  components=defaultdict(list)
  for i,m in enumerate(ss):components[d.find(i)].append(m)
  cc[key]=list(components.values());counts[len(cc[key])]+=1;excess+=len(cc[key])-1
  for j,c in enumerate(cc[key]):
   for m in c:ix[m]=(key,j)
 banks[k]=g;comps[k]=cc;where[k]=ix;groups[k]={key:DSU(len(c))for key,c in cc.items()}
 need(sum(map(len,g.values()))==data['counts'][str(k)]['states'],'independent full boundary universe')
 need(sum(v for a,v in counts.items()if a>1)==data['counts'][str(k)]['bad_fibres']and excess==data['counts'][str(k)]['extra_components'],'independent small-move component excess')
 report.append({'margin':k,'states':sum(map(len,g.values())),'component_counts':dict(sorted(counts.items())),'component_excess':excess})
moves=data['moves'];need(len(moves)==220 and len({tuple(sorted(z))for z in moves})==220,'220 distinct unoriented move types')
reductions=0;costs=Counter()
for a,b in moves:
 need(a.bit_count()==b.bit_count()==4 and not a&b and q(a)==q(b),'disjoint four-per-side balanced native labels')
 need(work(a^b)<=22,'every bridge actual clock bound')
 key,ia=where[4][a];key2,ib=where[4][b];need(key==key2,'same complete margin4 fibre')
 need(groups[4][key].join(ia,ib),'every bridge removes one margin4 quotient component')
 reductions+=1;costs[work(a^b)]+=1
 for i in range(28):
  bit=1<<i
  if(a|b)&bit:continue
  c,d=a|bit,b|bit;key,ia=where[5][c];key2,ib=where[5][d];need(key==key2,'correct occupied background lift')
  groups[5][key].join(ia,ib)
for k in (4,5):
 for key,c in comps[k].items():need(len({groups[k][key].find(i)for i in range(len(c))})==1,'selected alphabet connects every boundary fibre')
need(reductions==report[0]['component_excess']==220,'sharp additional-type lower bound attained')
need(f(123,11)==6,'sharp witness coordinate reconciliation')
result={'status':'PASS','checks':checks,'scope':'Independent audit of concurrent sixslot570 sparse220 bridge completion, plus exact Q-coordinate isomorphism; full middle-grade support6 census not rerun here.','coordinate_map':'incoming packed p = ((x mod2) xor old_y_bit3) | ((old_y &7)<<1); cyclic x unchanged','boundary_census':report,'bridge_clock_histogram':dict(sorted(costs.items())),'sharp_extra_alphabet_size':220,'input_sha256':{str(z.relative_to(args.repo)):hashlib.sha256(z.read_bytes()).hexdigest()for z in [incoming/'bridges.json',incoming/'labels.txt',base/'native_packet_decomposition_20261023/native/domains.json']}}
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if args.output:args.output.write_text(s)
else:print(s,end='')
