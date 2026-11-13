"""Exact endpoint-transport profiles for the 46 inherited core-changing rules."""
from pathlib import Path
from collections import defaultdict
from itertools import product
from math import comb
import gzip,json,sys,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OLD=ROOT/'05-knowledge/results/birth_grammar_repairs_20261112'
sys.path.insert(0,str(OLD))
from sparse_covers import kernel

def need(b,msg):
 if not b:raise RuntimeError(msg)
def read(name):return json.loads(gzip.decompress((OLD/name).read_bytes()))
def exponential(K):
 E=defaultdict(int)
 for r,c in K.items():
  for bs in product(*(range(1,x+1) for x in r)):
   z=c
   for x,b in zip(r,bs):z*=(-1)**(x-b)*comb(x,b)
   E[bs]+=z
 return {b:c for b,c in E.items() if c}
def diagonal(E):
 D=defaultdict(int)
 for (a,b),c in E.items():D[a*b]+=c
 return {b:c for b,c in D.items() if c}
def difference(A,B):
 return {b:A.get(b,0)-B.get(b,0) for b in A.keys()|B.keys() if A.get(b,0)!=B.get(b,0)}
def stirling_positive(E):
 # b^m=sum_r binom(b,r) r!S(m,r); r0 vanishes for m>=1.
 R=max(E)
 out=[sum(c*comb(b,r) for b,c in E.items() if b>=r) for r in range(1,R+1)]
 return out

def run():
 rows=read('mixed_grammar.json.gz')['records'];out=[];even=[];begin=time.time()
 for i,r in enumerate(read('synchrony_kernels.json.gz')['records']):
  ports=r['ports'][:2];O=[v for v in range(len(r['rows'])) if r['selected']>>v&1]
  d,u=kernel(r['rows'],ports,O);de=diagonal(exponential(d));ue=diagonal(exponential(u))
  need(de=={1:1},'even whole diagonal defect')
  growth={b:c*(b-1) for b,c in ue.items() if b!=1};cert=stirling_positive(growth)
  need(all(c>=0 for c in cert) and cert[0]>0,'even strict growth')
  lower=difference(ue,{2:1});lc=stirling_positive(lower)
  need(all(c>=0 for c in lc),'even U(m)>=2^m')
  even.append({'family':i,'rows':r['rows'],'selected':r['selected'],'ports':ports,'U':sorted(ue.items()),'growth_Stirling':cert,'lower_bound_Stirling':lc})
 for i,r in enumerate(rows):
  rr=r['rows'];O=[v for v in range(len(rr)) if r['selected']>>v&1];ports=r['coupled_ports']
  need(len(rr)==13 and not r['selected']>>12,'unselected last vertex')
  parent=[x&((1<<12)-1) for x in rr[:12]]
  base=even[r['parent_synchrony_family']]
  labels=[v for v in range(10) for _ in range(2 if v in base['ports'] else 1)]
  literal=[sum(1<<j for j,w in enumerate(labels) if (v==w and i<j) or (v!=w and base['rows'][v]>>w&1)) for i,v in enumerate(labels)]
  need(parent==literal,'actual parent graph equals even length2 expansion')
  need({labels[p] for p in ports}==set(base['ports']),'one changed representative in each inherited coupled block')
  need(r['selected']==sum(1<<j for j,v in enumerate(labels) if base['selected']>>v&1),'actual selected labels retained')
  d,u=kernel(rr,ports,O);pd,pu=kernel(parent,ports,O)
  de=diagonal(exponential(d));pe=diagonal(exponential(pd))
  need(de==pe=={1:1},'both complete diagonal seed identities')
  ue=diagonal(exponential(u));pue=diagonal(exponential(pu));diff=difference(ue,pue);cert=stirling_positive(diff)
  ep=dict(even[r['parent_synchrony_family']]['U'])
  need(pue=={b:b*c for b,c in ep.items()},'exact parent parameter m+1 transport')
  growth={b:c*(b-1) for b,c in ue.items() if b!=1};gc=stirling_positive(growth)
  need(all(c>=0 for c in gc) and gc[0]>0,'odd strict growth')
  lc=stirling_positive(difference(ue,{2:1}));need(all(c>=0 for c in lc),'odd U(m)>=2^m')
  need(all(x>=0 for x in cert) and cert[0]>0,'strict all-positive-length U descent')
  # Parent same coupled blocks; decrement synchronized length to its m=1 base.
  pdiff={b:c*(b-1) for b,c in pue.items() if b!=1};pcert=stirling_positive(pdiff)
  need(all(x>=0 for x in pcert) and pcert[0]>0,'strict diagonal clock decrement')
  out.append({'family':i,'parent_synchrony_family':r['parent_synchrony_family'],'rows':rr,'selected':r['selected'],'ports':ports,'child_U':sorted(ue.items()),'parent_U':sorted(pue.items()),'decrement':sorted(diff.items()),'positive_Stirling':cert,'parent_growth_Stirling':pcert,'growth_Stirling':gc,'lower_bound_Stirling':lc})
  print(json.dumps({'family':i,'decrement':sorted(diff.items()),'seconds':round(time.time()-begin,2)}),flush=True)
 data={'scope':'All46 inherited mixed-family presentations; exact all-positive synchronized length m. Deleting the retained last unselected vertex preserves delta1, reduces U by positive certified profile, and retains actual marked endpoint source data. Not arbitrary target coverage.','records':out,'even_records':even}
 name='profiles'+('.optimized' if not __debug__ else '')+'.json.gz'
 (HERE/name).write_bytes(gzip.compress((json.dumps(data,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
if __name__=='__main__':run()
