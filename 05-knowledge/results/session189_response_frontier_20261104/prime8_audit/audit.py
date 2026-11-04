from pathlib import Path
from functools import lru_cache
import hashlib,itertools,json,math,subprocess,sys,tempfile
D=Path(__file__).resolve().parent
P=D.parent/'prime8'
if not P.exists():P=D.parent/'session189-prime8'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def valid(T):
 n=len(T);return all(0<=r<(1<<n) and not(r>>i&1) for i,r in enumerate(T)) and all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))
def prime(T):
 n=len(T)
 for s in range(1,(1<<n)-1):
  if s.bit_count()<2:continue
  if all((T[i]&s) in (0,s) for i in range(n) if not s>>i&1):return False
 return True
def converse(T):return tuple(((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T))
# Independently recurse from a fixed first symbol by remaining multiplicities.
def starts(T,repeat=-1):
 n=len(T);counts=tuple(1+(i==repeat) for i in range(n))
 @lru_cache(None)
 def complete(left,last):
  if not any(left):return 1
  return sum(complete(left[:j]+(left[j]-1,)+left[j+1:],j) for j in range(n) if left[j] and T[last]>>j&1)
 return [complete(counts[:i]+(counts[i]-1,)+counts[i+1:],i) for i in range(n)]
def extend(T,c):
 n=len(T);return tuple(r|(((c>>i)&1)<<n) for i,r in enumerate(T))+(((1<<n)-1)^c,)

def main():
 cover=json.loads((P/'cover.json').read_text());levels=cover['levels'];need([x['order'] for x in levels]==list(range(9)),'all rooted levels0..8');need(levels[0]['representatives']==[[]],'empty base')
 counts=[len(l['representatives']) for l in levels];need(counts==[1,1,1,2,4,12,56,456,6880],'specified level sizes')
 for l in levels:
  need(all(valid(T) and len(T)==l['order'] for T in l['representatives']),'valid representative tournaments')
 transitions=0
 for n in range(1,9):
  prev=levels[n-1]['representatives'];current=levels[n]['representatives'];records=levels[n]['extension_records'];need(len(records)==len(prev)*(1<<(n-1)),'complete rooted extension domain')
  for k,(target,packed) in enumerate(records):
   parent,mask=divmod(k,1<<(n-1));T=prev[parent];G=[r|((1-((mask>>j)&1))<<(n-1)) for j,r in enumerate(T)]+[mask];p=[(packed>>(3*j))&7 for j in range(n)]
   need(sorted(p)==list(range(n)) and packed>>(3*n)==0,'valid physical permutation');need(0<=target<len(current),'target index')
   U=current[target];need(all(((G[p[i]]>>p[j])&1)==((U[i]>>j)&1) for i in range(n) for j in range(n)),'explicit rooted extension isomorphism');transitions+=1
 all8=json.loads((P/'all_representatives.json').read_text());need(all8==levels[8]['representatives'],'output representative identity')
 Q=json.loads((P/'quotients.json').read_text());primes=[T for T in all8 if prime(T)];need([q['rows'] for q in Q]==primes and len(Q)==4008,'independent exact prime filter')
 with tempfile.TemporaryDirectory(prefix='prime8-independent-audit-') as td:
  exe=Path(td)/'check';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'independent.cpp'),'-o',str(exe)],check=True)
  def run(rows,mode=None):
   data=str(len(rows))+'\n'+''.join(str(len(r))+' '+' '.join(map(str,r))+'\n' for r in rows);cmd=[str(exe)]+([mode] if mode else []);return [list(map(int,l.split())) for l in subprocess.check_output(cmd,input=data.encode()).decode().splitlines()]
  allrows=[T for l in levels for T in l['representatives']];canon=run(allrows);at=0;orbit=[];canon_records=[]
  for n,l in enumerate(levels):
   rows=canon[at:at+len(l['representatives'])];at+=len(rows);need(len(set(r[0] for r in rows))==len(rows),'independent canonical integer codes are disjoint')
   need(all(a>0 and math.factorial(n)%a==0 for c,a in rows),'automorphism orders divide factorial');total=sum(math.factorial(n)//a for c,a in rows);need(total==1<<(n*(n-1)//2),'complete labeled orbit mass');orbit.append(total);canon_records.append(dict(order=n,columns=['canonical_edge_code','automorphism_order'],records=rows))
  J=run(primes,'joint');need(len(J)==4008,'all prime joint endpoint matrices');pzero=[]
  for q,j in zip(Q,J):
   need(len(j)==64,'complete joint endpoint array');S=[sum(j[8*i:8*i+8]) for i in range(8)];E=[sum(j[i::8]) for i in range(8)];p=[(a+b)%2 for a,b in zip(S,E)];need(S==q['S'] and E==q['E'] and p==q['p'],'independent complete endpoints')
   if not any(p):pzero.append(q)
  need(len(pzero)==2,'exactly two endpoint-zero prime types')
 # The full retained finite anchors, independently counted by recursive completions.
 A=json.loads((P/'anchor_exploration.json').read_text());need([a['rows'] for a in A]==[q['rows'] for q in pzero],'precise residual representative set');anchor_records=[];primitive_checks=0;patterns=0
 for q,a in zip(pzero,A):
  T=tuple(q['rows']);U=converse(T);sk=[[0]*8 for _ in range(8)]
  for j in range(8):
   s=starts(T,j);u=starts(U,j)
   for i in range(8):sk[i][j]=s[i]-u[i]
  b=[];d=[]
  for i in range(8):
   R=extend(T,1<<i);S=extend(U,1<<i);b.append(sum(starts(R))-sum(starts(S)));d.append([sum(starts(R,j))-sum(starts(S,j)) for j in range(8)])
  need(sk==a['start_run_skews'] and b==a['one_port_base'] and d==a['one_port_double'],'all independently recounted anchor primitives');primitive_checks+=136
  expected=[]
  for V in range(256):
   first=[(b[i]+2*sum(d[i][j] for j in range(8) if V>>j&1))%4 for i in range(8)];t=[(q['S'][i]-q['E'][i]+2*sum(sk[i][j] for j in range(8) if V>>j&1))%4 for i in range(8)];second=[(2*t[i]+2*first[i])%8 for i in range(8)];need(any(first) or any(second),'every physical singleton-position pattern is separated');expected.append(dict(V=V,first=first,second=second));patterns+=1
  need(expected==a['patterns'] and a['blind']==[],'all512 complete anchor records');anchor_records.append(dict(rows=T,start_run_skews=sk,one_port_base=b,one_port_double=d,patterns=expected))
 # These complete independently generated arrays supplement, rather than replace, producer data.
 outputs={'canonical_orbits.json':canon_records,'joint_endpoints.json':dict(rows=primes,joint_endpoint_matrices=J),'anchor_recount.json':anchor_records}
 for name,obj in outputs.items():
  data=json.dumps(obj,separators=(',',':'))+'\n';path=D/name
  if path.exists():need(path.read_text()==data,'audit artifact is stable: '+name)
  else:path.write_text(data)
 return {'status':'PASS','level_sizes':counts,'rooted_extension_isomorphisms_checked':transitions,'order8_extensions_checked':58368,'labeled_orbit_totals':orbit,'prime8_representatives':4008,'joint_endpoint_coordinates':4008*64,'marginal_endpoint_coordinates':4008*16,'endpoint_zero_representatives':[q['rows'] for q in pzero],'independent_anchor_primitive_coordinates':primitive_checks,'singleton_position_patterns':patterns,'source_sha256':hashlib.sha256((D/'audit.py').read_bytes()).hexdigest(),'independent_cpp_sha256':hashlib.sha256((D/'independent.cpp').read_bytes()).hexdigest(),'input_pins':{name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['cover.json','all_representatives.json','quotients.json','anchor_exploration.json']},'output_pins':{name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
