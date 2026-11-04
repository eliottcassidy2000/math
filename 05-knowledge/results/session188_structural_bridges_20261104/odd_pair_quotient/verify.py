"""Odd paired quotient trace, full coarse identities and named quartic control."""
from pathlib import Path
from tempfile import TemporaryDirectory
from functools import lru_cache
from itertools import combinations,product
import json,hashlib,subprocess
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def op(R):return tuple(((1<<len(R))-1)^r^(1<<i) for i,r in enumerate(R))
def inflate(Q):
 R=[0]*(2*len(Q))
 for i,row in enumerate(Q):
  R[2*i]|=1<<(2*i+1)
  for j in range(len(Q)):
   if row>>j&1:R[2*i]|=3<<(2*j);R[2*i+1]|=3<<(2*j)
 return tuple(R)
def ext(R,a,c):
 n=len(R);U=list(R)+[1<<(n+1),0]
 for u in range(n):
  for mask,v in ((a,n),(c,n+1)):
   if mask>>u&1:U[u]|=1<<v
   else:U[v]|=1<<u
 return tuple(U)
def hp(R):
 N=(1<<len(R))-1
 @lru_cache(None)
 def f(left,u):
  if not left:return 1
  available=left&R[u];total=0
  while available:
   bit=available&-available;available-=bit;total+=f(left^bit,bit.bit_length()-1)
  return total
 return sum(f(N^(1<<u),u) for u in range(len(R)))
def mobius(values):
 a=list(values);n=(len(a)-1).bit_length()
 for i in range(n):
  for mask in range(len(a)):
   if mask>>i&1:a[mask]-=a[mask^(1<<i)]
 return {m:v for m,v in enumerate(a) if v}
def allpaths(Q):
 paths=[]
 def walk(p,used):
  if len(p)>=2:paths.append((tuple(p),used))
  for v in range(len(Q)):
   if not used>>v&1 and Q[p[-1]]>>v&1:walk(p+[v],used|1<<v)
 for v in range(len(Q)):walk([v],1<<v)
 return paths
def pairpoly(Q):
 n=len(Q);poly={};paths=allpaths(Q);pairs=0
 for P,A in paths:
  for R,B in paths:
   if A&B:continue
   s,t=P[0],P[-1];u,v=R[0],R[-1];pairs+=1
   masks=[(1<<t)|(1<<(n+v)),(1<<s)|(1<<(n+u)),
          (1<<s)|(1<<t)|(1<<(n+v)),(1<<s)|(1<<t)|(1<<(n+u)),
          (1<<(n+u))|(1<<(n+v))|(1<<t),(1<<(n+u))|(1<<(n+v))|(1<<s)]
   for mask in masks:poly[mask]=poly.get(mask,0)^1
 return {m:v for m,v in poly.items() if v},pairs

def main():
 names=['PROOF.md','verify.py','query_bank.cpp','trace_mod8.cpp','trace_mod8_samples.cpp'];pins={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names}
 with TemporaryDirectory() as tmp:
  tmp=Path(tmp);binaries={}
  for source in names[2:]:
   out=tmp/source.removesuffix('.cpp');subprocess.run(['g++','-O3','-std=c++17',str(HERE/source),'-o',str(out)],check=True);binaries[source]=str(out)
  exhaustive=subprocess.check_output([binaries['trace_mod8.cpp']],text=True)
  samples=subprocess.check_output([binaries['trace_mod8_samples.cpp']],text=True)
  need('failures=0' in exhaustive and 'FAIL ' not in exhaustive,'exhaustive trace')
  need('failures=0' in samples and 'FAIL ' not in samples,'sample trace')
  def bank(R,queries):
   text=f'{len(R)} {len(queries)}\n'+ ' '.join(map(str,R))+'\n'+'\n'.join(f'{a} {c}' for a,c in queries)+'\n'
   return json.loads(subprocess.check_output([binaries['query_bank.cpp']],input=text,text=True))
  cases=[]
  for Q in [(0,),(2,0),(2,4,1),(20,1,2,7,14),(18,8,3,5,14),(24,17,3,6,12)]:
   q=len(Q);N=1<<q;T,U=inflate(Q),inflate(op(Q));coarse=[(m&(N-1),m>>q) for m in range(N*N)]
   expanded=[(sum(3<<(2*i) for i in range(q) if a>>i&1),sum(3<<(2*i) for i in range(q) if c>>i&1)) for a,c in coarse]
   ft,fu=bank(T,expanded),bank(U,expanded);delta=mobius([a-b for a,b in zip(ft,fu)])
   core=mobius([(hp(ext(Q,a,c))-hp(ext(op(Q),a,c)))//2 for a,c in coarse]);B,pairs=pairpoly(Q)
   for mask in range(N*N):
    if mask&(N-1) and mask>>q:need(delta.get(mask,0)%8==4*((core.get(mask,0)+B.get(mask,0))%2),'coefficientwise mixed mod8 identity')
   D=lambda i,j:delta.get((1<<i)|(1<<(q+j)),0)
   h=lambda i,j:delta.get((1<<i)|(1<<j)|(1<<(q+j)),0)
   L=[sum(D(i,j)+(h(i,j) if Q[j]>>i&1 else 0) for i in range(q) if i!=j) for j in range(q)]
   need(all(x%4==0 for x in L),'individual trace divisible4')
   if q>=3:need(sum(L)%8==(4*q)%8,'trace residue')
   need(max((m.bit_count() for m in delta),default=0)<=3,'coarse converse degree3')
   if q%2 and q>=3:need(any(a!=b for a,b in zip(ft,fu)),'odd quotient separator')
   # Independent backward-DP checks of actual expanded physical attachments.
   for idx in sorted(set(min(j,len(coarse)-1) for j in [0,len(coarse)-1,N+1,2*N+1,3*N+2])):
    a,c=expanded[idx];need(hp(ext(T,a,c))==ft[idx] and hp(ext(U,a,c))==fu[idx],'independent actual physical query')
   cases.append({'Q':Q,'T':T,'U':U,'physical_pair_blocks':[[2*i,2*i+1] for i in range(q)],'coarse_order':'mask=a+(c<<q)','actual_source_values':ft,'actual_target_values':fu,'coarse_nonzero_contrast_coefficients':[{'mask':m,'coefficient':v} for m,v in sorted(delta.items())],'all_unlisted_coefficients_zero':True,'L':L,'trace':sum(L),'disjoint_old_path_pairs':pairs})
  # Retain complete named cubic/quartic control for six old names, not just coarse masks.
  T,U=inflate((2,4,1)),inflate((4,1,2));n=len(T);N=1<<n;queries=[(m&(N-1),m>>n) for m in range(N*N)];ft,fu=bank(T,queries),bank(U,queries);coeff=mobius([a-b for a,b in zip(ft,fu)])
  need(max(m.bit_count() for m in coeff)<=4 and any(m.bit_count()==4 for m in coeff),'full named quartic retained')
  for mask in [0,1,5,1|(5<<n),N*N-1]:
   a,c=queries[mask];need(hp(ext(T,a,c))==ft[mask] and hp(ext(U,a,c))==fu[mask],'full named backward control')
  full={'T':T,'U':U,'mask_order':'a+(c<<6)','source_values':ft,'target_values':fu,'contrast_coefficients':[{'mask':m,'coefficient':v} for m,v in sorted(coeff.items())],'all_unlisted_coefficients_zero':True,'degree_counts':{str(d):sum(m.bit_count()==d for m in coeff) for d in range(5)}}
  (HERE/'coarse_cases.json').write_text(json.dumps(cases,separators=(',',':'))+'\n');(HERE/'full_named6.json').write_text(json.dumps(full,separators=(',',':'))+'\n')
 need(pins=={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names},'source stayed frozen')
 result={'status':'PASS','exhaustive_trace_stdout':exhaustive,'exhaustive_quotients':1098,'sample_trace_stdout':samples,'sample_quotients':16,'complete_coarse_cases':len(cases),'complete_coarse_queries':sum(len(c['actual_source_values']) for c in cases),'full_named_old_order':6,'full_named_queries':len(ft),'full_named_degree_counts':full['degree_counts'],'source_sha256':pins,'data_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ('coarse_cases.json','full_named6.json')}}
 (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
