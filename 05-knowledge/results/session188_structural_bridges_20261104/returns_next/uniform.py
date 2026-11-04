"""Exact compressed-anchor and alias-cubic checks; standard library only."""
import itertools,json,hashlib
from functools import lru_cache
from pathlib import Path
from verify import rows_from_mask,modules,reverse,need
from probe import hp,ext
@lru_cache(None)
def value(R,a,c):return hp(ext(R,a,c))
@lru_cache(None)
def delta(R,a,c):
 O=reverse(R,[(1<<len(R))-1]);d=value(R,a,c)-value(O,a,c);need(d%2==0,'integral half difference');return (d//2)%2

def endpoint_b(R):
 n=len(R);b=[0]*n
 for p in itertools.permutations(range(n)):
  if all(R[u]>>v&1 for u,v in zip(p,p[1:])):b[p[0]]^=1;b[p[-1]]^=1
 return b

def packings(Ms):
 out=[]
 def go(i,covered,chosen):
  if i==len(Ms):out.append(tuple(chosen));return
  go(i+1,covered,chosen)
  if not covered&Ms[i]:go(i+1,covered|Ms[i],chosen+[Ms[i]])
 go(0,0,[]);return out

def cubic(R,u,w,v):
 z=0
 for i,j,k in itertools.product((0,1),repeat=3):z^=delta(R,((1<<u) if i else 0)|((1<<w) if j else 0),(1<<v) if k else 0)
 return z

def run():
 counts={'labelled_sources':0,'module_anchor_cubic_sums':0,'compressed_anchor_batteries':0,'actual_uniform_queries':0,'compressed_q2_batteries':0};examples={}
 for n in range(3,6):
  N=(1<<n)-1
  for code in range(1<<(n*(n-1)//2)):
   R=rows_from_mask(n,code);O=reverse(R,[N]);Ms=[M for M in modules(R) if M!=N];b=endpoint_b(R);counts['labelled_sources']+=1
   for M in Ms:
    vertices=[u for u in range(n) if M>>u&1]
    for v in range(n):
     if M>>v&1:continue
     need(sum(cubic(R,u,w,v) for u,w in itertools.combinations(vertices,2))%2==0,'alias cubic cancellation')
     counts['module_anchor_cubic_sums']+=1
   for pack in packings(Ms):
    covered=0
    for M in pack:covered|=M
    free=N^covered
    if not free:continue
    blocks=list(pack)+[1<<v for v in range(n) if free>>v&1];q=len(blocks);restored=reverse(O,pack)
    for v in range(n):
     if not free>>v&1:continue
     qs=[];d=0
     for C in blocks:
      if C==1<<v:continue
      rep=(C&-C).bit_length()-1
      if R[v]>>rep&1:qs.append(C|(1<<v));d+=1
      else:qs.append(C)
     if d%2:qs.append(1<<v)
     if (q-1+d)%2:qs.append(0)
     need(len(qs)==len(set(qs)),'distinct battery masks')
     need(sum(delta(R,a,1<<v) for a in qs)%2==(1+b[v])%2,'compressed anchor parity')
     need(delta(R,1<<v,0)==b[v],'endpoint certificate')
     need(any(delta(R,a,1<<v) for a in qs) or delta(R,1<<v,0),'actual uniform separator')
     for a,c in [(a,1<<v) for a in qs]+[(1<<v,0)]:
      need(value(O,a,c)==value(restored,a,c),'each uniform query invariant under module restoration')
      counts['actual_uniform_queries']+=1
     counts['compressed_anchor_batteries']+=1
     if q==2:counts['compressed_q2_batteries']+=1
     if 'q2' not in examples and q==2:examples['q2']={'rows':R,'blocks':blocks,'anchor':v,'d':d,'a_queries':qs,'b_anchor':b[v]}
 result={'status':'PASS','counts':counts,'examples':examples,'source_sha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest() for p in ('UNIFORM.md','uniform.py','verify.py','probe.py')}}
 Path(__file__).with_name('uniform.receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
