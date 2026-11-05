from compiler import *
from itertools import combinations,permutations
from collections import Counter
import hashlib,json


def graphs(n):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  rows=[0]*n
  for k,(u,v) in enumerate(pairs):
   if bits>>k&1:u,v=v,u
   rows[u]|=1<<v
  yield tuple(rows)
def starts(r):
 n=len(r);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1 and S!=(1<<v):D[S][v]=sum(D[S^(1<<v)][u] for u in range(n) if S>>u&1 and r[v]>>u&1)
 return D[-1]
def cycles(r):
 if len(r)<3:return []
 return [(0,)+p for p in permutations(range(1,len(r))) if all(r[u]>>v&1 for u,v in zip((0,)+p,p+(0,)))]
def main():
 one=two=regular=graphs_checked=0;digest=hashlib.sha256();maxbirth=0
 for n in range(1,7):
  for r in graphs(n):
   graphs_checked+=1
   if len(scc(r,range(n)))!=1:continue
   if n==1:
    compile_one(r,{0});continue
   cs=cycles(r);st=starts(r)
   if len(cs)==1:
    roots=[u for u,x in enumerate(st) if x==1]
    for size in (1,2):
     for marked in combinations(roots,size):
      result=compile_one(r,marked);one+=1
      need(set(marked)<=set(result['ports']),'requested original vertices not retained')
      need(len(result['steps'])==n,'nonlinear birth count')
      digest.update(json.dumps(result,sort_keys=True).encode())
   for root,x in enumerate(st):
    if x==2:
     result=compile_two(r,root);two+=1;need(root in result['ports'],'root not retained')
     need(len(result['steps'])==n,'two-root birth count');digest.update(json.dumps(result,sort_keys=True).encode())
   if len(cs)==2 and all(x.bit_count()!=1 for x in r):
    result=compile_regular_five(r);regular+=1;digest.update(json.dumps(result,sort_keys=True).encode())
 # Audited width-two H67 word, independently reconstructed as FIFO two-port births.
 word=[0,0,1,1,1,3,2];r=[0]*len(word);ports=[];protocol=[]
 for v,mask in enumerate(word):
  for u in range(v):
   back=u in ports and bool(mask>>ports.index(u)&1)
   r[v if back else u]|=1<<(u if back else v)
  ports.append(v);drop=0 if len(ports)==3 else None
  if drop is not None:del ports[drop]
  protocol.append(dict(vertex=v,mask=mask,drop=drop,ports=ports[:]))
 r=tuple(r);need(r==(122,116,105,82,68,24,32),'67 graph')
 paths=[p for p in permutations(range(7)) if all(r[u]>>v&1 for u,v in zip(p,p[1:]))]
 need(len(paths)==67 and len(cycles(r))==4,'67 literal count and cycle scope')
 return dict(status='PASS',complete_labelled_tournaments=graphs_checked,
             start_one_requested_set_compilations=one,start_two_root_compilations=two,
             regular_five_labelled_compilations=regular,physical_protocol_sha256=digest.hexdigest(),
             regular_five_protocol=compile_regular_five(),H67=dict(rows=list(r),H=67,C=4,protocol=protocol))

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
