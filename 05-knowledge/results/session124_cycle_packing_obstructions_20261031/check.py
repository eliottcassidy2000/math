"""Independent complete presentations of tournaments with triangle-packing <= 1.
No canonical-class inventory is assumed: choose an actual triangle, order its
transitive complement, and enumerate all 8^m cross-edge words.
"""
from itertools import combinations, permutations
from collections import Counter
from pathlib import Path
import json
P=Path(__file__).resolve().parent
checks=0
def need(c,s):
 global checks
 checks+=1
 if not c: raise RuntimeError(s)
def graph(m,code):
 n=m+3; out=[0]*n
 out[0]|=2;out[1]|=4;out[2]|=1
 for i in range(m):
  v=i+3;letter=code>>(3*i)&7
  for a in range(3):
   if letter>>a&1:out[a]|=1<<v
   else:out[v]|=1<<a
  for j in range(i+1,m):out[v]|=1<<(j+3)
 return out
def triangles(out):
 return [sum(1<<x for x in vs) for vs in combinations(range(len(out)),3)
         if all(out[vs[i]]>>vs[(i+1)%3]&1 for i in range(3)) or all(out[vs[(i+1)%3]]>>vs[i]&1 for i in range(3))]
def hp(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,c in enumerate(dp[mask]):
   if not c:continue
   nxt=out[v]&~mask
   while nxt:
    bit=nxt&-nxt;nxt-=bit;w=bit.bit_length()-1
    dp[mask|bit][w]+=c
 return sum(dp[-1])
def oddcycles(out):
 # Each literal directed cycle starts at its minimum vertex, removing rotations.
 n=len(out);cnt=Counter();supports=[]
 def extend(path,mask):
  v=path[-1];root=path[0]
  if len(path)>=3 and len(path)%2 and out[v]>>root&1:
   cnt[len(path)]+=1;supports.append(mask)
  nxt=out[v]&~mask&~((1<<(root+1))-1)
  while nxt:
   bit=nxt&-nxt;nxt-=bit;w=bit.bit_length()-1;extend(path+[w],mask|bit)
 for v in range(n):extend([v],1<<v)
 return cnt,supports
def main():
 rows=[];spectrum={1};witnesses={};total=0
 for m in range(6):
  hist=Counter();max_tau=0
  for code in range(8**m):
   total+=1;out=graph(m,code);ts=triangles(out)
   if any(not(a&b) for a,b in combinations(ts,2)):continue
   h=hp(out);cycles,supports=oddcycles(out)
   need(h==1+2*sum(cycles.values()),'literal OCF')
   need(all(a&b for a,b in combinations(supports,2)),'all odd cycle pairs intersect')
   need(all(t&7 for t in ts),'triangle core is feedback set')
   tau=next(k for k in range(1,4) if any(all(t&sum(1<<v for v in vs) for t in ts) for vs in combinations(range(m+3),k)))
   max_tau=max(tau,max_tau);hist[h]+=1;spectrum.add(h)
   witnesses.setdefault(h,dict(n=m+3,code=code,out=out,cycle_lengths=dict(cycles),feedback_number=tau))
  rows.append(dict(n=m+3,presentations=8**m,packing_free_presentations=sum(hist.values()),spectrum=sorted(hist),max_feedback_number=max_tau))
 need(total==37449,'presentation total')
 need(5**8>73**3,'Busch all-order cutoff')
 allowed=[h for h in range(1,75,2) if h not in (7,21)]
 misses=[h for h in allowed if h not in spectrum]
 need(misses==[29,35,41,43,55,57,59,67,71],'all-order relative obstructions')
 # Strong witness showing29 genuinely occurs, found independently on six vertices.
 w29=None
 for code in range(1<<15):
  out=[0]*6
  for k,(u,v) in enumerate(combinations(range(6),2)):
   if code>>k&1:out[u]|=1<<v
   else:out[v]|=1<<u
  if hp(out)==29:w29=out;break
 need(w29 is not None,'29 global witness')
 cycles,supports=oddcycles(w29);pairs=sum(not(a&b) for a,b in combinations(supports,2))
 need(1+2*sum(cycles.values())+4*pairs==29,'29 OCF')
 need(pairs>0,'29 needs disjoint cycles')
 return dict(status='PASS',checks=checks,total_presentations=total,rows=rows,
  all_order_packing_free_spectrum_through73=sorted(spectrum&set(range(1,75,2))),
  globally_attainable_but_packing_free_impossible=misses,
  witness29=dict(out=w29,cycles=dict(cycles),disjoint_odd_cycle_pairs=pairs),
  scope='All-order claims use the proved representation lemma, one nontrivial SCC, and cited Busch bound; finite computations independently enumerate actual graphs.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
