from itertools import permutations,combinations
from collections import Counter
import json

def edge(a,i,j): return (a[i]>>j)&1
def paths(a): return [p for p in permutations(range(len(a))) if all(edge(a,i,j) for i,j in zip(p,p[1:]))]
def iso(a,b): return [p for p in permutations(range(len(a))) if all(edge(a,i,j)==edge(b,p[i],p[j]) for i in range(len(a)) for j in range(len(a)))]
def converse(a): return [sum(1<<j for j in range(len(a)) if i!=j and not edge(a,i,j)) for i in range(len(a))]
def cycles(a,k):
 return [p for p in permutations(range(len(a)),k) if p[0]==min(p) and all(edge(a,p[i],p[(i+1)%k]) for i in range(k))]
def modules(a):
 return [s for k in range(2,len(a)) for s in combinations(range(len(a)),k) if all(len({edge(a,x,v) for v in s})==1 for x in range(len(a)) if x not in s)]
def rows(n,m):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if m>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
A=[14,28,24,16,1]; B=[10,28,25,16,1]
out={}
for name,a in [('A',A),('B',B)]:
 out[name]={'rows':a,'paths':paths(a),'scores':sorted(x.bit_count() for x in a),'automorphisms':iso(a,a),'anti_automorphisms':iso(a,converse(a)),'cycles':{k:cycles(a,k) for k in range(3,6)},'modules':modules(a)}
out['A_to_B_isomorphisms']=iso(A,B)
out['small_labelled_H_counts']={n:dict(sorted(Counter(len(paths(rows(n,m))) for m in range(1<<(n*(n-1)//2))).items())) for n in range(1,6)}
def require(condition, message):
 if not condition: raise RuntimeError(message)
for name in ('A','B'):
 r=out[name]
 require(len(r['paths'])==9, name+' H=9')
 require(len(r['automorphisms'])==1, name+' Aut=1')
 require(len(r['anti_automorphisms'])==1, name+' self-converse')
 require(r['scores']==[1,1,2,3,3], name+' score sequence')
require(not out['A_to_B_isomorphisms'], 'nonisomorphism')
require(len(out['A']['cycles'][4])==3 and len(out['B']['cycles'][4])==2, '4-cycle separation')
require(out['A']['anti_automorphisms']==[(4,3,2,1,0)], 'A specified anti-map')
require(out['B']['anti_automorphisms']==[(0,4,3,2,1)], 'B specified anti-map')
require(set(out['small_labelled_H_counts'][4])=={1,3,5}, 'order-four spectrum')
require(all(9 not in out['small_labelled_H_counts'][n] for n in range(1,5)), 'H9 not earlier')
require(out['small_labelled_H_counts'][5][9]==240, 'exact H9 fiber exhaustion')
require(out['small_labelled_H_counts'][1]=={1:1}, 'H1 birth')
require(out['small_labelled_H_counts'][3]=={1:6,3:2}, 'H3 birth')
require(out['small_labelled_H_counts'][4][5]==24, 'H5 unique birth')
out['verification']='PASS'
print(json.dumps(out,indent=2))
