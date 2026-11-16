from functools import lru_cache
from itertools import permutations,product,combinations
from collections import Counter
import csv,json,sys
from pathlib import Path
from program import dp,prof

def subset(t,ps):
 return tuple(sum(((t[u]>>v)&1)<<i for i,v in enumerate(ps)) for u in ps)
@lru_cache(None)
def canon(t):
 n=len(t)
 cells=[tuple(i for i,a in enumerate(t) if a.bit_count()==d) for d in range(n)]
 return min(subset(t,sum(v,())) for v in product(*(tuple(permutations(c)) for c in cells if c)))
def op(t):
 n=len(t)
 return tuple(((1<<n)-1)^v^(1<<i) for i,v in enumerate(t))
def packs(t):
 cycles=[];nums=[]
 for m in (3,5,7):
  k=0
  for S in combinations(range(len(t)),m):
   for u in permutations(S[1:]):
    p=(S[0],)+u
    if all(t[i]>>j&1 for i,j in zip(p,p[1:]+p[:1])):
     cycles.append(sum(1<<i for i in S));k+=1
  nums.append(k)
 return tuple(nums),sum(not (a&b) for j,a in enumerate(cycles) for b in cycles[:j])

def main(root):
 data={7:{canon(tuple(map(int,s.split()))) for s in (root/'n7_representatives.txt').read_text().splitlines()}}
 for n in range(6,0,-1):
  data[n]={canon(subset(t,tuple(i for i in range(n+1) if i!=j))) for t in data[n+1] for j in range(n+1)}
 rows=lambda z:tuple((z>>(8*j))&255 for j in range(8))
 data[8]={rows(int(c['code'])) for c in csv.DictReader((root/'n8_classes.csv').open())}
 assert [len(data[n]) for n in range(1,9)]==[1,1,2,4,12,56,456,6880]
 birth={};allcounts={}
 for n in range(1,9):
  allcounts[n]={t:dp(t) for t in sorted(data[n])}
  for D in allcounts[n].values():birth.setdefault(D[-1],n)
 out=['h\torder\toutgoing_masks']
 for n in range(1,9):
  counts=set();hist=Counter();orphan=rescued=0
  for t,D in allcounts[n].items():
   h=D[-1];N=(1<<n)-1
   if birth[h]!=n:continue
   out.append(f'{h}\t{n}\t'+','.join(map(str,t)))
   counts.add(h)
   p=sum(birth.get(D[N^(1<<i)])==n-1 for i in range(n)) if n>1 else 0
   hist[p]+=1
   if not p and n>2:
    orphan+=1
    if any(birth[D[N^(1<<i)^(1<<j)]]==n-2 for i in range(n) for j in range(i)):
     rescued+=1
  print(n,len(counts),sum(hist.values()),dict(sorted(hist.items())),orphan,rescued)
 print('values',len(birth),'classes',len(out)-1)
 assert len(birth)==320 and len(out)==5828
 assert sorted(h for h in range(1,610,2) if h not in birth)==[7,21]
 out[1:]=sorted(out[1:],key=lambda line:(int(line.split()[1]),int(line.split()[0]),line))
 Path(__file__).with_name('births.tsv').write_text('\n'.join(out)+'\n')
 for h in (63,125,127,255,511):
  n=birth[h]
  F=[t for t,d in allcounts[n].items() if d[-1]==h]
  st=set()
  for t in F:
   for x in range(n):
    _,g,r,s=prof(t,x)
    if r and s-r==g:st.add(2*r)
  print(h,n,len(F),'positive slopes at birth', sorted(st))
 for h,w,expected in ((125,(6,12,56,97,11,83,23),((12,28,10),6)),(127,(10,20,33,70,77,27,39),((11,23,11),9))):
  assert allcounts[7][canon(w)][-1]==h and packs(w)==expected
 print('Odd cycle packs at 125/127:',packs((6,12,56,97,11,83,23)),packs((10,20,33,70,77,27,39)))
 fan=tuple((((1<<9)-1)^((1<<(i+1))-1)) for i in range(9))
 q=list(fan)
 for t in (3,5,8):q[0]^=1<<t;q[t]^=1
 assert prof(q,0)==(127,1,107,108)
 for h in (125,127):
  F=[t for t,d in allcounts[8].items() if d[-1]==h]
  stable={2*r for t in F for x in range(8) for _,g,r,s in (prof(t,x),) if r and s-r==g}
  print(h,len(F),'order-8 slopes',sorted(stable))
if __name__=='__main__':main(Path(sys.argv[1]))
