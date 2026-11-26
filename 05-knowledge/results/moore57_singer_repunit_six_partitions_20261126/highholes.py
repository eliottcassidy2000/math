from itertools import combinations
from pathlib import Path
from collections import Counter
from singer_seed import C, N, read_tables
N54=(1<<N)-1
def holes(t):
 h={ij:[N54^m for m in b] for ij,b in t.items()}
 h.update({(j,i):[sum((m>>y&1)<<x for x,m in enumerate(h[i,j])) for y in range(N)]
           for i,j in C})
 return h
def bits(m):
 while m:
  z=m&-m;m^=z;yield z.bit_length()-1
def high(t):
 H=holes(t);qualified={}
 for I in combinations(range(6),4):
  extra=tuple(i for i in range(6) if i not in I)
  r=[-1]*6
  def go(i):
   if i==4:
    a,b=extra
    for x in range(N):
     k=sum(H[a,j][x]>>r[j]&1 for j in I)
     if k<2:continue
     r[a]=x
     for y in range(N):
      s=6+k+sum(H[b,j][y]>>r[j]&1 for j in I)+(H[a,b][x]>>y&1)
      if s>=13:
       r[b]=y;qualified[tuple(r)]=s
    return
   j=I[i];m=N54
   for k in I[:i]:m&=H[k,j][r[k]]
   for x in bits(m):
    r[j]=x;go(i+1)
  go(0)
 return Counter(qualified.values())
if __name__=='__main__':
 for i in range(4):
  t=read_tables(Path(__file__).with_name(f'b{i}.txt'));h=high(t);E=sum(N-a.bit_count() for b in t.values() for a in b)
  # Every six-vertex graph with 13 edges has a K4 (Turan, or exact 2^15 check).
  assert max(len(f) for mask in range(1<<15)
      if not any(all(mask>>C.index(e)&1 for e in combinations(I,2)) for I in combinations(range(6),4))
      for f in [tuple(e for k,e in enumerate(C) if mask>>k&1)])==12
  surplus=sum((k-12)*v for k,v in h.items())
  print(i, sorted(h.items()), 'remaining pairs', E,
        'minimum deleted OA rows',(E-surplus+11)//12,
        'maximum retained', N*N-(E-surplus+11)//12)
