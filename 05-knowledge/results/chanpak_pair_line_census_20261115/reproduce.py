#!/usr/bin/env python3
from itertools import combinations
from collections import defaultdict
import random
def stats(n,code):
  out=[0]*n
  for j,(i,k) in enumerate(combinations(range(n),2)):
    if code>>j&1: i,k=k,i
    out[i]|=1<<k
  path=[[0]*n for _ in range(1<<n)]
  H=[0]*(1<<n);H[0]=1
  cycle=[0]*(n+1)
  for m in range(1,1<<n):
    if m&(m-1)==0: path[m][m.bit_length()-1]=1
    else:
      for z in range(n):
        if m>>z&1:
          pre=path[m^(1<<z)]
          path[m][z]=sum(pre[i] for i in range(n) if out[i]>>z&1)
    H[m]=sum(path[m])
    if m.bit_count()>=3:
      start=(m&-m).bit_length()-1
      prev=[[0]*n for _ in range(1<<n)]
      prev[1<<start][start]=1
      # Traverse only subsets containing the smallest vertex, in ascending order.
      for sub in range(1,m+1):
        if sub&m==sub and sub>>start&1 and sub!=1<<start:
          for z in range(n):
            if sub>>z&1 and z!=start:
              pre=prev[sub^(1<<z)]
              prev[sub][z]=sum(pre[i] for i in range(n) if out[i]>>z&1)
      cycle[m.bit_count()]+=sum(prev[m][i] for i in range(n) if out[i]>>start&1)
  full=(1<<n)-1
  Q=sum(H[s]*H[full^s] for s in range(1<<n))
  return H[-1],Q,cycle
def main():
  tested=0
  for n in range(2,8):
    k=n*(n-1)//2
    rng=random.Random(115203+n)
    codes=range(1<<k) if n<=5 else [rng.randrange(1<<k) for _ in range(200)]
    for c in codes:
      h,q,cyc=stats(n,c);tested+=1
      assert q%3 == ((-1)**n*h)%3
      assert q%4==0
      if n%2==0 and n>=4:
        assert q%8==0 and q//8%2==cyc[n-1]%2
      if n%2:
        assert q//4%2==cyc[n]%2
  print('profile tests',tested)
  for n in (3,5,6,7,8,9,10):
    # Reverse just (0,n-1) at odd n and (0,n-1),(1,n-1) at even n.
    code=(1<<(n-2)) + ((1<<(2*n-4)) if n%2==0 else 0)
    h,q,cyc=stats(n,code)
    x=n if n%2 else n-1
    assert q%(8 if n%2 else 16)==(4 if n%2 else 8)
    assert cyc[x]==(1 if n%2 else n-1)
    print('sharp witness',n,code,h,q,cyc[x])
  data=defaultdict(dict)
  for c in range(1<<10):
    h,q,_=stats(5,c)
    data[h].setdefault(q%8,(c,q))
  print('same-order H with distinct Q mod 8:',{h:d for h,d in data.items() if len(d)>1})
if __name__=='__main__': main()
