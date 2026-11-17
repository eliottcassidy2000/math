from pathlib import Path
from math import isqrt
from replay import table
def ck(v,*x):
    if not v:raise ArithmeticError(x)
for line in Path(__file__).with_name('deep_targets.out').read_text().splitlines():
    prefix,full=line.split(' literal ')
    n,N,k,q=map(int,prefix.split()[:4])
    rows=tuple(map(int,full.split(',')))
    ck(n==len(rows) and rows[-1]==q,n,N,rows)
    ck(all(x>>i & 1==0 for i,x in enumerate(rows)))
    ck(all(((rows[i]>>j)^(rows[j]>>i))&1 for i in range(n) for j in range(i)),rows)
    ck(table(rows)[-1]==N,rows,N)
    print(N,'literal order',n,'birth interval',(12,n))
for p in (100003,30011,30013):
    ck(all(p%d for d in range(2,isqrt(p)+1)),p)
