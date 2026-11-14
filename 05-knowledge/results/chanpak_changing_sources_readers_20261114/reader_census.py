"""Compute exact longest consecutive even-forest coordinate in marked reader fibers."""
from collections import defaultdict, Counter
from pathlib import Path
from sys import argv

D = defaultdict(dict)
sizes = Counter()
for line in Path(argv[1]).read_text().splitlines():
    n,q,h,a,b,*R = row = tuple(map(int,line.split()))
    assert n==len(R)
    f=h+(a+b)//2
    assert f>0
    D[h,f//(f&-f)][b//2]=row
    sizes[n]+=1
windows = []
for (h,f),profiles in D.items():
    for e in profiles:
        if e-1 in profiles:
            continue
        k=1
        while e+k in profiles:
            k+=1
        windows.append((k,h,f,e))
print('Number of reader rows by order:',sorted(sizes.items()))
for k,h,f,e in sorted(windows,reverse=True)[:3]:
    print('W',k,h,f,e)
    for i in range(k):
        print(D[h,f][e+i])
