# Shallow-cut horizon X(m): least odd N > 67 that is not h + 2^j L (j>=1) for a rooted tournament of order <= m.
import sys
LIM=500000
pairs=set((1,2**t) for t in range(40))
res={}
for n in range(3,10):
    for line in open(f"hl{n}.txt"):
        parts=line.split(); H=int(parts[1])
        for item in parts[2:]:
            v,rest=item.split(':'); h,L,E,S=map(int,rest.split(','))
            pairs.add((H,L))
    cov=bytearray(LIM//2+2)
    for (H,L) in pairs:
        j=1
        while H+(L<<j)<=LIM:
            cov[(H+(L<<j))>>1]=1; j+=1
    miss=[N for N in range(69,LIM+1,2) if not cov[N>>1]]
    first=miss[0] if miss else None
    # density of misses in [first, 2*first]
    print(f"order<={n}: pairs={len(pairs):8d}  horizon X={first}  misses<=LIM: {len(miss)}  first few: {miss[:8]}")
