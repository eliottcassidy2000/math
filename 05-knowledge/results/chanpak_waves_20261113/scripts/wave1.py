import itertools, subprocess
# (a) complete-Omega spectrum: H-values of tournaments with no two vertex-disjoint 3-cycles (orders <= 9)
def parse(bits,n):
    a=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[k]=='1': a[i]|=1<<j
            else: a[j]|=1<<i
            k+=1
    return a
comp=set(); allH=set()
for n in range(1,10):
    if n<3: comp.add(1); allH.add(1); continue
    for line in open(f'hl{n}.txt'):
        b,H=line.split()[:2]; H=int(H); allH.add(H)
        if H in comp: continue
        a=parse(b,n)
        tri=[(1<<x)|(1<<y)|(1<<z) for x,y,z in itertools.combinations(range(n),3)
             if (a[x]>>y&1 and a[y]>>z&1 and a[z]>>x&1) or (a[x]>>z&1 and a[z]>>y&1 and a[y]>>x&1)]
        if not any(s&t==0 for s,t in itertools.combinations(tri,2)): comp.add(H)
LIM=400
miss_all=[N for N in range(1,LIM,2) if N not in allH]
miss_comp=[N for N in range(1,LIM,2) if N not in comp]
print("odd N<400 not attained at order<=9:",miss_all[:30])
print("odd N<400 not attained with pairwise-meeting odd cycles (order<=9):",miss_comp[:40])
print("attained but only with disjoint triangles (order<=9), below 400:",[N for N in miss_comp if N in allH][:40])
