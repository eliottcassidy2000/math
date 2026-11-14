import itertools
def parse(bits,n):
    a=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[k]=='1': a[i]|=1<<j
            else: a[j]|=1<<i
            k+=1
    return a
tau1=set(); tau2=set(); allH=set()
for n in range(1,10):
    if n<3: tau1.add(1); tau2.add(1); allH.add(1); continue
    for line in open(f'hl{n}.txt'):
        b,H=line.split()[:2]; H=int(H); allH.add(H)
        if H in tau2 and H in tau1: continue
        a=parse(b,n)
        ts=[(1<<x)|(1<<y)|(1<<z) for x,y,z in itertools.combinations(range(n),3) if (a[x]>>y&1 and a[y]>>z&1 and a[z]>>x&1) or (a[x]>>z&1 and a[z]>>y&1 and a[y]>>x&1)]
        if any(all(t>>v&1 for t in ts) for v in range(n)) or not ts: tau1.add(H); tau2.add(H); continue
        if H not in tau2 and any(all(t&((1<<u)|(1<<v)) for t in ts) for u,v in itertools.combinations(range(n),2)): tau2.add(H)
LIM=125
print("odd N<125 not attained (order<=9):",[N for N in range(1,LIM,2) if N not in allH])
print("odd N<125 not attained with a 1-vertex triangle transversal:",[N for N in range(1,LIM,2) if N not in tau1])
print("odd N<125 not attained with a <=2-vertex triangle transversal:",[N for N in range(1,LIM,2) if N not in tau2])
big=[N for N in range(125,700,2) if N not in tau2]
print("odd N in [125,700) not attained with tau<=2 at order<=9 (incomplete range):",len(big),big[:20])
