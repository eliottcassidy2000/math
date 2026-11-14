# Busch staircase: b(n) = min H over strong n-tournaments (3,5,9,15,25,45,75,125,...).
# Test: for each n, which odd N in [b(n), b(n+1)) are counts of strong tournaments of order exactly n / at most n.
import collections
def strong(bits,n):
    adj=[[False]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[k]=='1': adj[i][j]=True
            else: adj[j][i]=True
            k+=1
    def reach(s,a):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for v in range(n):
                if a[u][v] and v not in seen: seen.add(v); st.append(v)
        return seen
    if len(reach(0,adj))<n: return False
    rev=[[adj[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0,rev))==n
strongH=collections.defaultdict(set)
for n in range(3,10):
    for line in open(f"hl{n}.txt"):
        p=line.split()
        if strong(p[0],n): strongH[n].add(int(p[1]))
b={n:min(strongH[n]) for n in strongH}
print("Busch floor b(n):",[b[n] for n in sorted(b)])
cum=set()
for n in range(3,10):
    cum|=strongH[n]
    lo=b[n]; hi=b[n+1] if n+1 in b else 3*b[n]//1
    if n+1 not in b:  # n=9: next floor b(10)=125? use known Busch value
        hi=125
    tgt=list(range(lo|1,hi,2))
    miss_exact=[N for N in tgt if N not in strongH[n]]
    miss_cum=[N for N in tgt if N not in cum]
    print(f"n={n} step [{lo},{hi}): {len(tgt)} odd targets; missed at order exactly n: {miss_exact[:20]}; missed by orders<=n: {miss_cum}")

print()
for n in range(3,10):
    sp=sorted(strongH[n])
    # longest run of consecutive odd values
    best=(0,0,0); start=sp[0]; prev=sp[0]
    for v in sp[1:]+[None]:
        if v is not None and v==prev+2: prev=v; continue
        L=(prev-start)//2+1
        if L>best[0]: best=(L,start,prev)
        if v is not None: start=v; prev=v
    gaps_low=[N for N in range(sp[0],best[1],2) if N not in strongH[n]]
    print(f"n={n}: |spectrum|={len(sp)} range [{sp[0]},{sp[-1]}]; longest odd run [{best[1]},{best[2]}] ({best[0]} values); gaps below run: {gaps_low[:12]}{'...' if len(gaps_low)>12 else ''} ({len(gaps_low)})")
