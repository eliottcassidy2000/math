"""Independent native-complement seam identity and exact infinite-padding controls."""
from itertools import permutations
import json

def need(x,m):
    if not x: raise RuntimeError(m)
def hp(t):
    n=len(t); d=[[0]*n for _ in range(1<<n)]
    for v in range(n): d[1<<v][v]=1
    for s in range(1,1<<n):
        for v,c in enumerate(d[s]):
            if c:
                a=t[v]&~s
                while a:
                    b=a&-a;a-=b; w=b.bit_length()-1
                    d[s|b][w]+=c
    return sum(d[-1])
def strong(t):
    n=len(t)
    for v in range(n):
        seen=1<<v
        while True:
            new=seen
            for w in range(n):
                if seen>>w&1:new|=t[w]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True
def partner(t,p):
    n=len(t); r=[((1<<n)-1)^ (1<<v)^t[v] for v in range(n)]
    for u,v in zip(p,p[1:]):
        need(t[u]>>v&1,'Hamilton path guard')
        r[u]|=1<<v;r[v]&=~(1<<u)
    return tuple(r)
def pad(t,p):
    n=len(t); r=p[0]
    return tuple(row|((1<<n) if v!=r else 0) for v,row in enumerate(t))+(1<<r,), (n,)+p

def main():
    t=tuple(sum(1<<((v+d)%5) for d in (1,2)) for v in range(5))
    paths=[p for p in permutations(range(5)) if p[0]==0 and all(t[u]>>v&1 for u,v in zip(p,p[1:]))]
    counts={p:hp(partner(t,p)) for p in paths}
    p=paths[0];q=next(q for q in paths if counts[q]!=counts[p])
    original=(t,p,q); rows=[]
    for k in range(5):
        a,b=partner(t,p),partner(t,q)
        need(strong(t),'source strong')
        need((hp(a),hp(b))==(counts[original[1]],counts[original[2]]),'padding preserves exact H')
        changed={(i,j) for i in range(len(t)) for j in range(i+1,len(t)) if ((a[i]^b[i])>>j)&1}
        seams={tuple(sorted(e)) for e in zip(p,p[1:])}^{tuple(sorted(e)) for e in zip(q,q[1:])}
        need(changed==seams,'exact seam symmetric difference')
        rows.append(dict(n=len(t),H=[hp(a),hp(b)],changed_edges=sorted(changed)))
        t2,p=pad(t,p); _,q=pad(t,q);t=t2
    print(json.dumps(dict(source_rows=original[0],paths=original[1:],initial_partners=[partner(original[0],v) for v in original[1:]],padding=rows),indent=2))
if __name__=='__main__':main()
