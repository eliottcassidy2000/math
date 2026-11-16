from functools import lru_cache

def cuts(rows):
    n=len(rows); full=(1<<n)-1
    h=[0]*(full+1); h[0]=1
    all_d=[]
    for b in range(n):
        d=[[0]*n for _ in h]; d[1<<b][b]=1
        for s in range(1,full+1):
            if not s>>b&1: continue
            for a,k in enumerate(d[s]):
                if not k: continue
                h[s]+=k
                t=rows[a] & (full^s)
                while t:
                    bit=t&-t; t-=bit
                    d[s|bit][bit.bit_length()-1]+=k
        all_d.append(d)
    w=[[sum(d[s][a]*h[full^s] for s in range(1,full+1) if not s.bit_count()%2) for d in all_d] for a in range(n)]
    k=[sum(w[a][b] for a in range(n) if mask>>a&1 for b in range(n) if not mask>>b&1) for mask in range(full+1)]
    return h[full],w,k

def extend(rows,x):
    n=len(rows)
    return tuple(row|((x>>a&1)<<n) for a,row in enumerate(rows))+(((1<<n)-1)^x,)

def direct(rows):
    n=len(rows); d=[[0]*n for _ in range(1<<n)]
    for a in range(n): d[1<<a][a]=1
    for s in range(1,1<<n):
        for a,k in enumerate(d[s]):
            if k:
                for b in range(n):
                    if rows[a]>>b&1 and not s>>b&1: d[s|1<<b][b]+=k
    return sum(d[-1])

cases=[
 ('locked49',(2,4,73,19,39,15,59),7),
 ('released49',(2,12,17,37,11,87,31),7),
 ('rigid63',(102,252,248,241,225,192,128,1),3),
 ('rigid63',(102,252,248,241,225,192,128,1),7),
 ('delay343',(24,25,67,84,100,15,163,63),7),
]
if __name__=='__main__':
    for name,rows,p in cases:
        h,w,k=cuts(rows)
        n=len(rows)
        assert all(direct(extend(rows,x))==h+2*v for x,v in enumerate(k))
        lo=min([(v,x) for x,v in enumerate(k) if (h+2*v)%p],default=None)
        if lo:
            v,x=lo
            print(name,p,h,'min',h+2*v,'mask',x,'attaining',sum(t==v for t in k),'best ≤2',min(h+2*k[y] for y in range(1<<n) if min(y.bit_count(),n-y.bit_count())<=2 and (h+2*k[y])%p))
        else: print(name,p,h,'none')

def below_triple(rows):
    h=direct(rows); cap=3*h
    front=[tuple(rows)]; layers=[]; seen=set(front)
    while front:
        layer={}; nxt=[]
        for t in front:
            d,w,k=cuts(t)
            layer[d]=layer.get(d,0)+1
            for x in range(1,len(k)-1):
                f=d+2*k[x]
                if f>=cap: continue
                assert f>d
                u=extend(t,x)
                if u not in seen:
                    seen.add(u); nxt.append(u)
        layers.append(layer)
        front=nxt
    return layers

if __name__=='__main__':
    for name,rows,p in cases:
        if name=='rigid63' and p==7: continue
        got=below_triple(rows)
        print(name,'depth layer sizes',[sum(r.values()) for r in got], 'totals',sorted(set().union(*got)),flush=True)

def wrap(rows):
    n=len(rows)
    # T -> a -> b -> T
    return tuple(r|(1<<n) for r in rows)+(1<<(n+1),(1<<n)-1)

def colours(rows):
    n=len(rows); full=(1<<n)-1
    def induced(s):
        v=[j for j in range(n) if s>>j&1]
        return tuple(sum((rows[a]>>b&1)<<i for i,b in enumerate(v)) for a in v)
    hs=[direct(induced(s)) if s else 1 for s in range(full+1)]
    @lru_cache(None)
    def q(t,s):
        if not t: return int(s==0)
        total=0; x=s
        while True:
            total+=hs[x]*q(t-1,s^x)
            if x==0: return total
            x=(x-1)&s
    return lambda t:q(t,full)

if __name__=='__main__':
    for name,rows,_ in cases+ [('trap75',(2,12,9,113,39,71,23),7)]:
        q=colours(rows); t=wrap(rows); qq=colours(t)
        assert [qq(j) for j in range(1,5)]==[j*j*q(j)+j*(q(j+1)-q(j-1)) for j in range(1,5)]
        print(name,'wrap counts',[q(1),qq(1),direct(wrap(t))], 'initial mod7',[q(j)%7 for j in range(1,4)])
    for M in (8,9,18,35,97,100,200):
        W=105
        if __import__('math').gcd(M,W)!=1: continue
        a=-2*pow(M,-1,W)
        for r in range(W):
            x=a*r%W
            assert (__import__('math').gcd(r*(M-r),W)==1)==(__import__('math').gcd(x*(x+2),W)==1)

def transitive(n):
    return tuple(((1<<n)-1)^((1<<(j+1))-1) for j in range(n))

if __name__=='__main__':
    t=transitive(5)
    h,w,k=cuts(t)
    assert h==1 and k[16]==8 and k[20]==9
    assert [w[4][b] for b in range(4)]==[4,2,1,1]
    assert sum(w[4][b] for b in (0,1,3))==7
    assert sum(w[2][b] for b in (0,1,3))==2
    for n in range(2,10):
        _,_,k=cuts(transitive(n))
        z=1<<(n-1)
        assert [i for i in range(n-1) if k[z|(1<<i)]==k[z]+1]==([2] if n in (4,5) else [])
    assert direct(extend(transitive(7),98))==73
