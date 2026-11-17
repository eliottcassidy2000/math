from itertools import product, permutations

def matrices(c, o, s, r):
    st = [s%3, s//3%3, s//9]
    f, g = (c>>o)&1, (c>>(1-o))&1
    if not st[r] and (r == 1 and not f or r == 2 and not g):
        return None
    return s + 3**r * ((2-f,2-g,2)[r]-st[r])

def count(w, guard=True):
    h = 0
    for o in range(2):
        a = {0: 1}
        for c in w:
            b = {}
            for s, v in a.items():
                for r in range(3):
                    u = matrices(c, o, s, r)
                    if u is not None:
                        b[u] = b.get(u,0)+v
            a = b
        h += sum(v for s,v in a.items() if s%3!=1 and s//3%3!=1 and
                 (s//3%3!=0 or o==0) and (guard or all((s%3,s//3%3,s//9))))
    return h

def brute(w):
    n=len(w)
    def edge(a,b):
        if a<n and b<n: return a<b
        if a>=n and b>=n: return a<b
        if a>=n: return w[b]>>(a-n)&1
        return not (w[a]>>(b-n)&1)
    return sum(all(edge(a,b) for a,b in zip(p,p[1:])) for p in permutations(range(n+2)))

if __name__=='__main__':
    for m in range(1,30):
        for a,s in ((3,-1),(2,1)):
            w=(a,)+(1,)*m+(3,0)
            h=4*3**m + 9*2**m + s
            assert count(w)==h
            if m<5: assert brute(w)==h
    for c in (0,2,3):
        for n in range(1,7):
            breaks=[w for w in product((1,c),repeat=n) if count((2,)+w+(3,0))-count((3,)+w+(3,0))!=2]
            if breaks: print(c,breaks[0]);break
        else: print(c,'arbitrary candidate')

def rows(w):
    n=len(w)
    a=[sum(1<<j for j in range(i+1,n))+sum(1<<(n+k) for k in range(2) if not w[i]>>k&1)
       for i in range(n)]
    return a+[sum(1<<i for i in range(n) if w[i]&1)+(1<<(n+1)),
              sum(1<<i for i in range(n) if w[i]&2)]

def subset_count(a):
    n=len(a)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,h in enumerate(d[s]):
            if not h: continue
            nxt=a[i]&~s
            while nxt:
                b=nxt&-nxt; nxt-=b
                d[s|b][b.bit_length()-1]+=h
    return sum(d[-1])

def append(a, r=None):
    n=len(a)
    return [x + (1<<n if i!=r else 0) for i,x in enumerate(a)]+[0 if r is None else 1<<r]

def put(a, u, v):
    a[u] |= 1 << v

def pair(q, mark, mask, plus=False):
    """Replace mark by b,c,y; then add ports x,a. Keep exterior in index order."""
    U = [i for i in range(len(q)) if i != mark]
    b,c,y,x,a = range(len(U),len(U)+5)
    out = [0]*(len(U)+5)
    for i,u in enumerate(U):
        for j,v in enumerate(U):
            if q[u]>>v&1:
                put(out,i,j)
        for z in (b,c,y):
            put(out,i,z) if q[u]>>mark&1 else put(out,z,i)
        for z in (x,a):
            put(out,z,i) if mask>>i&1 else put(out,i,z)
    for u,v in [(b,c),(c,y),(y,b),(x,b),(a,b),(c,x),(a,c),(x,y),(y,a),
                 (a,x) if plus else (x,a)]:
        put(out,u,v)
    return out

def cover(a):
    from itertools import permutations
    n=len(a); z=[0]*n
    for w in permutations(range(n)):
        for mask in range(1<<(n-1)):
            if all((mask>>i&1) or (a[u]>>v&1) for i,(u,v) in enumerate(zip(w,w[1:]))):
                z[mask.bit_count()]+=1
    return z

def changed_slot(w, block, i=1):
    a=rows(w)
    pos=[j for j in range(len(a)) if j!=i]
    M=len(pos)
    out=[0]*(M+len(block))
    for j,u in enumerate(pos):
        for k,v in enumerate(pos):
            if a[u]>>v&1: put(out,j,k)
        for k in range(len(block)):
            put(out,j,M+k) if a[u]>>i&1 else put(out,M+k,j)
    for k,r in enumerate(block): out[M+k]|=r<<M
    return out

def direct_tests():
    from random import Random
    rng=Random(1925)
    cases=0
    for n in range(1,7):
        for _ in range(50):
            q=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    put(q,i,j) if rng.randrange(2) else put(q,j,i)
            d=rng.randrange(n); mask=rng.randrange(1<<(n-1))
            assert subset_count(pair(q,d,mask,True))-subset_count(pair(q,d,mask))==2*subset_count(q)
            z=cover(q)+[0]*3
            for w,k in [((2,1,0,1),27),((3,1,3,0),29),((2,1,3,0),31)]:
                assert subset_count(changed_slot(w,q))==k*z[0]+21*z[1]+4*z[2]
            cases+=1
    for p in range(1,5):
        q=[sum(1<<j for j in range(i+1,p)) for i in range(p)]
        for d in range(p):
            for mask in range(1<<(p-1)):
                assert subset_count(pair(q,d,mask,True))-subset_count(pair(q,d,mask))==2
    print(cases,'arbitrary quotients and block substitutions passed')

if __name__=='__main__':
    direct_tests()
