from functools import lru_cache

WITNESSES=((4588920935503036672,212),(4579966650245382400,88),
           (4582218210346205440,210),(4579966650245382400,210))
def from_code(code, outmask):
    return tuple(((code>>(8*i))&255)|(((outmask>>i&1)^1)<<8) for i in range(8))+(outmask,)
def wrapper(t):
    n=len(t)
    return tuple(a|1<<n for a in t)+(1<<(n+1),(1<<n)-1)
def path_array(t):
    n=len(t)
    assert all(not(t[i]>>i&1) for i in range(n))
    assert all(((t[i]>>j)^(t[j]>>i))&1 for i in range(n) for j in range(i))
    @lru_cache(None)
    def anchored(s,i):
        s^=1<<i
        return sum(anchored(s,j) for j in range(n) if s&t[i]>>0 & (1<<j)) if s else 1
    return (1,)+tuple(sum(anchored(s,i) for i in range(n) if s>>i&1) for s in range(1,1<<n))
def c2(h):
    m=len(h)-1
    return sum(h[s]*h[m^s] for s in range(1,m))
if __name__=='__main__':
    rows=[]
    for code,m in WITNESSES:
        a=from_code(code,m)
        h=path_array(a)
        w=path_array(wrapper(a))
        w2=path_array(wrapper(wrapper(a)))
        assert w[-1]==3*h[-1]+c2(h)
        assert w2[-1]==3*w[-1]+c2(w)
        rows.append((h[-1],c2(h),w[-1],c2(w),w2[-1]))
    print(*rows,sep='\n')
    assert [a[0] for a in rows]==[367,369,371,373]
    assert [a[1] for a in rows]==[10122]*4
    assert [a[2] for a in rows]==[11223,11229,11235,11241]
    assert len(set(a[3] for a in rows))>1
    def slot_poly(t,j):
        h=path_array(t)
        @lru_cache(None)
        def q(r,s):
            if not r:return int(not s)
            v=s;ans=0
            while True:
                ans+=h[v]*q(r-1,s^v)
                if not v:return ans
                v=(v-1)&s
        return q(j,len(h)-1)
    for i,(code,m) in enumerate(WITNESSES):
        a=from_code(code,m)
        for j in range(10):
            assert slot_poly(a,j)==j**9+30*j**7+162*j**5+162*j**3+(12+2*i)*j
    t=[]
    for i in range(6):
        t.append(sum(1<<j for j in range(i+1,6)))
    t[0]^=1<<5;t[5]^=1
    assert path_array(tuple(t))[-1]==17
    assert slot_poly(tuple(t),2)==224
    assert path_array(wrapper(tuple(t)))[-1]==241
    assert pow(2,4,17)==16 and pow(2,12,241)==240 and pow(2,8,241)==15
    print('All slot polynomials, triangle wrapper transfer, and 17→241 replayed.')
    SMALL=((62,60,88,112,96,68,3),(62,28,88,112,32,70,19),
           (62,28,72,80,36,78,19),(62,76,72,80,38,78,17),(62,76,104,16,38,74,25))
    assert all(sum((x^y).bit_count() for x,y in zip(a,b))==4 for a,b in zip(SMALL,SMALL[1:]))
    for i,t in enumerate(SMALL):
        for j in range(8):
            assert slot_poly(t,j)==j**7+20*j**5+60*j**3+(10+2*i)*j
        assert c2(path_array(t))==1086
        a=path_array(wrapper(t))
        assert a[-1]==1359+6*i
        assert c2(a)==19566+12*i
        assert path_array(wrapper(wrapper(t)))[-1]==3*(1359+6*i)+19566+12*i
    print('Order-seven five-state palette replayed; each move flips two interior arcs.')
