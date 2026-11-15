from itertools import combinations, permutations
from math import factorial


def tt(k):
    return [(1<<k)-(1<<(i+1)) for i in range(k)]


def profile(rows, max_j):
    if any(a>>i&1 for i,a in enumerate(rows)) or any(((rows[i]>>j)^(rows[j]>>i))&1!=1 for i,j in combinations(range(len(rows)),2)):
        raise ValueError(('not a tournament', rows))
    n=len(rows); F=(1<<n)-1
    reach=[]
    for i in range(n):
        prev=0; cur=1<<i
        while prev!=cur:
            prev=cur
            for j,r in enumerate(rows):
                if cur>>j&1:cur|=r
        reach.append(cur)
    h=[1]+[0]*F
    d=[[0]*n for _ in range(F+1)]
    for s in range(1,F+1):
        for i,r in enumerate(rows):
            if not s>>i&1:continue
            t=s^(1<<i)
            d[s][i]=sum(d[t][j] for j in range(n) if (t&r)>>j&1) if t else 1
        h[s]=sum(d[s])
    out=[h[F]]
    a=h
    for j in range(2,max_j+1):
        b=[0]*(F+1)
        for s in range(F+1):
            t=s
            while True:
                b[s]+=a[s^t]*h[t]
                if not t:break
                t=(t-1)&s
        out.append(b[F]);a=b
    tri=[s for s in range(F+1) if s.bit_count()==3 and h[s]==3]
    return tuple(out), all(x==F for x in reach), any(x&y==0 for x in tri for y in tri)


def ordinal(u,w):
    a=len(u)
    return [r|(((1<<len(w))-1)<<a) for r in u]+[r<<a for r in w]


def wrap(u,w):
    # U -> singleton -> W -> U
    a,b=len(u),len(w)
    return [r|(1<<a) for r in u]+[((1<<b)-1)<<(a+1)]+[(r<<(a+1))|((1<<a)-1) for r in w]


def ports(k):
    # q_i -> p_i, p_j -> q_i for all j != i
    a=tt(k); b=[r<<k for r in tt(k)]
    for i in range(k):
        a[i]|=1<<(k+i)
        b[i]|=((1<<k)-1)^(1<<i)
    return a+b


def replace(u,k):
    # Module -> q_i; p_i -> module
    n=len(u); p=ports(k)
    return [a|(((1<<k)-1)<<n) for a in u]+[(r<<n)|(((1<<n)-1) if i>=k else 0) for i,r in enumerate(p)]


def skeleton(k):
    rr=ports(k); c=[0]*(k+2)
    for p in permutations(range(2*k)):
        compulsory=[]; optional=[]
        for i in range(2*k+1):
            left=p[i-1] if i else -1
            right=p[i] if i<2*k else -1
            allowed=(left<0 or left>=k) and (right<0 or right<k)
            adjacency=(left<0 or right<0 or rr[left]>>right&1)
            if not allowed and not adjacency:break
            if not adjacency:compulsory.append(i)
            elif allowed:optional.append(i)
        else:
            for m in range(1<<len(optional)):
                j=len(compulsory)+m.bit_count()
                if j:c[j]+=1
    return c[1:]


def old_profiles(bound):
    # Whole old alphabet: singleton SCCs and one-apex strong SCCs.
    # Lower floor Q2 >= 2*3^(n-3) limits each factor.
    factors={(1,2)}
    for n in range(3,40):
        if 2*3**(n-3)>bound:break
        k=n-3
        for bits in range(1<<k):
            w='1'+(f'{bits:0{k}b}' if k else '')+'0'
            h,q=1,2**n
            for i,x in enumerate(w):
                for j,y in enumerate(w[i+1:],i+1):
                    if x!='1' or y!='0':continue
                    d=j-i-1
                    h+=max(2,2**d)
                    q+=2**(n-d-2)*(3**d+1)
            if q<=bound:factors.add((h,q))
    seen={(1,1)}; todo=[(1,1)]
    for h,q in todo:
        for a,b in factors:
            key=h*a,q*b
            if key[1]<=bound and key not in seen:
                seen.add(key);todo.append(key)
    return seen


def guard_rows(k,w):
    # Hub i bit at spine position j selects hub -> spine.
    n=2*k+len(w)
    s=tt(n); h=tt(k)
    for i in range(k):
        for j,mask in enumerate([(1<<k)-1]*k+list(w)+[0]*k):
            if mask>>i&1:h[i]|=1<<(k+j)
    for j,m in enumerate([(1<<k)-1]*k+list(w)+[0]*k):
        s[j]=(s[j]<<k)|(((1<<k)-1)^m)
    return h+s


def main():
    # Full literal profiles and physical matched-port readouts.
    u=[92,141,32,228,206,19,38,101]
    w=[52,229,208,231,10,148,49,81]
    for b,want in [(u,45225),(w,44505)]:
        got,strong,_=profile(b,3)
        if got!=(345,5640,want) or not strong:raise AssertionError((got,strong))
    old=old_profiles(5640)
    if (345,5640) in old:raise AssertionError('old profile')
    print('Old alphabet distinct profiles through Q2=5640 (incl empty):',len(old))
    print('Literal new-slot profile switch: (345,5640,45225) vs (345,5640,44505)')
    for k in range(1,5):
        c=skeleton(k)
        if c[-1]!=factorial(k) or any(z%factorial(i) for i,z in enumerate(c)):
            raise AssertionError((k,c))
        print('Portal skeleton',k,c)
    p1,p2=(profile(replace(b,2),1)[0][0] for b in (u,w))
    if p1-p2!=1440:raise AssertionError((p1,p2))
    print('Matched-portal literal HP counts:',p1,p2)
    rr=[24,53,57,50,32,1]; r2=[4,1,2,7]
    if profile(rr,2)!=((17,224),True,True) or profile(r2,2)[0]!=(3,24):
        raise AssertionError('borrow profiles')
    for j in range(5,8):
        a=profile(wrap(ordinal(rr,tt(j-5)),tt(2)),1)[0][0]
        b=profile(wrap(tt(j),r2),1)[0][0]
        if (a,b)!=(21*2**j-17,21*2**j-15):raise AssertionError((j,a,b))
    print('Exact borrow pair literal at j=5,6,7; stronger core has two disjoint triangles.')
    for k in range(1,4):
        pieces=[]
        m=3
        for i in range(1<<(k*m)):
            w=[i>>(k*j)&((1<<k)-1) for j in range(m)]
            got,strong,rich=profile(guard_rows(k,w),2)
            n=2*k+m
            if not strong or (k>=2 and not rich):raise AssertionError((k,w,got,strong,rich))
            for z,a in enumerate(got,1):
                if not (2*z)**k*(k+z)**max(n-3*k,0)<=a<=(2*z*(k+1))**k*(k+z)**n:
                    raise AssertionError((k,w,z,a))
            pieces.append(got)
        print('Guarded k=%d: %d literal configurations, %d distinct (H,Q2) profiles'% (k,len(pieces),len(set(pieces))))

if __name__=='__main__':main()
