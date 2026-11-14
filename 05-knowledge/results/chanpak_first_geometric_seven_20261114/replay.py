from itertools import product, permutations
from functools import cache
from collections import Counter
from math import isqrt

def check(x):
    if not x: raise AssertionError("rotor fiber control")

def rotor(D):
    return tuple(sum(1<<((i+d)%7) for d in D) for i in range(7))

def table(T):
    n=len(T); ans=[1]*(1<<n)
    p=[[0]*n for _ in range(1<<n)]
    for S in range(1,1<<n):
        u=S
        while u:
            b=u&-u;u-=b;i=b.bit_length()-1;R=S^b
            if not R:p[S][i]=1
            else:
                v=R&T[i]
                while v:
                    c=v&-v;v-=c
                    p[S][i]+=p[R][c.bit_length()-1]
        ans[S]=sum(p[S])
    return ans

def quotient(r, triangle=False):
    n=3 if triangle else 4
    # 0=v,1=a,2=b,3=x; the old triangle is v->b->a->v.
    a=[4,1,2,0][:n]
    if not triangle:
        for i in range(3):
            if r>>i&1:a[3]|=1<<i
            else:a[i]|=8
    return tuple(a)

def inflate(D,Q):
    t=rotor(D)
    bins=[tuple(range(7))]+[(i,) for i in range(7,7+len(Q)-1)]
    v=list(t)+[0]*(len(Q)-1)
    for i,P in enumerate(bins):
        q=sum(1<<w for j,B in enumerate(bins) if Q[i]>>j&1 for w in B)
        for u in P:v[u]|=q
    return tuple(v)

def strong(T):
    n=len(T); mask=(1<<n)-1
    for sign in (0,1):
        reach=1
        for _ in T:
            for i,v in enumerate(T):
                if reach>>i&1:reach|=v if not sign else v^mask^(1<<i)
        if reach!=mask:return False
    return True

def paths(T):
    n=len(T); F=(1<<n)-1
    def rec(p,S):
        if S==F:
            yield p
        else:
            c=(F^S)&T[p[-1]]
            while c:
                b=c&-c;c^=b
                yield from rec(p+(b.bit_length()-1,),S|b)
    for a in range(n):
        yield from rec((a,),1<<a)

def flip(T,p):
    u=list(T)
    for i,j in zip(p,p[1:]):
        u[i]^=1<<j;u[j]^=1<<i
    return tuple(u)

def rotate_word(p,d):
    return tuple((i+d)%7 if i<7 else i for i in p)

orbit_groups={frozenset(frozenset(k*x%7 for x in D) for k in range(1,7)) for D in ((1,2,3),(1,2,4))}
check(len(orbit_groups)==2)
check(set().union(*orbit_groups)=={frozenset(j if b else 7-j for j,b in enumerate(z,1))
                                  for z in product((0,1),repeat=3)})
for i in range(8):
    g=quotient(i)
    nr=[sum(1 for w in set(permutations([0]*u+[1,2,3]))
            if all(g[x]>>y&1 for x,y in zip(w,w[1:]))) for u in (1,2)]
    check(nr==([3,1] if i in (0,7) else [5,2] if i in (4,5) else [5,3]))
for D,h,q,expect in [
    ((1,2,3),175,1988,(309,827,827,827,593,593,827,309)),
    ((1,2,4),189,2016,(315,837,837,837,603,603,837,315))]:
    l=table(rotor(D))
    check(l[-1]==h and sum(l[i]*l[127^i] for i in range(128))==q)
    check(table(inflate(D,quotient(0,True)))[-1]==h+q)
    for i,e in enumerate(expect):
        g=quotient(i)
        T=inflate(D,g)
        check(table(T)[-1]==7*e)
        check(strong(g)==(i not in [0,7]))
check(all(593%d for d in range(2,isqrt(593)+1)))
T=inflate((1,2,3),quotient(4))
p=list(paths(T))
cycs=[w for w in p if w[0]==9 and T[w[-1]]>>w[0]&1]
closed=[w for w in p if T[w[-1]]>>w[0]&1]
check((len(p),len(cycs),len(closed))==(4151,175,1750))
reps=[w for w in p if w==min(rotate_word(w,d) for d in range(7))]
check(len(reps)==593)
C=Counter();O=Counter()
N=[]
for w in reps:
    f=flip(T,w)
    if T[w[-1]]>>w[0]&1:
        C[strong(f)]+=1
        if not strong(f):N.append(table(f)[-1])
    else:
        O[strong(f)]+=1
check(C==Counter({True:175,False:75}) and O==Counter({True:343}))
print('Literal minimal-symmetry atlas, 4151 paths, 175 cycles, and 250/343 path-orbit split passed.')
print('Closed partner connectivity:',C,'; open partner connectivity:',O)
print('Nonstrong partners at closed lines (path count, frequency):',sorted(Counter(N).items()))
