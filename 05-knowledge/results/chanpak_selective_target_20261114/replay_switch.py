"""Literal tournament DP for two-tape selective cuts and their inverse."""
from functools import cache
from math import isqrt

def check(ok):
    if not ok:
        raise AssertionError('literal cut check failed')

def tt(n):
    return tuple((1 << n) - (2 << i) for i in range(n))

def join(a,b):
    return tuple(x | (((1 << len(b))-1) << len(a)) for x in a)+tuple(
        x << len(a) for x in b)

def apex(b,m):
    return tuple(x|((1 << len(b)) if not (m >> i)&1 else 0)
                 for i,x in enumerate(b))+(m,)

def endpoint_dp(rows):
    n=len(rows)
    full=(1<<n)-1
    check(all(not(x>>i&1) for i,x in enumerate(rows)))
    check(all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i in range(n) for j in range(i)))
    @cache
    def p(s,i,j):
        if s == 1 << i:
            return int(i==j)
        if i==j:
            return 0
        z=s^(1<<i)
        rem=z&rows[i]
        c=0
        while rem:
            u=(rem&-rem).bit_length()-1
            rem^=1<<u
            c+=p(z,u,j)
        return c
    @cache
    def h(s):
        if not s:
            return 1
        return sum(p(s,i,j) for i in range(n) if s>>i&1
                             for j in range(n) if s>>j&1)
    def marked(u):
        e0=o0=e1=o1=c=0
        for s in range(1,full+1):
            k=h(full^s)
            for i in range(n):
                if not s>>i&1:
                    continue
                for j in range(n):
                    if not s>>j&1:
                        continue
                    part=p(s,i,j)*k
                    if s.bit_count()%2:
                        o0+=part*((u>>i)&1)
                        o1+=part*(1-(u>>j&1))
                    else:
                        e0+=part*((u>>i)&1)
                        e1+=part*(1-(u>>j&1))
                        c+=part*((u>>i)&1)*(1-(u>>j&1))
        return h(full),sum(h(s)*h(full^s) for s in range(full+1)),e0,o0,e1,o1,c
    return h(full), marked

def inv(N,Q,b,d,offs):
    found=set()
    for z in offs:
        u=Q*(N-z)+b*d
        for f in range(1,isqrt(max(0,u))+1):
            if u%f:
                continue
            for x,y in [(f,u//f),(u//f,f)]:
                if (x+d)%Q or (y+b)%Q:
                    continue
                X,Y=(x+d)//Q,(y+b)//Q
                if X>=2 and Y>=2 and X&(X-1)==0 and Y&(Y-1)==0:
                    found.add((X.bit_length()-1,Y.bit_length()-1,z))
    return found

C5=tuple((1<<((i+3)%5))|(1<<((i+4)%5)) for i in range(5))
group6=[((40,37,17,6,11,28),23),((48,49,3,7,12,28),15),
        ((48,49,3,7,12,28),27)]
replays=0
for group,coeffs,offs in [
    ([(C5,17),(C5,9)],(120,57,78),[71,73]),
    (group6,(432,303,174),[225,227,229])]:
    for k,((core,mask),c0) in enumerate(zip(group,offs)):
        h,Q,e0,o0,e1,o1,c=endpoint_dp(core)[1](mask)
        check((Q,Q-e1-o1,Q-e0-o0)==coeffs)
        check(h+Q+2*(c-e0-e1)==c0)
        for l,m in [(1,1),(1,2),(2,1),(2,3),(3,2)]:
            old=join(join(tt(l),core),tt(m))
            actual=endpoint_dp(apex(old,((1<<l)-1)|(mask<<l)))[0]
            X,Y=1<<l,1<<m
            expected=Q*X*Y-coeffs[1]*X-coeffs[2]*Y+c0
            check(actual==expected)
            check((l,m,c0) in inv(actual,*coeffs,offs))
            replays+=1
        for l,m,lo,ro in [(1,1,0,1),(2,1,2,0),(1,2,1,1),(2,2,1,2)]:
            a0,m0=group[0]
            baseline=endpoint_dp(apex(join(join(tt(l),a0),tt(m)),lo|(m0<<l)|(ro<<(l+len(a0)))))[0]
            nxt=endpoint_dp(apex(join(join(tt(l),core),tt(m)),lo|(mask<<l)|(ro<<(l+len(core)))))[0]
            check(nxt-baseline==2*k)
            replays+=1
print(f'{replays} literal complete graph checks; factor inverses and arbitrary tape controls passed.')
