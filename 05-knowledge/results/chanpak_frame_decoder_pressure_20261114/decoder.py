from functools import cache
from math import isqrt
from itertools import product

FLAGS=((1,1),(47,56))
RESET=(0,13)

@cache
def kernel(q):
    if q<0: raise ValueError(q)
    e,f,g=(1,1,1) if q==0 else (1,2,2) if q==1 else (3,3,4) if q==2 else tuple(
        sum(t) for t in zip(kernel(q-1),kernel(q-2),kernel(q-3)))
    return e,f,g

def mul(m,n):
    a,b,c,d=m
    e,f,g,h=n
    return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)

def matrix(q):
    a,b,d=kernel(q)
    return a,b,b,d

def powexp(x):
    return x.bit_length()-1 if x>=2 and x&(x-1)==0 else None

def coeff(m,l,r):
    a,b,c,d=m
    lu,lv=FLAGS[l]
    ru,rv=FLAGS[r]
    return a*lu*ru,c*lv*ru,b*lu*rv,d*lv*rv

def value(c,x,y):
    a,b,f,d=c
    return a+b*x+f*y+d*x*y

def inv2(n,c):
    a,b,f,d=c
    if n<value(c,2,2): return
    # floor(log2(n/d)) without division by a floating-point number.
    t=(n//d).bit_length()-1
    for s in (t-1,t):
        if s<2: continue
        p=1<<s
        u=n-a-d*p
        if u<=0: continue
        delta=u*u-4*b*f*p
        if delta<0: continue
        j=isqrt(delta)
        if j*j!=delta: continue
        for z in set((u+j,u-j)):
            if z%(2*b): continue
            x=z//(2*b)
            a0=powexp(x)
            if a0 is not None and 1<=a0<s and value(c,x,1<<(s-a0))==n:
                yield a0,s-a0

def direct_word(gaps):
    assert gaps and gaps[0]>0 and gaps[-1]>0
    v=[1,1<<gaps[0]]
    q=0
    for g in gaps[1:]:
        if not g:
            q+=1
            continue
        a,b,d=kernel(q)
        v=[a*v[0]+b*v[1],(b*v[0]+d*v[1])<<g]
        q=0
    return sum(v)

@cache
def primitive(n,l,r):
    # Return best gap word, i.e. minimum cost then lex order under this boundary.
    best=None
    lu,lv=FLAGS[l]
    ru,rv=FLAGS[r]
    rem=n-lu*ru
    if rem>0 and rem%(lv*rv)==0:
        e=powexp(rem//(lv*rv))
        if e is not None: best=(e,)
    def visit(m,inner):
        nonlocal best
        for a,b in inv2(n,coeff(m,l,r)):
            best=better(best,(a,)+inner+(b,))
        # The former terminal positive site becomes an internal one of width e.
        e=1
        while True:
            pre=mul(m,(1,0,0,1<<e))
            if value(coeff(mul(pre,matrix(1)),l,r),2,2)>n: break
            q=1
            while True:
                new=mul(pre,matrix(q))
                if value(coeff(new,l,r),2,2)>n: break
                if q not in RESET:
                    visit(new,inner+(e,)+(0,)*q)
                q+=1
            e+=1

    q=1
    while value(coeff(matrix(q),l,r),2,2)<=n:
        if q not in RESET:
            visit(matrix(q),(0,)*q)
        q+=1
    return best

def cost(w): return len(w)+sum(w)

def better(a,b):
    if a is None: return b
    if b is None: return a
    return min(a,b,key=lambda w:(cost(w),w))

@cache
def divisors(n):
    div=set()
    for d in range(2,isqrt(n)+1):
        if n%d==0:
            div.update((d,n//d))
    return tuple(sorted(div))

@cache
def solve(n,l=0,r=0):
    if n<3 or not n&1: return None
    best=primitive(n,l,r)
    for d in divisors(n):
        for flag,q in enumerate(RESET):
            left=primitive(d,l,flag)
            if left is None: continue
            right=solve(n//d,flag,r)
            if right is None: continue
            best=better(best,left+(0,)*q+right)
    return best

def minimum(n):
    w=solve(n)
    return None if w is None else (1+cost(w),w)

def adjacency_rows(gaps):
    anchors=[0]
    for g in gaps: anchors.append(anchors[-1]+g+1)
    n=anchors[-1]+1
    back=set(zip(anchors,anchors[1:]))
    return tuple(sum(1<<j for j in range(n) if i!=j and ((i<j) != ((min(i,j),max(i,j)) in back)))
                 for i in range(n))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('target',nargs='+',type=int)
    for n in p.parse_args().target:
        z=minimum(n)
        print(n,'no first-terrace witness' if z is None else (*z,adjacency_rows(z[1])))
