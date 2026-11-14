from math import isqrt

def row(k):
    u,v=47,56
    for _ in range(k):
        u,v=u+4*v,2*u+4*v
    return u,v

def guard(q,a,b):
    return (a+2*b)%3!=0 if q==1 else q==2 and b%3!=0

def step(a,b,q,r):
    assert q in (1,2) and r>=1
    if q==1:return a+2*b,(a+b)<<(r+1)
    return 3*(a+b),(3*a+4*b)<<r

def decompose(a,b):
    """Exact inverse of the guarded q=1,2 language, or None."""
    letters=[]
    while a>0 and a&1 and b>=2*a and not (b&1):
        if a==1:
            if b&(b-1): return None
            return (b.bit_length()-1,tuple(reversed(letters)))
        val=(b&-b).bit_length()-1
        q=1 if a%3 else 2
        r=val-(q==1)
        if r<1:return None
        d=b>>val
        if q==1:
            x,y=2*d-a,a-d
        else:
            u=a//3
            x,y=4*u-d,d-3*u
        if x<=0 or x>=a or x%2!=1 or y<2*x or y&1 or not guard(q,x,y):
            return None
        assert step(x,y,q,r)==(a,b)
        letters.append((q,r))
        a,b=x,y
    return None

def choose(k,seed,letters,q):
    a,b=1,1<<seed
    for j,r in letters:
        assert guard(j,a,b)
        a,b=step(a,b,j,r)
    assert guard(q,a,b)
    A,D=step(a,b,q,1)
    D//=2
    u,v=row(k)
    W=u+2*v
    if 2*A>W:return None
    r=1
    while (D<<r)<3*W:r+=1
    B=D<<r
    assert 3*W<=B<6*W
    p=u*A+v*B
    return p,A,B,r

def greedy(p):
    assert p>0 and p&1
    k=0
    u,v=47,56
    if p<(u+2*v)**2:return None
    while True:
        a,b=u+4*v,2*u+4*v
        if (a+2*b)**2>p:break
        u,v=a,b;k+=1
    a=p*pow(u,-1,2*v)%(2*v)
    b=(p-u*a)//v
    if b<2*a or b&1:return None
    return k,a,b

def exact(p):
    z=greedy(p)
    if z is None:return None
    k,a,b=z
    z=decompose(a,b)
    if z is None:return None
    seed,letters=z
    w=(seed,)
    for q,r in letters+((1,1),)*k+((7,2),):
        w+=(0,)*q+(r,)
    return (k,a,b,w)

def packing(gaps):
    cum=[0]
    for g in gaps:cum.append(cum[-1]+g)
    P=[1]
    for j,x in enumerate(cum):
        P.append(P[-1]+sum((P[i]<<(x-cum[i])) if x>cum[i] else
                           2*P[i]*int((j-i)%2==0) for i in range(j)))
    return P[-1]

def prime(p):
    return p>=2 and all(p%i for i in range(2,isqrt(p)+1))

if __name__=='__main__':
    # Distinct changed and mixed seeds, verified via literal gap packing.
    examples={1:None,2:None,3:None,4:None}
    menu=[(1,(),q) for q in (1,2)]
    for w in range(1,9):
        for r in range(1,7):
            a,b=1,2**w
            for q1 in (1,2):
                if not guard(q1,a,b):continue
                x,y=step(a,b,q1,r)
                for q2 in (1,2):
                    if guard(q2,x,y):menu.append((w,((q1,r),),q2))
    for k in examples:
        for seed,letters,q in menu:
            if len(letters)<1:continue
            z=choose(k,seed,letters,q)
            if z is None:continue
            p,a,b,r=z
            if not prime(p):continue
            c=exact(p)
            assert c[:3]==(k,a,b)
            assert packing(c[3])==7*p
            examples[k]=(p,a,b,seed,letters+( (q,r), ),1+len(c[3])+sum(c[3]))
            break
    for k,x in examples.items():print(k,x)
