from math import isqrt
pure_paths = (
 ((1,2,1,2),(64,96,72,108,81),(64,17,8,73,96,4,21,28,72,49,51,13,108,36,45,19,81)),
 ((2,1,2,1),(64,48,72,54,81),(64,17,8,1,48,33,3,97,72,49,15,10,54,27,9,40,81)),
)
Q15=(8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
short=(64,36,28,72,9,40,81)

def sq(x):
    return isqrt(x)**2==x

def guard(x,y):
    k=0; u=3*x+1
    while u%2==0:k+=1;u//=2
    assert u==y
    return k

def glue(paths):
    U=set()
    for i,P in enumerate(paths):
        if i:assert U&set(P)=={P[0]}
        else:assert not U&set(P)
        U.update(P)
    return len(U)

assert (Q15[0], Q15[7], Q15[-1]) == (8,12,9)
assert all((p+7)%16 == q for p,q in zip(Q15[::2],Q15[2::2]))
for k in range(1,12):
    lam=tuple(2**(6*k-3*i)*3**(2*i) for i in range(k))
    path=tuple(tuple(d*x for x in Q15) for d in lam)
    assert glue(path)==14*k+1
    for i,P in enumerate(path):
        assert all(sq((x+y)//(1+i%2)) for x,y in zip(P,P[1:]))
        d=lam[i]
        assert [guard(8*d-5,12*d-7),guard(12*d-7,9*d-5)]==[1,2]
        assert [guard(8*d-7,6*d-5),guard(6*d-5,9*d-7)]==[2,1]
    mu=tuple(2**(6*(k-i))*3**(4*i) for i in range(k))
    assert glue(tuple(tuple(d*x for x in short) for d in mu))==6*k+1
    assert all(sq(x+y) for x,y in zip(short,short[1:]))
    for J,C,P in pure_paths:
        assert [P.index(x) for x in C]==[0,4,8,12,16]
        paths=tuple(tuple(d*x for x in P) for d in mu)
        assert glue(paths)==16*k+1
        for d,R in zip(mu,paths):
            assert all(sq(x+y) for x,y in zip(R,R[1:]))
            X=tuple(d*x-(5 if J[i%2]==1 else 7) for i,x in enumerate(C))
            assert tuple(guard(x,y) for x,y in zip(X,X[1:]))==J
for j in range(2,24,2):
    A=2**j+3
    m=pow(2,-j//2,A)
    if m%2:m+=A
    if m==0:m=2*A
    x=(2**j*m*m-1)//A;y=(3*m*m+1)//A
    assert guard(x,y)==j and x+y==m*m and x%4==1 and y%4==3
print("all exact guards, paths, intersections, and sample even raw edges pass")

# One common square carrier, both ordered words, and their common middle.
def shared():
    P=(64,105,16,48,73,96,100,21,28,72,9,27,54,46,3,61,108,88,81)
    assert len(set(P))==len(P) and P[len(P)//2]==72
    assert all(isqrt(a+b)**2==a+b for a,b in zip(P,P[1:]))
    left,right=(64,48,72,54,81),(64,96,72,108,81)
    for k in range(1,7):
        parts=[]
        for i in range(k):
            m=(1<<(6*(k-i)))*3**(4*i)
            parts.append({m*x for x in P})
            for coords,shifts in [(left,(7,5,7,5,7)),(right,(5,7,5,7,5))]:
                assert list(map(P.index,coords))==sorted(map(P.index,coords))
                vals=[m*x-c for x,c in zip(coords,shifts)]
                for a,b,j in zip(vals,vals[1:],(2,1,2,1) if coords==left else (1,2,1,2)):
                    assert (b&1) and 3*a+1==(1<<j)*b
        for i,S in enumerate(parts):
            for j,T in enumerate(parts[i+1:],i+1):
                assert len(S&T)==int(j==i+1)
        assert len(set.union(*parts))==18*k+1
    print('common square-sum carrier and paired integer guards pass')
shared()

P=(64,105,16,33,48,52,69,100,96,73,27,9,72,
   97,3,46,54,10,39,61,108,36,45,19,81)
assert len(P)==len(set(P))==25
assert P[::4]==(64,48,96,72,54,108,81)
assert all(isqrt(a+b)**2==a+b for a,b in zip(P,P[1:]))
assert [24-P.index(x) for x in (64,48,72,54,81)]==list(
    reversed([P.index(x) for x in (64,96,72,108,81)]))
print('regular mirror ports pass')
