from itertools import product

Q=(4,2,5,3,7)
table=[(1,1,1),(1,2,2),(3,3,4)]
while len(table)<=max(Q):
    table.append(tuple(map(sum,zip(*table[-3:]))))

def mul(a,b,m):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in (0,1))%m
                 for i,j in product(range(2),repeat=2))

def clock(Y,m):
    a=(1,0,0,1)
    for q,y in zip(Q,Y):
        e,f,g=table[q]
        a=mul((e,f,y*f,y*g),a,m)
    return a

def power(a,n,m):
    r=(1,0,0,1)
    while n:
        if n&1:r=mul(r,a,m)
        a=mul(a,a,m);n//=2
    return r

A=clock([1]*5,63)
E=power(A,6,63)
assert E==(12,51,53,52) and mul(E,E,63)==E
for r in range(6):
    Z=mul(E,mul(clock([2**r]+[1]*4,63),E,63),63)
    assert all(any(x%p for x in Z) for p in (3,7))

gamma=9*3*17*5*57
for k in range(1,12):
    m=1<<k
    baseline=clock([0]*5,m)
    Einv=power(baseline,1<<k,m)
    assert Einv==(1,11*pow(9,-1,m)%m,0,0)
    for r in range(1,k):
        t=(1+11*pow(9,-1,m)*(1<<r))%m
        U=mul(Einv,mul(clock([(1<<r)]+[0]*4,m),Einv,m),m)
        assert U[0]==gamma*t%m and U[2]==0
    for n in range(1,m,2):
        need=n*pow(pow(gamma,k,m),-1,m)%m
        u=1
        for r in range(1,k):
            if (u^need)>>r&1:
                u=u*(1+11*pow(9,-1,m)*(1<<r))%m
        assert u==need
print('six odd unit-preservation tests; all binary suffix repairs through modulus 2048')
