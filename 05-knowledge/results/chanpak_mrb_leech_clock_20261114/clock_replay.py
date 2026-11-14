from collections import deque

QS = (4,2,5,3,7)
K = [(1,1,1),(1,2,2),(3,3,4)]
for i in range(3,14):
    K.append(tuple(sum(K[j][a] for j in range(i-3,i)) for a in range(3)))

def step(v,q,y,p):
    a,b=v
    e,f,g=K[q]
    return ((e*a+f*b)%p,((f*a+g*b)*y)%p)

def reset_orbit(p):
    seeds={(1,1)}
    gens=((1,1,1,1,1),(1,1,1,1,2),(1,1,1,2,1))
    queue=deque(seeds)
    states=set(seeds)
    while queue:
        u=queue.popleft()
        for w in gens:
            v=u
            for q,y in zip(QS,w):
                v=step(v,q,y,p)
            if v not in states:
                states.add(v)
                queue.append(v)
    return states

def whole_cycle(p):
    # All positive exponent options, including residue 1 at full period.
    assert p%2==1
    Y={1};y=2%p
    while y!=1:
        Y.add(y)
        y=2*y%p
    S={(1,y) for y in Y}
    out=[]
    for q in QS:
        S={step(v,q,y,p) for v in S for y in Y}
        out.append(len(S))
    return S,out

if __name__=='__main__':
    from math import prod
    ds=[K[q][0]*K[q][2]-K[q][1]**2 for q in QS+(13,)]
    assert ds==[-4,3,8,-1,-7,0]
    assert prod(ds[:-1])==-672
    for p in (5,11,13,17,19,23,29,31,37,47,73,127):
        S=reset_orbit(p)
        assert len(S)==p*p-1,(p,len(S))
    print('3 explicit full-cycle generators reach all nonzero states at tested good and exceptional primes')
    for p in (3,7,21,9,27,49,63,147):
        S,history=whole_cycle(p)
        print(p,history,len({(a+b)%p for a,b in S}))
        assert {(a+b)%p for a,b in S}==set(range(p))

def matmul(A,B,p):
    a,b,c,d=A
    e,f,g,h=B
    return tuple(x%p for x in (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h))

def projective_order(p):
    def norm(A):
        z=next(x for x in A if x)
        zi=pow(z,-1,p)
        return tuple(x*zi%p for x in A)
    def inv(A):
        a,b,c,d=A
        return norm(tuple(x%p for x in (d,-b,-c,a)))
    (e,f,g)=K[7]
    B=(e%p,f%p,f%p,g%p); I=(1,0,0,1); D=(1,0,0,2)
    R=norm(matmul(matmul(B,D,p),inv(B),p))
    C0=I
    for q in QS:
        e,f,g=K[q]
        C0=matmul((e,f,f,g),C0,p)
    states={I};todo=deque(states)
    while todo:
        A=todo.popleft()
        for C in (D,R,C0):
            B=norm(matmul(C,A,p))
            if B not in states:
                states.add(B)
                todo.append(B)
    return len(states)

if __name__=='__main__':
    for p in (5,17,19,31):
        h=projective_order(p); l=p*(p*p-1)//2
        assert h in (l,2*l)
        print(f'projective monodromy p={p}: {h} elements; PSL size {l}')

def clock_matrices(m):
    G=[]
    for ys in ((1,1,1,1,1),(1,1,1,1,2),(1,1,1,2,1)):
        B=(1,0,0,1)
        for q,y in zip(QS,ys):
            e,f,g=K[q]
            B=matmul((e,f,y*f,y*g),B,m)
        G.append(B)
    return G

def power(A,j,m):
    X=(1,0,0,1)
    while j:
        if j&1:
            X=matmul(X,A,m)
        j>>=1
        A=matmul(A,A,m)
    return X

def terminal_discharge():
    from itertools import product
    M=63
    A=clock_matrices(M)[0]
    E=power(A,6,M)
    assert E==(12,51,53,52) and matmul(E,E,M)==E
    for p,lambdas in ((3,(2,1,2)),(7,(3,5,5))):
        e=tuple(x%p for x in E)
        for g,h in zip(clock_matrices(p),lambdas):
            assert matmul(matmul(e,g,p),e,p)==tuple(h*x%p for x in e)
    def P(w,m):
        a,b,c,d=E
        v=((a+b*w[0])%m,(c+d*w[0])%m)
        for q,y in zip(QS,w[1:]):
            v=step(v,q,y,m)
        return sum(v)%m
    witnessed={}
    for w in product((1,2,4,8,16,32),repeat=6):
        N=P(w,M)
        if N in witnessed: continue
        ds=[(P(w[:i]+(w[i]+1,)+w[i+1:],21)-N)%21 for i in range(6)]
        ix=next((i for i,d in enumerate(ds) if d%3),None)
        iy=next((i for i,d in enumerate(ds) if d%7),None)
        if ix is not None and iy is not None:
            witnessed[N]=(w,ix,iy)
    assert set(witnessed)==set(range(M))
    print('the idempotent and compressed units check; all 63 terminal residues carry both bad-prime derivative guards')
    return witnessed

if __name__=='__main__':
    terminal_discharge()
