from replay_power_three import ends,validate
from itertools import product

def g(w):
    if w<=0:return 1
    a,b=divmod(w,2)
    return g(a)+(0 if b else a+a%2)

def apex(w):
    n=w.bit_length()
    F=(1<<n)-1
    a=tuple(F^((1<<(i+1))-1) for i in range(n))
    U=sum((w>>(n-1-i)&1)<<i for i in range(n))
    return tuple(t|(int(not(U>>i&1))<<n)for i,t in enumerate(a))+(U,),n

def converse(R):
    F=(1<<len(R))-1
    return tuple(F^r^(1<<i) for i,r in enumerate(R))

def state(R,u):
    validate(R)
    F=(1<<len(R))-1
    h=[sum(a)for a in ends(R)]
    h[0]=1
    st=ends(converse(R))
    e=[0,0]
    for q in range(1,F+1):
        e[q.bit_count()%2]+=st[q][u]*h[F^q]
    return (h[F],*e)

def glue(A,u,B,w,up):
    n,m=len(A),len(B)
    out=[r|(((1<<m)-1)<<n) for r in A]+[r<<n for r in B]
    out[u]^=1<<(n+w)
    out[n+w]|=1<<u
    return tuple(out),n+up

def words(s):
    n=int(len(s)**0.5)
    return tuple(sum(int(s[i*n+j])<<j for j in range(n))for i in range(n))

def mult(M,V):
    return tuple(sum(M[i*3+j]*V[j]for j in range(3))for i in range(3))

C=((0,1),0,1,(1,2,2,1,1,0,1,2,1))
S=(words("0001100011000110"),1,2,(5,10,8,4,5,4,5,8,5))
J=((0,),0,0,(1,2,0,0,0,1,1,1,0))
R=(words("000100110"),0,2,(1,4,4,2,3,2,2,4,3))

def norm(W):return g(W),(W+1)//2,1+W//2

if __name__=='__main__':
    v=(2,-2,1);w=(2,-1,0)
    assert mult(S[-1],v)==tuple(-t for t in v)
    assert mult(J[-1],v)==tuple(-t for t in w)
    assert mult(R[-1],w)==tuple(-t for t in w)
    for r in [2,3]:
        for off in (-2,1):
            W=(3<<r)+off
            A,u=apex(W)
            V=norm(W)
            assert state(A,u)==V
            for B,ent,out,M in (C,S,J,R):
                A,u=glue(A,u,B,ent,out)
                V=mult(M,V)
                assert state(A,u)==V
    for w,r in product(range(1,4),range(2,5)):
        a=(1<<(w+1))-1
        B=tuple(((1<<(w+1))-1)^((1<<(i+1))-1) for i in range(w+1))
        for off,d in ((-2,1),(1,3)):
            A,u=apex((a<<r)+off)
            T,t=glue(A,u,B,w,w)
            assert sum(ends(T)[-1])==((a*(a+2))<<(r-1))+d
    k,r=3,2
    y=(27<<k)-1;z=(y-6)//19;m=2*z+1
    B,t=apex(y)
    for off,d in ((-2,0),(1,2)):
        A,u=apex((m<<r)+off)
        T,_=glue(A,u,converse(B),t,t)
        assert sum(ends(T)[-1])==(m*(19*m+14)<<(r-1))+19*g(z)+12+d
    print('PASS literal physical reader and return samples.')
