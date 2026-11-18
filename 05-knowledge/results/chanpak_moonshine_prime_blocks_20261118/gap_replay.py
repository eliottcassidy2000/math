"""Exact certificate for a variable unmarked ordinal gap."""
from random import Random
from fractions import Fraction
from frame_replay import join, add, cut
from padding_replay import slots, rooted
from replay import profile
from prime_support_replay import A as L, B as R

def coefficients(A,I,B,J):
    a,b=([x[-1] for x in slots(T,3)] for T in (A,B))
    X,Y=rooted(A,I),rooted(B,J)
    c=[X[q][0]*b[q]-Y[q][1]*a[q] for q in range(3)]
    d=[Y[q][0]*a[q+1]-X[q][1]*b[q+1] for q in range(3)]
    return (c[1],c[2]+d[1],d[2]), (c[1]+d[0],d[1])

def marked_key(A,I):
    _,h,b,g=[x[-1] for x in slots(A,3)]
    (_,s0),(e1,s1),(e2,s2),_=rooted(A,I)
    return (Fraction(e1,h),Fraction(s2,g),Fraction(e2-s1,b)),(Fraction(s0,h),Fraction(s1,b))

def build(A,I,Q,B,J):
    D=join(A,join(Q,B))
    T=add(D,I|(J<<(len(A)+len(Q))))
    return T,cut(T,1<<(len(T)-1))

def rand_t(n,r):
    A=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if r.getrandbits(1):A[i]|=1<<j
            else:A[j]|=1<<i
    return tuple(A)

def run():
    r=Random(725958)
    for _ in range(100):
        A,Q,B=(rand_t(r.randrange(1,5),r) for j in range(3))
        I,J=r.randrange(1<<len(A)),r.randrange(1<<len(B))
        u,v=coefficients(A,I,B,J)
        T,Tp=build(A,I,Q,B,J)
        h,b=profile(T);hp,bp=profile(Tp)
        q=[x[-1] for x in slots(Q,3)]
        assert bp-b==2*sum(k*z for k,z in zip(u,q[1:]))
        assert hp-h==sum(k*z for k,z in zip(v,q[1:]))
    print('literal subset DP agrees on 100 marked flank/gap choices')

    A=join(L,(0,));B=join((0,),R);I,J=1<<len(L),1
    u,v=coefficients(A,I,B,J)
    c,d,q=u
    assert c*5+q*13//3==0 and q%3==0
    assert (c,d,q)==(-134784,-165888,155520)
    assert sum(k*z for k,z in zip(u,[1,2,3]))==0
    print('347/349 middle-gap defect coefficients',u,'Hamiltonian shift coefficients',v)

    W=(20,5,24,19,2); weights=(-9,-5,8,2,4)
    xyz=[]
    for i in range(len(W)):
        (_,s0),(e1,s1),(e2,s2),_=rooted(W,1<<i)
        xyz.append((e1,s2,e2-s1,s0,s1))
    w=tuple(sum(a*b for a,b in zip(weights,col)) for col in zip(*xyz))
    assert w==(0,0,0,12,15)
    print('weighted, NOT physical, collision changes hidden readings',w[3:])

if __name__=='__main__':run()
