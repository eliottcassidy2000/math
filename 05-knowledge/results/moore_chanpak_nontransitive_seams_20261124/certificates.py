"""Exact profile, run matrices, rigidity and target lattice for three witnesses."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations, product
from math import gcd, isqrt
from bounded_template_profile import profile

BASE=(122,100,25,82,34,76,20)
FLIP=(*BASE[:5],BASE[5]^64,BASE[6]^32)
PRIME_CASE=(
 (258,388,481,247,263,211,19,81,232),(2,132,97,119,7,83,19,125),
 (1139,22912,218313),(137,2888,27513),85587188881,85587188879)
CASES=(
 ((66,284,409,209,33,463,150,275,89),(62,124,56,112,32,64,21),
  (1739,28384,244641),(31,644,5541),19116233339,19116233337),
 ((338,444,201,369,356,325,130,313,68),(126,452,504,338,34,202,400,281,49),
  (1747,30344,265161),(1753,30068,262827),984193130605,984193130607),
 ((94,476,344,464,416,463,144,261,65),(446,476,312,336,224,202,389,268,48),
  (1001,19288,185787),(1061,21328,205407),544332106255,544332106257))

def check_tournament(T):
    assert all(not (T[i]>>i&1) and T[i]<(1<<len(T)) for i in range(len(T)))
    assert all(((T[i]>>j)^(T[j]>>i))&1 for i,j in combinations(range(len(T)),2))

def auto(T):
    cells={}
    for i,w in enumerate(T):cells.setdefault(w.bit_count(),[]).append(i)
    keys=[i for C in cells.values() for i in C]
    return sum(all(((T[i]>>j)^(T[p[i]]>>p[j]))&1==0 for i,j in combinations(range(len(T)),2))
        for P in product(*(permutations(C) for C in cells.values()))
        for p in [dict(zip(keys,sum((tuple(q) for q in P),())))])

def matrix(T):
    W=sum(1<<i for i in (0,1,2,5,6))
    @cache
    def go(U,last,a,b):
        out=Counter()
        if not U and a and b:out[a,b]=1
        for i in range(7):
            if last>=0 and not (T[last]>>i&1):continue
            if i==3 and a<3:out.update(go(U,i,a+1,b))
            if i==4 and b<3:out.update(go(U,i,a,b+1))
            if U>>i&1:out.update(go(U^(1<<i),i,a,b))
        return out
    cnt=go(W,-1,0,0)
    return tuple(tuple(cnt[i,j] for j in (1,2,3)) for i in (1,2,3))

def runs(t):
    h,q,r=t
    return h,q-2*h,r-3*q+3*h

def mv(M,y):
    return tuple(sum(a*b for a,b in zip(row,y)) for row in M)

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def det(u,v):
    return u[0]*v[1]-u[1]*v[0]

def lattice(M,Mp,y,N,D):
    """Exact membership in the integer image (not requiring nonnegative run counts)."""
    w=mv(M,y); d=tuple(b-a for a,b in zip(w,mv(Mp,y)))
    columns=tuple(zip(w,d))
    g=gcd(*(det(a,b) for a,b in combinations(columns,2)))
    assert g
    return all(det(c,(N,D))%g==0 for c in columns),g

def order_ok(A,n):
    h,q,r=A
    return (2*n+4)*q*(n+1) <= (4*r+5*h)*(n+1)+5*q+2*(h+1)

def max_partner(A):
    n=1
    while order_ok(A,n):n+=1
    return n-1

def prime(n):
    if n<2:return False
    return all(n%i for i in range(2,isqrt(n)+1))

def main():
    M,Mp=matrix(BASE),matrix(FLIP)
    assert M==((109,117,32),(94,117,40),(24,36,16))
    assert Mp==((111,100,24),(111,117,36),(32,40,16))
    for X,Y,a,b,H,J in (*CASES,PRIME_CASE):
        check_tournament(X);check_tournament(Y)
        assert tuple(profile(X,k) for k in (1,2,3))==a
        assert tuple(profile(Y,k) for k in (1,2,3))==b
        assert (auto(X),auto(Y))==(1,1)
        assert (gcd(*a),gcd(*b))==(1,1)
        xr,yr=runs(a),runs(b)
        assert dot(xr,mv(M,yr))==H and dot(xr,mv(Mp,yr))==J
        assert lattice(M,Mp,yr,H,J-H)[0]
        assert max_partner(a)>=len(Y) and max_partner(b)>=len(X)
        if H==PRIME_CASE[-2]:
            assert prime(H) and J==37*167*2099*6599
        else:
            assert not prime(H) and not prime(J)
        print(len(X),len(Y),a,b,H,J,'partner bounds',max_partner(a),max_partner(b),
              'fixed-Y lattice index',lattice(M,Mp,yr,H,J-H)[1])

if __name__=='__main__':main()
