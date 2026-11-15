#!/usr/bin/env python3
from fractions import Fraction

P,Q='1101010010','1010'
W,V='0001010101','0000000001'
S,Z='0110','0110011100'

def words(k):
    return P+W*k+S, Q+V*k+Z

def state(word,seed=(1,0,0)):
    n,s,d=seed
    for b in map(int,word):
        n,s,d=n+(1-b)*(s+d),2*s+b,b
    return n,s,d

def coeff(p,w,s):
    q=1<<len(w)
    n,u,_=state(w,(0,0,0))
    t=sum(1<<i for i,b in enumerate(w) if b=='0')
    v=sum(1<<i for i,b in enumerate(s) if b=='0')
    return (Fraction(int(p,2))+Fraction(u,q-1))*(v+Fraction(t,q-1)), n+1-Fraction(t*u,q-1)

def adjacent(k):
    A,B=words(k)
    n=len(A)
    rows=[0]*(n+2)
    for i in range(n):
        for j in range(i+1,n):
            rows[i]|=1<<j
        for apex,c in ((n,A[i]),(n+1,B[i])):
            if c=='1': rows[apex]|=1<<i
            else: rows[i]|=1<<apex
    rows[n]|=1<<(n+1)
    return rows

def compare_old(k):
    a,b=adjacent(k),adjacent(k+1)
    image=[i if i<10 else i+10 for i in range(len(a))]
    for i in range(len(a)):
        for j in range(len(a)):
            assert (a[i]>>j)&1==(b[image[i]]>>image[j])&1

if __name__=='__main__':
    C,D=coeff(P,W,S)
    E,F=coeff(Q,V,Z)
    assert C==E==Fraction(8305014250,1046529)
    assert D-F==2
    for k in range(121):
        a,b=map(lambda x: state(x)[0],words(k))
        assert a==8259+C*((1<<10*k)-1)+D*k
        assert b==8259+C*((1<<10*k)-1)+F*k
        assert a-b==2*k
    for k in range(1,8): compare_old(k)
    print('exact states (121 heights), 7 literal nesting controls passed')

