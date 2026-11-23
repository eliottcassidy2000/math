from math import gcd,isqrt
from itertools import combinations,permutations

def tri(a,b):
    assert a%2 and b%2 and 0<a<b and gcd(a,b)==1
    return (a*b,(b*b-a*a)//2,(b*b+a*a)//2)

def A(a,b):return a,b+2*a
def B(a,b):return b,2*b+a
def C(a,b):return b,2*b-a

def height_word(w):
    H,P,e=1,0,0
    for bit in w:
        if bit=="O":P,e=2*P+1,1
        else:H,P,e=H+P+e,2*P,0
    return H

def pairs(n):return list(combinations(range(n),2))
def rows(n,mask):
    out=[0]*n
    for i,(a,b) in enumerate(pairs(n)):
        if mask>>i&1:out[a]|=1<<b
        else:out[b]|=1<<a
    return out
def apex(w):
    n=len(w)+1
    return rows(n,sum(1<<i for i,(a,b) in enumerate(pairs(n))
                      if b<n-1 or w[a]=="I"))
def count(R):
    n=len(R)
    return sum(all(R[p[i]]>>p[i+1]&1 for i in range(n-1))
               for p in permutations(range(n)))
def canon(R):
    n=len(R)
    return min(sum((R[p[a]]>>p[b]&1)<<i for i,(a,b) in enumerate(pairs(n)))
               for p in permutations(range(n)))
def induce(R,I):
    return [sum(((R[u]>>v)&1)<<j for j,v in enumerate(I)) for u in I]
def strong(R):
    n=len(R)
    for i in range(n):
        S=1<<i
        for _ in range(n):
            S|=__import__("functools").reduce(int.__or__,(r for j,r in enumerate(R) if S>>j&1),0)
        if S!=(1<<n)-1:return False
    return True
def split(M):
    return [(d,M//d) for d in range(1,isqrt(M)+1) if M%d==0 and gcd(d,M//d)==1]

if __name__=="__main__":
    assert tri(*B(1,3))==(21,20,29) and tri(*C(1,11))==(231,160,281)
    assert [height_word(x) for x in ("OII","OIOIO","OOOI")]==[5,9,9]
    for j in range(1,11):
        u=height_word("OI"+"OIO"*j);v=height_word("O"+"OOI"*j)
        assert 49*(u-3)==38*(8**j-1)+28*j
        assert 49*(v-1)==52*(8**j-1)+28*j
        assert 7*(v-u)==2*(8**j-8)
    T_B,T_C=apex("OIOI"),apex("OOOI")
    assert count(T_B)==count(T_C)==9 and strong(T_B) and strong(T_C)
    cB,cC=canon(T_B),canon(T_C)
    assert (cB,cC)==(10,8)
    lower=703
    x,y=[1<<pairs(5).index(e) for e in ((0,4),(0,1))]
    rr=[rows(5,r) for r in (lower,lower^x,lower^y,lower^x^y)]
    assert [count(r) for r in rr]==[5,9,9,9]
    aa=apex("OII")
    assert any(canon(induce(rr[0],I))==canon(aa) for I in combinations(range(5),4))
    assert (canon(rr[1]),canon(rr[2]))==(cC,cB)
    E,F,G=[1<<pairs(6).index(e) for e in ((3,4),(0,5),(2,4))]
    for g,table in [(0,[5,9,19,23]),(G,[9,9,29,23])]:
        assert [count(rows(6,4^g^e^f)) for f in (0,F) for e in (0,E)]==table
    for x,c in ((E,cC),(G,cB)):
        R=rows(6,4^x)
        ks={canon(si) for I in combinations(range(6),5)
            for si in (induce(R,I),) if strong(si)}
        assert ks=={c}
    assert [tri(*p) for p in split(21)]==[(21,220,221),(21,20,29)]
    assert len(split(231))==4
    assert [(b*b-a*a,a*a+b*b) for a,b in split(66)]==[
        (4355,4357),(1085,1093),(475,493),(85,157)]
    print("Berggren factors, divisor-leg quartet, nonisomorphic twins and physical cube PASS")
