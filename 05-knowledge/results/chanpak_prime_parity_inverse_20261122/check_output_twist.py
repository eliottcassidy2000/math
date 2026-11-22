#!/usr/bin/env python3
from itertools import product
from math import gcd, isqrt

M = dict(A=(3,3,6,8), B=(9,9,30,32), C=(3,3,12,16))


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def state(word):
    h,r=1,2
    for letter in word:
        a,b,c,d=M[letter]
        h,r=a*h+b*r,c*h+d*r
    return h,r


def divisor_list(n):
    q=2
    ans=[]
    while q*q<=n:
        if n % q == 0:
            ans.append(q)
            n//=q
        else:
            q+=1
    if n>1:
        ans.append(n)
    return ans


def decode(N1,N2):
    if (N2-N1)%2:
        raise ValueError("odd difference")
    r,h=(N2-N1)//2,2*N1-N2
    word=[]
    while (h,r)!=(1,2):
        cand=[]
        if min(h,r)<=0:
            raise ValueError("nonpositive state")
        for letter,lo,ld,hi,hd in (("A",22,9,8,3),
                                    ("B",94,27,32,9),
                                    ("C",44,9,16,3)):
            if lo*h<=ld*r and hd*r<hi*h:
                cand.append(letter)
        if len(cand)!=1:
            raise ValueError("slope")
        l=cand[0]
        a,b,c,d=M[l]
        det=a*d-b*c
        x,y=d*h-b*r,-c*h+a*r
        if x%det or y%det or not 0<x//det<h:
            raise ValueError("inverse")
        h,r=x//det,y//det
        word.append(l)
    return ''.join(word[::-1])


def verify_depths():
    z=0
    for n in range(1,9):
        for w in map(''.join,product("ABC",repeat=n)):
            h,r=state(w)
            N=[h+2**j*r for j in range(1,4)]
            require(decode(*N[:2])==w, "pair inverse")
            require(all(gcd(x,6)==1 for x in N), "bad prime")
            require(all(gcd(a,b)==1 for i,a in enumerate(N) for b in N[i+1:]),"coprimality")
            # v2(det)=#A+#B+2#C; v3(det)=#A+2#B+#C
            delta=(w.count('A')+w.count('B')+2*w.count('C')
                   +w.count('A')+2*w.count('B')+w.count('C'))
            require(delta%2==(w.count('B')+w.count('C'))%2,"det parity")
            z+=1
    return z


if __name__=='__main__':
    for words, modulus, checks in [
        (('CAABBBB','CBBBBAA'),385,[[19,1488116341],[27107951539]]),
        (('AAB','BAA'),11,[[23,1009],[22811]])]:
        states=[state(w) for w in words]
        require(len({tuple(w.count(c) for c in "ABC") for w in words})==1, "histogram")
        require(len({tuple(x%modulus for x in s) for s in states})==1,"collision modulus")
        for (h,r),f in zip(states,checks):
            require(divisor_list(h+2*r)==f,"literal factorization")
        print("aliased modulus",modulus,"outputs",[h+2*r for h,r in states])
    print("physical word pairs inverted:",verify_depths())
