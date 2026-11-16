#!/usr/bin/env python3
"""Exact count for the g=3 decoder, with a linear scan of each contact word."""
from itertools import permutations
from check import strong_g3_contacts,state,from_words,require

def tile(seq,arcs):
    n=len(seq);H=[0]*(1<<n);P=[]
    for S in range(1<<n):
        for p in permutations([i for i in range(n) if S>>i&1]):
            if all((a,b) in arcs for a,b in zip(p,p[1:])):
                H[S]+=1
                if p:P.append((S,p))
    menu=[[0,0] for _ in range(4)]
    for S,p in P:
        for j in range(4):
            if (not(j&1) or seq[p[0]]) and (not(j&2) or not seq[p[-1]]):
                menu[j][len(p)&1]+=H[((1<<n)-1)^S]
    return H[-1],tuple(tuple(x) for x in menu)

TRI=tuple(tile((1,)*k+(0,)*(3-k),{(0,1),(1,2),(2,0)}) for k in range(4))
VERT=tuple(tile((b,),set()) for b in range(2))

def read(tiles):
    z=1; A=[0,0]; B=[0,0]
    for h,(wild,heads,tails,both) in tiles:
        def conv(w,i):return A[i]*w[0]+A[i^1]*w[1]
        A,B=([h*A[i]+z*heads[i]+conv(wild,i) for i in range(2)],
             [h*B[i]+z*both[i]+conv(tails,i) for i in range(2)])
        z*=h
    return z+2*B[0]

def H3(prefix,k,suffix):
    return read([*(VERT[b] for b in prefix),TRI[k],*(VERT[b] for b in suffix)])

def H1(J):
    return read([VERT[int(b)] for b in bin(J)[2:]])

def offset(prefix,k,suffix):
    U=0
    for c in prefix:U=2*U+c
    L=U+U//2
    e=suffix[0] if suffix else 1
    return L+e if k==1 else -L-1-e if k==2 else 0

def verify_offsets(limit):
    for J in range(1,limit+1):
        h=H1(J)
        D={p:set() for p in (-3,-2,-1)}
        for p in D:
            for w in strong_g3_contacts(J,p):
                d=offset(*w)
                require(H3(*w)==3*h+2*d,"signed offset")
                D[p].add(d)
            require(all(not({d,d+1,d+2}<=D[p]) for d in D[p]),"three-digit control")
        if J%2==0:
            # The g=1 terminal case can additionally supply offset 0.
            D[-3].add(0)
            require(all(not({d,d+1,d+2}<=D[-3]) for d in D[-3]),"padded control")

def study(limit):
    three=None; pair=None; mixed=None
    for J in range(1,limit+1):
        for D in (-3,-2,-1):
            C=strong_g3_contacts(J,D)
            menu={H3(*r):r for r in C}
            if not mixed and len({h%4 for h in menu})>1:mixed=(J,D,menu)
            for h in menu:
                if not pair and h+2 in menu:pair=(J,D,h,menu)
                if h+2 in menu and h+4 in menu:
                    three=(J,D,h,menu)
                    return mixed,pair,three
    return mixed,pair,three

if __name__=="__main__":
    for J in range(1,49):
        for D in (-3,-2,-1):
            for r in strong_g3_contacts(J,D):
                require(H3(*r)==state(*from_words(*r))[0],"transfer mismatch")
    verify_offsets(4096)
    print("closed signed offsets agree with the independent component scan")
