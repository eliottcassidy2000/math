#!/usr/bin/env python3
from check_permutation import counts, validate, need

def suffix(old,b):
    # Last vertex is the root. Relabel it last after inserting the block.
    r=len(old)-1
    k=r+b
    G=tuple((a&((1<<r)-1))|(((1<<b)-1)<<r)|((a>>r&1)<<k) for a in old[:-1])
    B=tuple(((1<<(b-i))-1)<<(r+i+1) for i in range(b)) # later block vertices and root
    return G+B+(old[-1],),k

banks=[
 (8943,((45,2934,(2046,4092,2040,4080,1376,3776,1408,3856,3616,3152,0,1109)),
        (45,2933,(2046,4092,2040,4080,1248,1216,3968,3840,3632,3120,0,1077)),
        (45,2932,(2046,4092,2040,4080,1248,3264,3968,1792,3632,3120,0,1173)))),
 (189,((37,38,(958,1020,952,944,992,448,909,768,0,288)),
       (37,37,(30,108,24,48,98,69,13)),
       (49,42,(190,252,152,240,32,68,21,112))))
]
G=(30,12,56,16,34,11)
swap=[]
for A,K in ((40,8),(48,7),(20,6)):
    swap.append((29,K,tuple(v+((1<<6) if not A>>i&1 else 0) for i,v in enumerate(G))+(A,)))
banks.append((81,tuple(swap)))
for F,bank in banks:
    for g,K,a in bank:
        validate(a)
        h,e,o,d=counts(a,len(a)-1)
        need((d,e+o-d,2*o-d-h)==(g,F,2*K),"three signature coordinates")
        print("signature",g,F,K)
        for b in (1,2,3):
            a2,root=suffix(a,b)
            validate(a2)
            want=(F*(1<<b)+g-2*K,F*(1<<(b-1)),F*(1<<(b-1))+g,g)
            need(counts(a2,root)==want,"suffix subset DP")
            print(b,*want)
print("All target and swap templates verified.")

bits="11101001000";r=len(bits)
apex=tuple(sum(1<<j for j in range(i+1,r))|((b=="0")<<r) for i,b in enumerate(bits))+(sum((b=="1")<<i for i,b in enumerate(bits)),)
need(counts(apex,r)==(1729,932,933,1),"literal 1729")
print("Literal 1729 verified outside the mirror architecture.")
