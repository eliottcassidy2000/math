#!/usr/bin/env python3
"""All nonconstant bases through length 8, repetitions 2..5, every cyclic unequal swap."""
from itertools import product
from math import gcd
import json

def need(ok,message):
    if not ok: raise RuntimeError(message)

def carry(w):
    B=0
    for i,b in enumerate(w): B=(3 if b else 1)*B+b*2**i
    return B

count=crossing=0
examples=[]
for n in range(2,9):
    for v in product((0,1),repeat=n):
        p=sum(v)
        if p in (0,n):continue
        D0=2**n-3**p
        for d in range(2,6):
            c=list(v)*d
            q=len(c)
            D=2**q-3**(d*p)
            H=sum(2**(n*(d-1-j))*3**(p*j) for j in range(d))
            need(D==H*D0 and carry(c)==H*carry(v),'repetition law')
            need(H>1 and gcd(H,6)==1,'clock unit')
            for i in range(q):
                j=(i+1)%q
                if c[i]==c[j]:continue
                w=c.copy();w[i],w[j]=w[j],w[i]
                B=carry(w);g=gcd(abs(D),B)
                need(gcd(H,B)==1,'repetition factor not destroyed')
                need(g==gcd(abs(D0),B),'exact gcd reduction')
                need(g<abs(D),'integer cycle')
                if g==abs(D0) and g>1 and len(examples)<5:
                    examples.append({'base':''.join(map(str,v)),'d':d,'swap_start':i,'edited':''.join(map(str,w)),'D0':D0,'H':H,'g':g})
                count+=1;crossing+=int(i==q-1)
print(json.dumps({'status':'PASS','universe':'Every nonconstant binary base length2..8, repetition count2..5, every unequal cyclic adjacent swap occurrence','swap_checks':count,'crossing_checks':crossing,'retained_base_clock_examples':examples},indent=2))
