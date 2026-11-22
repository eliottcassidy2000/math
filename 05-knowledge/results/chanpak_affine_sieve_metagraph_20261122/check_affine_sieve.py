#!/usr/bin/env python3
"""Finite checks for the fixed-alphabet transplantation, not the expander theorem."""
from collections import deque
from fractions import Fraction as Q
from itertools import product
from math import gcd

AA=((3,3),(6,8))
BB=((9,9),(30,32))
CC=((3,3),(12,16))
DD=((1,0),(0,2))
BASIS=((3,0),(8,1))
ALPH=((3,2),(5,8),(4,2))


def need(b, label):
    if not b:
        raise ArithmeticError(label)


def mm(a,b):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)


def mv(a,x):
    return tuple(sum(t*z for t,z in zip(row,x)) for row in a)


def inv(a):
    (x,y),(z,w)=a
    det=x*w-y*z
    return ((Q(w,det),Q(-y,det)),(Q(-z,det),Q(x,det)))


def add(T,k,u):
    """Literal inherited migration: old marked terminal is the last label."""
    n=len(T)
    rows=list(T)
    rows[-1]|=1<<(n+k)
    bits=f'{u:0{k}b}'
    for j,b in enumerate(bits):
        r=((1<<n)-1)|(((1<<k)-(1<<(j+1)))<<n)
        if b=='0':
            r|=1<<(n+k)
        rows.append(r)
    rows.append(((1<<(n-1))-1)|(int(bits[::-1],2)<<n))
    return tuple(rows)


def hp(T):
    n=len(T)
    states=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        states[1<<v][v]=1
    for mask,row in enumerate(states):
        for v,x in enumerate(row):
            if x:
                dest=T[v]&~mask
                while dest:
                    bit=dest&-dest
                    dest^=bit
                    states[mask|bit][bit.bit_length()-1]+=x
    return sum(states[-1])


def main():
    G=mm(BB,inv(AA))
    J=mm(mm(mm(DD,G),inv(DD)),inv(G))
    need(mm(DD,AA)==CC and G==((3,0),(8,1)), 'relative matrix')
    need(mm(mm(inv(BASIS),J),BASIS)==((1,0),(8,1)), 'lower shear')
    own=mm(mm(AA,J),inv(AA))
    need(mm(mm(inv(BASIS),own),BASIS)==((1,Q(-4,3)),(0,1)), 'upper shear')
    literal=0
    for w in product(range(3), repeat=2):
        for length in (1,2):
            if length==1 and w[1]:
                continue
            T=(2,0)
            x=(1,2)
            for j in w[:length]:
                T=add(T,*ALPH[j])
                x=mv((AA,BB,CC)[j],x)
            out=hp(add(T,1,1))
            need(out==x[0]+2*x[1] and gcd(out,6)==1,'literal input and strong reader count')
            literal+=1
    for q in (5,7,11,25,35,49,55,77):
        seed=(1,2)
        seen={seed}
        queue=deque([seed])
        while queue:
            x=queue.popleft()
            for a in (AA,BB,CC):
                t=tuple(v%q for v in mv(a,x))
                if t not in seen:
                    seen.add(t)
                    queue.append(t)
        need(seen=={x for x in product(range(q),repeat=2) if gcd(*x,q)==1}, 'primitive class support')
    print(f'{literal} literal finished path counts; rational shear identities; 8 complete primitive residue closures')


if __name__=='__main__':
    main()
