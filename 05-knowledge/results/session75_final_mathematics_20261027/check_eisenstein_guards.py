#!/usr/bin/env python3
"""Small exact tests of the proved potentials, pairing identities and guards.
This does not replay the451-state catalyst certificate.
"""
from itertools import product


def need(ok,message):
    if not ok: raise ValueError(message)


def Q(z):
    a,b=z; return a*a-a*b+b*b


def add(x,y): return (x[0]+y[0],x[1]+y[1])


def potential():
    points=list(product(range(-8,9),repeat=2))
    for x in points:
        for y in points:
            z=add(x,y); beta=Q(z)-Q(x)-Q(y)
            need(max(beta,0)+Q(x)//3+Q(y)//3-Q(z)//3>=0,'base potential')
    shell=[z for z in product(range(-22,23),repeat=2) if Q(z)==121]
    need(len(shell)==6,'norm121 shell')
    need({z for z in shell if z[0]%2==1 and z[1]%2==0}=={(11,0),(-11,0)},'parity shell')
    need(not any(Q(z)==123 for z in product(range(-23,24),repeat=2)),'123 nonnorm')
    need(pow(11,820,41*41)==41*41-1,'infinite-family congruence')
    x,y=(-1,-4),(12,4); beta=Q(add(x,y))-Q(x)-Q(y)
    slack=max(beta,0)+Q(x)//3+Q(y)//3-Q(add(x,y))//3
    need(slack==1,'two-unit bump hostile control')
    need(sum(Q(z)<=124 for z in product(range(-23,24),repeat=2))==451,'finite disk size')
    print('PASS:83521 base-slack pairs; norm/parity shell;451-state disk count; hostile two-unit bump.')


def pairings():
    pairs=[((0,2),(1,3)),((0,1),(2,3)),((0,3),(1,2))]
    for x in product(range(-2,3),repeat=4):
        parent=sum(x)**2; children=sum(a*a for a in x)
        middle=[sum((x[i]+x[j])**2 for i,j in pairing) for pairing in pairs]
        need(sum(middle)==parent+2*children,'pairing average')
        a,b,c,d,e,f=[x[i]+x[j] for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))]
        need(a+f==b+e==c+d and (a+b-d)%2==0,'integral image')
        u=(a+b-d)//2
        need((u,a-u,b-u,c-u)==x,'reconstruction')
    x=(-2,-1,0,1); bills=[]
    for pairing in pairs:
        children=sum(a*a for a in x)
        mid=sum((x[i]+x[j])**2 for i,j in pairing)
        increments=[(x[i]+x[j])**2-x[i]**2-x[j]**2 for i,j in pairing]+[sum(x)**2-mid]
        need(sum(increments)==sum(x)**2-children==-2,'signed invoice')
        bills.append(sum(max(z,0) for z in increments))
    need(bills==[0,4,2],'positive ordering control')
    seeds=(33,41,49,57); targets=[]; clocks=[]
    for s in seeds:
        a=2 if s%8==1 else 1; p=(3*s+1)//2**a; j=0
        while p%8==5:
            p=(p-1)//4; j+=1
        t=a+2*j
        need(2**t*(3*p+1)==3*(3*s+1)+2**a,'native carry')
        targets.append(p); clocks.append(t)
    need(targets==[25,31,9,43] and clocks==[2,2,4,2],'native seed examples')
    print('PASS:625 pairing/reconstruction examples; invoices0,4,2; four native carries.')


if __name__=='__main__':
    potential(); pairings()
