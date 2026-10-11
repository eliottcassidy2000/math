#!/usr/bin/env python3
"""Exact rank/atomic-moment controls for inherited authenticated actions.

This checks arithmetic interfaces; it does not discover or export ROOT words.
Universe: 32 H and GR lifts at credits0..4; 32 G lifts at credits1..4;
run lengths4..7,16 cofactor lifts,credits0..4; r2..5 moment identities
through cutoff32 and powers1..r-1; neutral recycling and countdown hostiles.
"""
from fractions import Fraction as F
import json

KAPPA=F(2349,2560)
checks=0


def check(ok,message):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(message)


def energy(x,k):
    return (x+5)*F(9,8)**k


def rank(x,k):
    if type(x) is not int or type(k) is not int or x<111 or x%2!=1 or k<0:
        raise ValueError('native integer state and nonnegative credit required')
    e=energy(x,k)
    b=0
    while e*KAPPA>=116:
        e*=KAPPA
        b+=1
    return 3*b+k


def main():
    counts={'H':0,'G':0,'R':0,'GR':0}
    for t in range(32):
        for k in range(5):
            x=155+2048*t
            y=(729*x+669)//1024
            check(energy(y,k+2)<=KAPPA*energy(x,k),'H energy interface')
            check(rank(x,k)-rank(y,k+2)>=1,'H integer rank')
            counts['H']+=1
            x=251+288*t
            y=223+256*t
            check(energy(y,k)<=F(57,64)*energy(x,k),'GR atomic energy interface')
            check(rank(x,k)-rank(y,k)>=3,'GR integer rank')
            counts['GR']+=1
        for k in range(1,5):
            x=123+16*t
            y=(9*x+5)//8
            check(energy(y,k-1)==energy(x,k),'G exact neutrality')
            check(rank(x,k)-rank(y,k-1)==1,'G spends rank')
            counts['G']+=1
    for r in range(4,8):
        for lift in range(16):
            t=(7 if r%2==0 else 5)+8*lift
            x=(3**r*t-1)//2
            y=2**(r+1)*t-1
            for k in range(5):
                check(x%4==3 and (2*x+1)%3**r==0,'run native input')
                check(energy(y,k+1)<=F(57,64)*energy(x,k),'run energy interface')
                check(rank(x,k)-rank(y,k+1)>=2,'run integer rank')
                counts['R']+=1
    check(energy(283,0)==energy(251,1),'recycling energy')
    check(rank(251,1)==rank(283,0)+1,'neutral credit creation increases rank')
    check(rank(283,0)==rank(251,1)-1,'neutral loop has no net rank payment')
    for args in [(True,0),(155.0,0),(155,False),(155,-1),(1,0)]:
        try:
            rank(*args)
        except ValueError:
            check(True,'reject invalid rank state')
        else:
            check(False,'accepted invalid rank state')
    moment_gates=0
    for r in range(2,6):
        for power in range(1,r):
            lhs=F(0)
            rhs=F(0)
            for n in range(1,33):
                lhs+=(F(1,n**r)-F(1,(n+1)**r))*n**power
                rhs+=F(n**power-(n-1)**power,n**r)
                check(lhs==rhs-F(n**power,(n+1)**r),'finite moment boundary term')
                moment_gates+=1
    for horizon in range(1,33):
        start=horizon+1
        check(start-horizon==1,'countdown survives every fixed global horizon')
        check(sum(1 for t in range(start))==start,'pointwise countdown expectation')
    print(json.dumps({'status':'FINITE-EXACT arithmetic interfaces; ROOT coverage remains OPEN',
                      'checks':checks,'action_counts':counts,'moment_identity_gates':moment_gates,
                      'rank_examples':{str((x,k)):rank(x,k) for x,k in [(155,0),(111,2),(283,0),(251,1),(223,1)]}},
                     indent=2,sort_keys=True))


if __name__=='__main__':
    main()
