#!/usr/bin/env python3
"""Exact first-departure phase for the periodic shortcut-27 comparison.

Standard library only. Generic parity isometry is inherited; this checks
its explicit 59-bit specialization and both native clocks.
"""
from fractions import Fraction
import json

WORD = tuple(map(int, '11011111010110111011110100111011011111100111100010101000100'))
Q, A = 2**59, 3**37
B, D = 23*Q-27*A, Q-A
RHO = Fraction(B,D)
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def v2(n):
    if not n:
        raise ValueError('zero has no finite valuation')
    n=abs(n)
    return (n & -n).bit_length()-1


def shortcut(n):
    return (3*n+1)//2 if n%2 else n//2


def rat_step(n):
    need(n.denominator%2==1,'odd-denominator rational domain')
    return (3*n+1)/2 if n.numerator%2 else n/2


def ordinary(n):
    return 3*n+1 if n%2 else n//2


need(len(WORD)==59 and sum(WORD)==37,'fixed word clocks')
need(D%2==1 and B%2==1 and 8<RHO<9,'nonintegral template')
r=RHO
for bit in WORD:
    need(r.numerator%2==bit,'rational template bit')
    r=rat_step(r)
need(r==RHO,'full rational period')


def inspect(n):
    need(isinstance(n,int) and n>0,'positive source')
    ell=v2(D*n-B)
    cur=n
    odd=0
    for i in range(ell):
        need(cur%2==WORD[i%59],'common native bit')
        odd+=cur%2
        cur=shortcut(cur)
    need(cur%2==1-WORD[ell%59],'first native discrepancy')
    k,r=divmod(ell,59)
    clock=ell+37*k+sum(WORD[:r])
    need(clock==ell+odd,'ordinary clock combinatorics')
    direct=n
    for _ in range(clock):
        direct=ordinary(direct)
    need(direct==cur,'independent ordinary replay')
    need(shortcut(cur)==((3*cur+1)//2 if 1-WORD[r] else cur//2),'first departure branch')
    return {'source':str(n),'ell':ell,'blocks':k,'phase':r,'ordinary_clock':clock,'departure':str(cur),'actual_bit':cur%2,'template_bit':WORD[r]}


for n in range(1,4097):
    inspect(n)

# Every residual phase at five block depths, three positive lifts per exact
# valuation class. Residue construction is independent of trajectory replay.
phase_sources=0
for k in range(5):
    for phase in range(59):
        ell=59*k+phase
        modulus=2**(ell+1)
        res=((B+2**ell)*pow(D,-1,modulus))%modulus
        for lift in (1,2,7):
            n=res+lift*modulus
            need(v2(D*n-B)==ell,'exact cylinder source')
            row=inspect(n)
            need(row['blocks']==k and row['phase']==phase,'phase and full-block projection')
            phase_sources+=1

# Every low cylinder count checked directly against the arithmetic progression
# formula; this tests the finite boundary behind the stated density.
for j in range(13):
    modulus=2**(j+1)
    res=((B+2**j)*pow(D,-1,modulus))%modulus
    for X in (0,1,27,100,4096):
        direct=sum(v2(D*n-B)==j for n in range(1,X+1))
        formula=(X-res)//modulus-(-res)//modulus
        need(direct==formula,'exact finite cylinder count')

special=[inspect(n) for n in (27,27+Q,27+2*Q,209923714958077886127437625648218139)]
need(special[0]['ell']==61 and special[0]['departure']=='53' and special[0]['ordinary_clock']==100,'27 sharp phase')
need(special[0]['blocks']==special[1]['blocks']==1 and special[0]['actual_bit']!=special[1]['actual_bit'],'whole-block count loses next branch')

print(json.dumps({
 'status':'PASS','checks':CHECKS,
 'scope':'Inherited parity-isometry mechanism, explicit all-source phase proof, finite replay controls only.',
 'universe':{'positive_sources':[1,4096],'block_depths':[0,4],'phases':[0,58],'lifts':[1,2,7],'phase_sources':phase_sources,'count_depths':[0,12],'count_caps':[0,1,27,100,4096]},
 'word':''.join(map(str,WORD)),
 'special_controls':special,
 'limitation':'No global stopping-time, recurrence-entry or convergence theorem.'
},indent=2,sort_keys=True))
