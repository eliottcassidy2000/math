#!/usr/bin/env python3
"""Narrow independent exact audit of native 27-block exhaustion and refill."""
from fractions import Fraction
import json

checks = 0


def check(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)


def step(n):
    return (3*n+1)//2 if n % 2 else n//2


def ordinary(n):
    return 3*n+1 if n % 2 else n//2


def v2(n):
    if not n:
        raise RuntimeError('valuation of zero')
    n=abs(n)
    return (n & -n).bit_length()-1


# Derive the word and affine coefficients from actual native dynamics,
# rather than importing the producer's bit word or block formula.
n=27
word=[]
aa,bb,qq=1,0,1
for j in range(59):
    bit=n%2
    word.append(bit)
    n=step(n)
    aa,bb,qq=(3 if bit else 1)*aa, (3 if bit else 1)*bb+bit*qq, 2*qq
check(n==23 and sum(word)==37, 'derived block endpoint and odd count')
q,a,b=qq,aa,bb
d=q-a
check((q,a,b)==(2**59,3**37,23*2**59-27*3**37),'affine identity')
rho=Fraction(b,d)
check(8<rho<9 and d%2==1 and b%2==1,'comparison domain')
r=rho
for bit in word:
    check(r.numerator%2==bit,'actual rational word')
    r=(3*r+1)/2 if bit else r/2
    check(r.denominator%2==1,'odd denominator preserved')
check(r==rho,'rational block closes')

# First mismatch and ordinary clock, without a precomputed departure endpoint.
ell=v2(d*27-b)
n=27
clock=0
for j in range(ell):
    check(n%2==word[j%59],'27 common prefix')
    clock+=1+n%2
    n=step(n)
check((ell,n,clock)==(61,53,100),'27 exact departure')
check(n%2!=word[ell%59],'27 actual first mismatch')
nn=27
for _ in range(100):
    nn=ordinary(nn)
check(nn==53,'independent ordinary-clock replay')

c=23*d-2*b
check(c==31*a-23*q,'reset cofactor')
fixtures=[]
for h in (1,3):
    modulus=2**(59*h+2)
    th=((2**(59*h+1)-c)*pow(a*d,-1,modulus))%modulus
    check(th>0 and th%2==1,'positive odd reset residue')
    for lift in (0,1):
        t=th+modulus*lift
        source=27+q*t
        actual=source
        for j in range(59):
            check(actual%2==word[j],'reset source actual itinerary')
            actual=step(actual)
        middle=actual
        after=step(middle)
        check(middle==23+a*t and middle%2==0,'actual even reset endpoint')
        check(after==middle//2 and source>middle>after>0,'reset positivity/descent')
        check(tuple(v2(d*x-b) for x in (source,middle,after))==(59,0,59*h),
              'exact 1,0,h fuel')
        for j in range(59*h):
            check(after%2==word[j%59],'refilled native suffix')
            after=step(after)
        check(v2(d*after-b)==0,'suffix exhaustion endpoint')
        fixtures.append({'h':h,'lift':lift,'source_ell':59,'middle_ell':0,
                         'refilled_ell':59*h,'suffix_terminal_ell':0})
    # Same visible current state can have different next exact fuel.
    lift=(-th*pow(modulus,-1,7))%7
    t=th+modulus*lift
    middle=23+a*t
    check(t%7==0 and (middle%7,(middle//2)%7)==(2,1), 'same modulo-seven state')
    check((v2(d*middle-b),v2(d*(middle//2)-b))==(0,59*h),
          'same-state distinct refill')

print(json.dumps({'passed':True,'checks':checks,
  'status':'FINITE-EXACT supporting independent all-h proof audit',
  'universe':'derived 59-step word; exact rational period; n=27 departure; h=1,3 with two lifts each',
  'departure_27':{'ell':ell,'value':n,'ordinary_clock':clock},
  'fixtures':fixtures},indent=2,sort_keys=True))
