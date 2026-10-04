#!/usr/bin/env python3
"""Compile rational growth anchors into guarded first-descent families.

Complete small-word universe: primitive growth necklaces with at most six
odd letters; all m=m0..m0+19 and four coefficient lifts. The specific
(-19/11) family is compared to the frozen 171-row / 65-cylinder old bank.
"""
from fractions import Fraction
from pathlib import Path
import re
import json


def check(test,message='check failed'):
    if not test: raise AssertionError(message)


def compose(word):
    A=p=B=0
    for a in word:
        B=3*B+2**A; A+=a; p+=1
    return 3**p,2**A,B


def step(n):
    z=3*n+1
    check(z>0)
    a=(z&-z).bit_length()-1
    return z>>a,a


def compositions(total,length):
    if length==1:
        yield (total,)
    else:
        for a in range(1,total-length+2):
            for tail in compositions(total-a,length-1):yield (a,)+tail


def growth_prefixes(word):
    A=0
    for i,a in enumerate(word,1):
        A+=a
        if 3**i<=2**A:return False
    return True


def growth_necklaces(bound):
    for p in range(1,bound+1):
        A=p
        while 2**A<3**p:
            for w in compositions(A,p):
                if any(p%j==0 and w==w[:j]*(p//j) for j in range(1,p)):continue
                rotations=[w[j:]+w[:j] for j in range(p)]
                if w!=min(rotations):continue
                yield min(v for v in rotations if growth_prefixes(v))
            A+=1


def anchor(word):
    P,Q,B=compose(word);check(P>Q and growth_prefixes(word))
    r=Fraction(B,P-Q)
    return P,Q,r.numerator,r.denominator


def cylinder(word,m):
    P,Q,h,d=anchor(word)
    check(Q**m>h)
    t=1
    while 2**t*(Q**m-h)<=P**m-h:t+=1
    beta=h*pow(P**m,-1,2**t)%2**t
    K=m*sum(word)+t
    residue=(beta*Q**m-h)*pow(d,-1,2**K)%2**K
    check(residue>0 and residue%2 and beta%2)
    return residue,K,t


def verify(word,m,n):
    P,Q,h,d=anchor(word);start=n
    check((d*n+h)%Q**m==0)
    b=(d*n+h)//Q**m;check(b>0 and b%2==1)
    actual=[]
    for i in range(m*len(word)):
        n,a=step(n);actual.append(a)
        check((n<start)==(i==m*len(word)-1),'not exact first descent')
    z=b*P**m-h; extra=(z&-z).bit_length()-1
    expected=word*(m-1)+word[:-1]+(word[-1]+extra,)
    check(tuple(actual)==expected,'wrong word')
    check(z%(d*2**extra)==0 and n==z//(d*2**extra),'wrong endpoint')
    MP,MQ,MB=compose(actual)
    check(MP*start+MB==MQ*n,'independent affine replay')
    return n


def intersects(a,b):
    return (a[0]-b[0])%2**min(a[1],b[1])==0


def log2_mod3(target,power):
    """Lift one discrete log by three candidates at each additional digit."""
    target%=3**power;check(target%3!=0)
    t=0 if target%3==1 else 1
    for k in range(1,power):
        candidates=[t+j*2*3**(k-1) for j in range(3)
                    if pow(2,t+j*2*3**(k-1),3**(k+1))==target%3**(k+1)]
        check(len(candidates)==1)
        t=candidates[0]
    return t


def sealed_certificate(m,target):
    """Power-expression certificate; its source need not be expanded."""
    cycles={1:[2],-1:[1],-5:[1,2],-17:[1,1,1,2,1,1,4]}
    check(target in cycles,'this seal requires a named terminal cycle')
    modulus=27**m
    t=log2_mod3(-19*pow(11*target,-1,modulus),3*m)
    period=2*3**(3*m-1)
    if t==0:t=period
    check((19+11*target*pow(2,t,modulus))%modulus==0)
    return {'chart':{'word':[1,1,2],'anchor_numerator':-19,'anchor_denominator':11},
            'repetitions':m,'terminal':target,'terminal_exponent':t,
            'terminal_cycle_word':cycles[target],
            'exponent_period':period,
            'source_expression':'(16^m*(19+11*r*2^t)-19*27^m)/(11*27^m)',
            'parameters':{'m':m,'r':target,'t':t},
            'actual_word':'(1,1,2)^(m-1) followed by (1,1,2+t)'}


def replay_sealed(cert):
    m,r,t=(cert['parameters'][s] for s in ('m','r','t'))
    numerator=16**m*(19+11*r*2**t)-19*27**m
    denominator=11*27**m
    check(numerator%denominator==0)
    source=numerator//denominator;current=source
    word=(1,1,2)*(m-1)+(1,1,2+t)
    for a in word:
        z=3*current+1
        check(z!=0 and (abs(z)&-abs(z)).bit_length()-1==a)
        current=z//2**a
    check(current==r and (source>0)==(r>0))
    P,Q,B=compose(word)
    check(P*source+B==Q*r)
    current=r
    for a in cert['terminal_cycle_word']:
        z=3*current+1
        check(z!=0 and (abs(z)&-abs(z)).bit_length()-1==a)
        current=z//2**a
    check(current==r,'terminal cycle does not close')
    cert['replayed_source_bit_length']=abs(source).bit_length()
    if abs(source).bit_length()<150:cert['replayed_source']=source
    return source


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificates',help='optional JSON destination for 12 complete signed certificates')
    args=parser.parse_args()
    words=list(growth_necklaces(6));cases=0
    for word in words:
        P,Q,h,d=anchor(word);m0=1
        while Q**m0<=h:m0+=1
        for m in range(m0,m0+20):
            r,K,t=cylinder(word,m)
            for lift in (0,1,7,29):
                verify(word,m,r+lift*2**K);cases+=1
    print('General compiler:',len(words),'primitive growth necklaces; exact first-descent controls',cases)

    word=(1,1,2);check(anchor(word)==(27,16,19,11))
    repo=Path(__file__).resolve().parents[2]
    text=(repo/'05-knowledge/results/reset_20260926_swaplift.out').read_text()
    rows=[]
    for line in text.splitlines():
        if re.fullmatch(r'\d+(?: \d+){8}',line):
            row=list(map(int,line.split()));rows.append((row[-2],row[-1]))
    check(len(rows)==171,'frozen bank universe')
    minimal=[]
    for r,k in sorted(set(rows),key=lambda p:(p[1],p[0])):
        if not any(k>=j and r%2**j==s for s,j in minimal):minimal.append((r,k))
    check(len(minimal)==65)
    check(sum((Fraction(1,2**k) for r,k in minimal),Fraction())==Fraction(6985206796614369409,2**65))
    cone=(743,11)
    check(cone[0]==-19*pow(11,-1,2**11)%2**11)
    check(all(not intersects(cone,row) for row in rows))
    print('All 171 old rows / 65 minimal cylinders avoid n=743 mod2048')
    mass=Fraction()
    for m in range(2,31):
        r,K,t=cylinder(word,m); endpoint=verify(word,m,r)
        if m>=3:
            check(r%2**11==743 and all(not intersects((r,K),old) for old in rows))
            mass+=Fraction(1,2**K)
        check(r%16==7,'must be disjoint from -5 sources 3 mod8 and -17 sources 15 mod16')
        if m<=10:print('rational return m',m,'t',t,'source',r,'modulus 2^',K,'endpoint',endpoint)
    check(cylinder(word,3)==(21223,15,3))
    tail=Fraction(1,26*27**30)
    print('Added density lower (m=3..30):',mass)
    print('Added density tail upper:',tail)
    print('Added density decimal display:',format(float(mass),'.17g'))
    # A deliberately weakened exit guard preserves the growing prefix but
    # leaves the selected endpoint above the original source.
    n=487;path=[n]
    for _ in range(6):n,a=step(n);path.append(n)
    check(n==695 and all(x>487 for x in path[1:]))
    print('One-bit-weakened exit hostile:',path)

    # Completed certificates: a source-decodable repeat chart, followed by
    # a prescribed final power-of-two division directly to the odd root 1.
    for m in range(1,7):
        power=3*m;modulus=3**power
        t=log2_mod3(-19*pow(11,-1,modulus),power)
        check(t>=1 and (19+11*pow(2,t,modulus))%modulus==0)
        print('Sealed source m',m,'terminal exponent t',t,'mod',2*3**(power-1))
        if m<=4:
            b=(19+11*2**t)//27**m
            source=(b*16**m-19)//11
            check(11*source+19==b*16**m)
            check(verify(word,m,source)==1)
            if m==1:check(source==151)
            print('  finite replay: source bit length',source.bit_length(),'odd steps',3*m,'target 1')
    certificates=[]
    for r in (1,-1,-5,-17):
        for m in range(1,4):
            cert=sealed_certificate(m,r);replay_sealed(cert);certificates.append(cert)
            print('Signed seal:',r,'m',m,'t',cert['terminal_exponent'],
                  'source bits',cert['replayed_source_bit_length'])
    if args.certificates:
        Path(args.certificates).write_text(json.dumps(certificates,indent=2)+'\n')
    print('ALL CHECKS PASSED')


if __name__=='__main__':main()
