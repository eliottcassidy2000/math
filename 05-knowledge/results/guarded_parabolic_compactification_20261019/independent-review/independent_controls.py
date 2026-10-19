#!/usr/bin/env python3
"""Independent exact controls. No imports from the reviewed implementation."""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd
import json

count = 0

def check(value, label):
    global count
    count += 1
    if not value:
        raise RuntimeError(label)

def step(x):
    return (3*x + 1)//2 if x % 2 else x//2

def orbit(x, length):
    bits = []
    for _ in range(length):
        bits.append(str(x % 2))
        x = step(x)
    return ''.join(bits), x

def affine(word):
    # Closed sum over positions, rather than forward affine accumulation.
    positions = [i for i, z in enumerate(word) if z == '1']
    h = len(positions)
    return (3**h, sum(3**(h-k-1) * 2**i
                     for k, i in enumerate(positions)), 2**len(word))

def cylinder(word):
    # Lift the canonical input one binary digit at a time using actual parity.
    # No modular inverses or repeated-word formula enter this construction.
    r = y = 0
    modulus = slope_numerator = 1
    for bit in map(int, word):
        if y % 2 != bit:
            r += modulus
            y += slope_numerator
        check(y % 2 == bit, 'binary lift parity')
        if bit:
            y = (3*y + 1)//2
            slope_numerator *= 3
        else:
            y //= 2
        modulus *= 2
    return r, y

digest = sha256()
words = squares = lifted_repetitions = 0
for n in range(1, 13):
    for bits in product('01', repeat=n):
        word = ''.join(bits)
        a, b, M = affine(word)
        r, y = cylinder(word)
        check(a*r + b == M*y, 'closed carry agrees with binary lift')
        check(orbit(r, n) == (word, y), 'actual canonical trajectory')
        check(4*b <= M*(a-1), 'universal carry estimate')
        words += 1
        for k in (1, 2, 3):
            rr, yy = cylinder(word*k)
            modulus = M**k
            check(((a-M)*rr+b) % modulus == 0, 'repetition congruence')
            check(gcd(a-M, modulus) == 1, 'unique repeated congruence')
            aa, bb, MM = affine(word*k)
            check(aa*rr+bb == MM*yy, 'repeated closed affine formula')
            check((a-M)*bb == b*(a**k-M**k), 'repeated carry identity')
            for lift in (0, 1, 4):
                x = rr + lift*modulus
                if x == 0:
                    continue
                replayed, endpoint = orbit(x, n*k)
                check(replayed == word*k, 'positive lifted repetition')
                one_word, one_endpoint = orbit(x, n)
                check((one_endpoint-x) % M**(k-1) == 0,
                      'displacement lattice')
                mismatch = (a-M)*x+b
                if mismatch:
                    v = (abs(mismatch) & -abs(mismatch)).bit_length()-1
                    budget = v//n
                    check(orbit(x, budget*n)[0] == word*budget,
                          'all repetitions through exact budget are lawful')
                    check(orbit(x, (budget+1)*n)[0] != word*(budget+1),
                          'valuation is exact repetition budget')
                if a<M and k>=2:
                    check(one_endpoint <= x, 'contracting repeated descent')
                if a>M and k>=2:
                    check(one_endpoint-x >= M**(k-1), 'expansion gap')
                digest.update(f'{word},{k},{x},{one_endpoint},{endpoint}\n'.encode())
                lifted_repetitions += 1
            if a<M and k==2:
                # This affine inequality proves all positive lifts, not samples.
                minimal_positive = rr if rr else modulus
                check((M-a)*minimal_positive >= b, 'all-lift square descent')
                squares += 1

K = '110110110110110101011011101010'
F = '10'
H = '110101000'
V = K*4+F
Kplus = '001111111111111101010110000100'
Fplus = '01'
Hplus = '001101001'
Vplus = Kplus*4+Fplus
P, C0, Q = affine(K)
D = P-Q
c = Fraction(C0, D)
check((P,C0,Q,D) == (3**19,4912087879,2**30,88519643), 'bank constants')
va,vb,vM=affine(V)
check(Fraction(vb,va-vM)+1 == Fraction(
    83986881009841395726881741282576539890,
    157489106280555890461927692448373459), 'inherited episode carry constant')

for lower, upper in ((K,Kplus),(F,Fplus),(H,Hplus),(V,Vplus)):
    a,b,M = affine(lower)
    aa,bb,MM = affine(upper)
    r,y = cylinder(lower)
    rr,yy = cylinder(upper)
    check((aa,MM,bb)==(a,M,b+M-a), 'companion affine identity')
    check(rr==(r+1)%M, 'separate companion cylinder')
    check(orbit(r+1,len(upper))==(upper,y+1), 'actual adjacent companion')
for word in (K,V):
    for j in range(1,len(word)+1):
        a,b,M=affine(word[:j])
        check(a>M, 'every expansion prefix expands')
for word in (F,H):
    a,b,M=affine(word)
    for j in range(len(word)):
        aa,bb,MM=affine(word[:j])
        check(Fraction(aa,MM)>Fraction(a,M), 'reset proper-prefix slope')
        check(Fraction(aa+bb,MM)>=Fraction(a+b,M), 'reset comparison at one')

short_words = [K*l+reset for l in (1,2,3) for reset in (F,H)]
short_words += [K*4+H]
short_words += [V+K*i+reset for i in (0,1,2) for reset in (F,H)]
for word in short_words:
    a,b,M=affine(word)
    r,y=cylinder(word)
    check(a<M and r>0 and (M-a)*r>b, 'thirteen all-lift source descents')

# Enumerate the entire phase orbit, rather than only a group-order certificate.
inverse_Q = pow(Q,-1,D)
initial = C0*inverse_Q % D
phase = initial
seen = set()
phase_digest=sha256()
while phase not in seen:
    seen.add(phase)
    phase_digest.update(f'{phase}\n'.encode())
    phase=phase*inverse_Q%D
check(phase==initial and len(seen)==146768, 'complete modular phase cycle')
check(gcd(C0,D)==1 and D==193*458651, 'unit phase and denominator')
check(all(1<=t<D for t in seen), 'all phase limits distinct and nonzero')
for L in range(1,21):
    r,y=cylinder(K*L)
    modulus=Q**L
    t=(D*r+C0)//modulus
    check(D*r+C0==t*modulus, 'K exact compact scaling')
    check(Fraction(r,modulus)==Fraction(t,D)-c/modulus, 'real mantissa identity')
    check(D*r>=modulus-C0, 'uniform minimum K source height')

near_neutral=[]
for L,j in ((11,3),(40,11),(69,19)):
    word=K*L+F*j
    a,b,M=affine(word)
    s,e=Fraction(a,M),Fraction(b,M)
    check(s>1 and e==c*s+1-(c+1)*Fraction(3,4)**j, 'near-neutral branch')
    r,y=cylinder(word*2)
    _,one=orbit(r,len(word))
    check(one-r>=M, 'near-neutral repeated source displacement gap')
    check(Fraction(b,(a-M)*r)<=Fraction(b,M*M-b), 'nonlinear crossover exclusion')
    near_neutral.append({'K_blocks':L,'F_blocks':j,'length':len(word)})

offsets=[4,12852,14196,15364,18436,49324,128308,215276,266092,288932,
         295660,310916,390388,394148,480032,522404,529652]
anchors=[]
for d in offsets:
    z=Fraction(d)-c
    bits=[]
    for j in range(1,10001):
        bit=z.numerator%2
        bits.append(str(bit))
        z=(3*z+1)/2 if bit else z/2
        a=3**bits.count('1')
        M=2**j
        if a<M:
            break
    else:
        raise RuntimeError('anchor control exhausted declared bound')
    word=''.join(bits)
    a,b,M=affine(word)
    check(z==(a*(Fraction(d)-c)+b)/M, 'rational orbit and closed branch')
    threshold=Fraction(a*d+b,M-a)
    L=(j+29)//30
    while Q**L<=D*threshold+C0:
        L+=1
    r,y=cylinder(K*L)
    check(r>threshold, 'all K lifts exceed original-height threshold')
    for lift in (0,2,11):
        x=r+lift*Q**L
        ww,yy=orbit(x+d,j)
        check(ww==word and yy<x,'exact shifted rational-anchor repair')
    anchors.append({'d':d,'coefficient_stop':j,'sufficient_K_blocks':L})
check(anchors[:2]==[{'d':4,'coefficient_stop':92,'sufficient_K_blocks':4},
                   {'d':12852,'coefficient_stop':135,'sufficient_K_blocks':5}],
      'declared short anchor controls')

# Check actual endpoint overshoots, then all candidate word/offset pairs below.
endpoint_pairs=endpoint_realizations=0
for n in range(1,11):
    for x in range(1,513):
        word,y=orbit(x,n)
        a,b,M=affine(word)
        if a<M and y>x:
            d=y-x
            check((b-M*d)==(M-a)*x, 'actual overshoot solves endpoint equation')
            check(cylinder(word)[0]==x%M, 'endpoint candidate authenticates parity')
            endpoint_realizations+=1
        endpoint_pairs+=1
for n in range(1,10):
    for bits in product('01',repeat=n):
        word=''.join(bits)
        a,b,M=affine(word)
        if a>=M:
            continue
        for d in range(1,17):
            num=b-M*d
            den=M-a
            if num>0 and num%den==0:
                x=num//den
                check(orbit(x,n)==(word,x+d), 'positive integer filter sufficient')

# Boundary obstructions, and the shorter self-matching specialization of A2.
for d in (1,4,12852,529652):
    for J in range(1,97):
        x=2**J*(d+2)-d-1
        word,y=orbit(x+d,J)
        check(word=='1'*J and y>x,'all-odd boundary prevents uniform deadline')
for n in range(2,13):
    for x in range(1,1025):
        word,y=orbit(x,n)
        a,b,M=affine(word)
        if a<M and (y-x)%(2**(n-2))==0:
            check(y<=x,'n-minus-two matching-bit consequence')
check(orbit(7,8)==('11101001',8),'unguarded positive overshoot hostile')
check(orbit(1,4)==('1010',1) and orbit(2,4)==('0101',2),'return exceptions')

print(json.dumps({'status':'PASS','explicit_checks':count,
  'scope':'finite independent controls; written proofs carry all-size quantifiers',
  'all_words_through_length':12,'word_count':words,
  'contracting_square_all_lift_certificates':squares,
  'positive_lifted_repetitions':lifted_repetitions,
  'repetition_digest':digest.hexdigest(),
  'inherited_short_all_lift_exclusions':len(short_words),
  'enumerated_phase_period':len(seen),'phase_digest':phase_digest.hexdigest(),
  'near_neutral_examples':near_neutral,'rational_anchor_controls':anchors,
  'endpoint_source_length_pairs':endpoint_pairs,
  'positive_overshoot_endpoint_controls':endpoint_realizations},
  indent=2,sort_keys=True))
