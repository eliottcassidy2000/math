# Independent review of the bounded-cycle precision theorem

2026-10-25. **REVIEW PASS:** all stated mathematical claims in `README.md`
hold with their given scopes. **FINITE-EXACT:** the producer and independent
controls below pass normally and with `python -O`. No producer change was
required. This remains a precision filter of the odd-cycle Burnside sum.

## The valuation bound, including sharpness

For each partition, the denominator's binary valuation is precisely
sum_d v2(m_d!), since the cycle lengths d are odd. The numerator exponent
is (n-m)/2+G. Both comparisons in(P3) have the correct direction:
G>=binom(m,2), and multinomial integrality gives
sum_d v2(m_d!)<=v2(m!). The claim concerns rational summand valuations;
it does not assume individual summands are ordinary integers.

For odd n and odd m, subtracting B=(n-1)/2 gives exactly
(m-1)^2/2-v2(m!). The asserted sharp partition is valid for every odd
n>=m: its gcd sum is binom(m,2). If its large part is greater than1,
the only nontrivial multiplicity is m-1 and
v2((m-1)!)=v2(m!), because m is odd. If n=m, all parts are1 and the
denominator is m!, again giving equality. This verifies the quantifier
"every admissible odd n>=m", including m=1 and the all-ones boundary.

The relative bound equals0 at m=1,1 at m=3, and is positive for every
odd m>=5. Thus the n-cycle is the unique minimal-valuation summand;
its coefficient1/n is a binary unit. This reproves the odd valuation
without needing any cancellation assertion.

For even n>=4, m=2 must be handled separately. The coarse formula gives
-1 there, so using it as the baseline would be invalid. Distinct odd
parts have relative valuation gcd(a,b)-1>=0. Equal parts d,d have
valuation n-2, hence relative valuation d-2>=1, since their admissible
d is at least3. The separate n=2 exception is necessary and correctly
handled. For every even m>=4, the displayed coarse bound is at least1.

## Leading bit and equality classification

The only even-order terms of relative valuation0 are distinct coprime
odd pairs summing to n. Their odd units are each1 modulo2. Pairing the
units a modulo n with n-a gives exactly phi(n)/2 unordered pairs; there
is no fixed unit for n>=4. This proves(P6).

The equality v2(a(n))=n/2 is therefore equivalent to v2(phi(n))=1.
For n=2^s product p_i^a_i this valuation is
max(s-1,0)+sum v2(p_i-1). It equals1 precisely for n=4, or s=1 with
one odd prime divisor p congruent to3 modulo4. Arbitrary positive powers
of that prime are permitted. Thus the iff in(P7), including its prime-
power quantifier and all excluded even cases, is correct.

## Both normalized modulo-four formulas

For odd n, only m=1 and3 can contribute below relative valuation2.
The singleton gives1/n=n modulo4. For a distinct triple G>=3 and its
factorial valuation is0, giving relative valuation at least2. For a
repeated triple a,a,b the relative valuation is
a+2*gcd(a,b)-2; it equals1 only when a=1. For three equal parts it is
3a-2 and equals1 only when a=1. Therefore exactly the partition
(1,1,n-2), with its correct all-ones multiplicity at n=3, contributes2
modulo4. This proves(P8).

For even n, a distinct noncoprime pair has relative valuation at least2
because its gcd is odd and at least3. The distinct coprime pair products
are -1 modulo4 when n=0 modulo4, and +1 when n=2 modulo4; each inverse
has the same residue. Equal parts contribute at relative valuation1
only at (3,3), n=6.

For four cycles the relative valuation is G-2-sum v2(m_d!). The five
multiplicity patterns give exactly the claimed lower bounds:3 for3+1,
3 for2+1+1,4 for four distinct parts, at least4 for2+2, and6d-5 for
four equal d. In the2+2 case the parts a,b are distinct odd integers,
so a+b>=4 and G=a+b+4*gcd(a,b)>=8. The sole four-cycle relative-one
case is four ones at n=4. For even m>=6, beta_even(m)>=8: check m=6
directly, and use v2(m!)<=m-1 for m>=8. These facts prove exactly the
two exceptional corrections in(P9); neither is optional.

## Cycle budget and modular evaluation

The beta inequalities imply the stated O(sqrt(K)) cap. More explicitly,
beta_even(m)<K entails (m-2)^2<2K+2; the odd bound is slightly stronger.
Consequently2+floor(sqrt(2K+2)) is a safe upper cap even when the square
root is integral. The exact beta test removes any rounded overestimate.
The isolated even m=2 case is included separately and is always safe.
The table for K=1,...,6 is correct.

Each discarded summand individually has valuation at least B+K. Its
removal is therefore valid independently of signs or cancellation.
For a retained term, factoring z=2^t*z_odd before reduction gives
2^(e-t-B)/z_odd. The exponent is nonnegative by the proved baseline,
and z_odd is invertible modulo2^K. This establishes(P11). Inverting the
original even z, or reducing its numerator before cancelling2^t, would
not implement this calculation. The fixed-length recursion generates
all permitted partitions exactly once and never produces an excluded
cycle count. Low precision bounds the number of parts, not their sizes
or the complexity of recovering the complete integer.

## Reproduction and independent controls

Producer `check.py` passes normally and under `-O`:31,583 full partition
records,306 modular comparisons, n=0..50 and K=1..6. The following
additional control uses the separately proved integer transfer in
`counting/check.py` for complete A000568 counts, rather than the
precision producer's rational partition oracle. It checks both low-bit
formulas and modular residues through n=60, sharpness on every odd pair
1<=m<=n<=101, all triple/four-cycle low-valuation exceptions in the stated
finite ranges, totient pair counts through500, and budget rounding.
No tournament or Hamiltonian-path enumeration is used.

Run this embedded block from the repository root using:

    @'
    from pathlib import Path
    p=Path('05-knowledge/results/port_residual_integral_20261025/precision/INDEPENDENT_REVIEW.md')
    code=p.read_text(encoding='utf-8').split('```python\n',1)[1].split('\n```',1)[0]
    exec(compile(code,str(p),'exec'))
    '@ | python -

Use `python -O -` instead of `python -` for the optimized replay.

```python
from collections import Counter
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt
from pathlib import Path
import json

def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)

def load(name,path):
    spec=spec_from_file_location(name,path)
    module=module_from_spec(spec); spec.loader.exec_module(module)
    return module

root=Path('05-knowledge/results/port_residual_integral_20261025')
integers=load('independent_integral_counts',root/'counting/check.py')
precision=load('precision_under_review',root/'precision/check.py')

def fv(m):
    answer=0
    while m:
        m//=2; answer+=m
    return answer

def relative(parts):
    n,m=sum(parts),len(parts)
    return ((n-m)//2+sum(gcd(a,b) for a,b in combinations(parts,2))
            -sum(fv(k) for k in Counter(parts).values())-n//2)

def factors(n):
    out={}; p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1:
        out[n]=out.get(n,0)+1
    return out

def phi(n):
    answer=n
    for p in factors(n):
        answer=answer//p*(p-1)
    return answer

sharp=0
for n in range(1,102,2):
    for m in range(1,n+1,2):
        part=(1,)*(m-1)+(n-m+1,)
        need(relative(part)==(m-1)**2//2-fv(m),('sharp odd bound',n,m))
        sharp+=1

triple_cases=four_cases=0
for n in range(3,62,2):
    for part in combinations_with_replacement(range(1,n+1,2),3):
        if sum(part)==n:
            need((relative(part)==1)==(part==(1,1,n-2)),('triple equality',part))
            triple_cases+=1
for n in range(4,43,2):
    for part in combinations_with_replacement(range(1,n+1,2),4):
        if sum(part)==n:
            need((relative(part)==1)==(part==(1,1,1,1)),('four-cycle equality',part))
            shape=sorted(Counter(part).values())
            lower={(1,1,1,1):4,(1,1,2):3,(1,3):3,(2,2):4}.get(tuple(shape))
            if lower is not None:
                need(relative(part)>=lower,('four-cycle lower bound',part))
            else:
                need(relative(part)==6*part[0]-5,('four equal parts',part))
            four_cases+=1

for n in range(4,501,2):
    pairs=[(a,n-a) for a in range(1,n//2+1,2) if a<n-a and gcd(a,n-a)==1]
    need(len(pairs)==phi(n)//2,('totient pair count',n))
    need(all((a*b)%4==(3 if n%4==0 else 1) for a,b in pairs),('pair product residues',n))

for K in range(1,129):
    cap=2+isqrt(2*K+2)
    for parity,n in ((1,1001),(0,1000)):
        expected=[]
        for m in range(1 if parity else 2,n+1,2):
            beta=((m-1)**2//2 if parity else m*(m-2)//2)-fv(m)
            if m==2:
                beta=0
            if beta<K:
                need(m<=cap,('safe budget cap',K,m))
                expected.append(m)
        need(precision.cycle_counts(n,K)==expected,('exact budget selection',K,parity))

modular=0
for n in range(61):
    value,_=integers.integral_target(n,core=False)
    B=n//2 if n>=3 else 0
    need(value%(1<<B)==0,('baseline from integer engine',n))
    quotient=value>>B
    if n>=3 and n%2:
        need(quotient%4==(n+2)%4,('odd modulo four',n))
    elif n>=4:
        fact=factors(n); odd={p:e for p,e in fact.items() if p!=2}
        classified=(n==4 or (fact[2]==1 and len(odd)==1 and next(iter(odd))%4==3))
        need((quotient%2==1)==classified,('even iff class',n))
        expected=(-1)**(n//2+1)*(phi(n)//2)+2*int(n in (4,6))
        need(quotient%4==expected%4,('even modulo four',n))
    for K in (1,2,6):
        got,_=precision.normalized_residue(n,K)
        need(got==quotient% (1<<K),('independent integer/modular agreement',n,K))
        modular+=1
need(relative((1,1))==-1,'order two really violates even baseline')
need(relative((3,3))==1,'equal-three correction is real')
value50,_=integers.integral_target(50,core=False)
need((value50 & -value50).bit_length()-1==26,'order50 valuation consequence')
print(json.dumps({'status':'PASS','sharp_odd_cases':sharp,
                  'triple_cases':triple_cases,'four_cycle_cases':four_cases,
                  'totient_pair_orders':[4,500],'budget_precisions':[1,128],
                  'independent_integer_orders':[0,60],
                  'independent_modular_comparisons':modular},indent=2))
```
