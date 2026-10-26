# Independent review of the fixed-precision constant-tail quotient

2026-10-25. **ACCEPTED / PROVED:** the operation-preserving quotient in
[MODULAR.md](MODULAR.md), with the stated absolute-precision and
semantic-support boundaries. **FINITE-EXACT:** the producer and independent
checks below. No claim of minimum module rank, uniformly bounded metadata,
or improved ordinary scalar state count follows.

## Three different lattices remain distinct

The preceding integral theorem preserves the full product lattice in
the normalized Y coordinates. It does not impose congruences between
arbitrary coordinates in that lattice.

Actual moment banks satisfy

    Y_r=sum_(m>=b) x_m 2^(r(m-b)),   x_m in Z.

Modulo2^K, all m>b contributions vanish once r>=K, while the m=b
contribution survives. Thus Y_r=x_b for r>=K. This proves (M1),
including banks with negative integer weights or a zero residue at x_b.

The constant-tail module is a permissible larger domain, not a claim
to characterize actual moments. For K=2, the vector with Y_1=1 and
Y_r=0 for every r>=2 belongs to that module but cannot be an integer
moment bank: any such bank has Y_1=Y_2 modulo2. Conversely the full
product-lattice vector(Y_1,Y_2)=(0,1) modulo2 is outside the
constant-tail module and cannot be compressed this way. The independent
code tests its uncompressed identity transfer against the incorrectly
compressed transfer; they disagree exactly at channel2. The producer’s
revised hostile now tests this actual operation too.

## Transfer, mask deletion and terminal boundary

Write the inherited integral batch as

    Y'_r=T 2^(J0+r*delta) Y_(a+r),
    a=ell*k, J0>=0, delta=k+b-b'>=0.

The affordable source index satisfies a+r<=S. If a destination index
r>=K exists, then its source index also lies in the constant tail.
For delta=0 its whole contribution is independent of r; for delta>=1
it vanishes modulo2^K. Therefore every batch preserves the larger
constant-tail module, and addition preserves it. This argument covers
negative T, zero T, and all coordinate-mask drops. It does not need an
underlying sparse moment representation.

A mask drop can change the gauge. Three past9-cycles give
D(3)=0,D(5)=12, hence b=3 before the5 stage. Skipping length5 deletes
D(5) and gives b'=0. The shift has delta=3, not zero. A relaxed moment
weight x_3=1 therefore sends its nonzero all-one tail to zero modulo8.
The independent control includes this allowed defect history, rather
than assuming that dropping a coordinate is an unscaled identity.

S=0 has its separate channel Y_0. A positive-horizon source can reach
it only through an affordable positive index a=ell*k; k=0,e=0 cannot
exhaust S>0. A terminal source reads its actual channel0. No equation
identifies Y_0 with the positive-index tail. For instance the formal
integer moment bank x_1=1,x_2=-1 with b=1 has Y_0=0 and Y_1=1
modulo2. The implementation's source-index and terminal conventions
respect this distinction.

## Adaptive support and precision scope

The numerical closure proof applies to arbitrary constant-tail inputs.
The adaptive sparse-support proof uses additional exact counting semantics:
at fixed u,m all actual prefix terms have sign(-1)^(u-m), and a legal
batch shifts support injectively by m->m+k with a nonzero integer
multiplier. Those facts prevent integer cancellation within a semantic
cycle-count key. They would not hold for arbitrary externally supplied
signed histories.

The implementation retains zero-residue sparse keys and nonterminal
zero dense banks. Its dense conversion threshold h>min(S,K) is therefore
based on semantic support. Since S only decreases, such a dense state
never needs sparse reconstruction. It uses min(h,S,K) scalar weight
residues, with support labels and defect keys accounted for separately.
Ordinary enumeration already uses one scalar and gains no additional
state reduction here.

Every denominator in the modular coefficient and final O(N) removal is
odd. Reducing an exact numerator to K bits and multiplying by an odd-unit
inverse recovers the quotient modulo2^K without extra digits. This is a
statement about residue arithmetic; the order, defect coordinates,
support metadata and shift exponents are not themselves restricted to
K-bit integers.

Dividing instead by2^B requires B+K absolute digits to recover K quotient
digits by this route. The independent unmerged Burnside computation
below confirms the explicit hostile:

    a(50) modulo4 =0;
    a(50) modulo2^27 =67108864=2*2^25;
    a(50)/2^25 modulo4 =2.

This does not conflict with the separate normalized bounded-cycle
algorithm, which cancels the baseline before reducing.

## Reproduction

The producer passed normal and optimized Python with byte-identical
outputs:498 exact-count comparisons and1,489 arbitrary constant-tail
transfers plus1,489 sparse/dense squares. The order44 and order50
boundaries were included. The revised product-lattice and nonzero-tail
controls were read in the producer source.

The independent code below imports neither producer nor predecessor
modules. It checks74,592 images of unit generators of constant-tail
modules over Z/(2^K), K1..6,S1..9, through all affordable offsets and
horizon drops, with signed multipliers and zero/positive deltas.
Linearity makes those generator tests sufficient for every input in
the corresponding finite module for each tested operation. It also
checks5,376 explicit sparse/dense squares,2,460 odd-division cases,
and an independent3,658-partition calculation of a(50). Both normal
and optimized runs pass with identical outputs.

```python
from itertools import combinations
from collections import Counter
from fractions import Fraction
from math import comb, factorial, gcd, prod
import json

def need(ok,message):
    if not ok: raise RuntimeError(message)

def expanded(packed,S,K):
    if S==0:
        need(len(packed)==1,'terminal coordinate')
        return {0:packed[0]}
    need(len(packed)==min(S,K),'nonterminal coordinate count')
    return {r:packed[min(r,K)-1] for r in range(1,S+1)}

def step(full,offset,rest,j0,delta,factor,modulus):
    return {r:factor*(1 << (j0+r*delta))*full[offset+r] % modulus
            for r in (range(1,rest+1) if rest else (0,))}

basis_checks=0
for K in range(1,7):
    modulus=1<<K
    for S in range(1,10):
        cap=min(S,K)
        for coordinate in range(cap):
            packed=[int(j==coordinate) for j in range(cap)]
            source=expanded(packed,S,K)
            for offset in range(S+1):
                for rest in range(S-offset+1):
                    if rest==0 and offset==0:
                        continue  # a positive-horizon bank has no Y_0
                    for j0 in sorted({0,K-1,K}):
                        for delta in (0,1,2):
                            for factor in (-3,1):
                                full=step(source,offset,rest,j0,delta,factor,modulus)
                                kept=[full[r] for r in
                                      (range(1,min(rest,K)+1) if rest else (0,))]
                                need(expanded(kept,rest,K)==full,'constant-tail closure')
                                basis_checks+=1
# These module maps are linear over Z/(2^K), so unit-generator tests
# span every constant-tail input, including non-realizable moment banks.
need(step({0:3},0,0,0,0,-3,8)=={0:7},'separate terminal channel')

# Identity on the full product lattice does not impose a constant tail.
bad={1:0,2:1}
truth=step(bad,0,2,0,0,1,2)
pretended=step(expanded([bad[1]],2,1),0,2,0,0,1,2)
need(truth=={1:0,2:1} and pretended=={1:0,2:0},'product-lattice hostile')
# This is permitted in the larger constant-tail module mod4, but not
# in any integer moment image: actual Y_1 and Y_2 have equal parity.
larger=expanded([1,0],5,2)
need(larger=={1:1,2:0,3:0,4:0,5:0} and larger[1]%2!=larger[2]%2,
     'constant-tail module strictly exceeds the moment image')

def moment(support,b,r,modulus):
    return sum(x*(1 << (r*(m-b))) for m,x in support.items())%modulus

commutations=0
for K in range(1,7):
    modulus=1<<K
    for b in range(4):
        support={b:3,b+1:-2,b+3:5}
        S=9
        source={r:moment(support,b,r,modulus) for r in range(1,S+1)}
        for offset in range(5):
            for k in range(4):
                for delta in range(b+k+1):
                    bp=b+k-delta
                    for rest in (0,1,4):
                        if offset+rest>S or offset==rest==0:
                            continue
                        j0=1;factor=-3
                        transported={m+k:factor*(1 << (j0+offset*(m-b)))*x
                                     for m,x in support.items()}
                        direct={r:moment(transported,bp,r,modulus)
                                for r in (range(1,rest+1) if rest else (0,))}
                        require=step(source,offset,rest,j0,delta,factor,modulus)
                        need(direct==require,'sparse/dense square')
                        commutations+=1

# Actual defect history: three past9-cycles yield D(3)=0,D(5)=12.
# Skipping5 drops its coordinate: b=3 becomes b'=0, delta=3.
# The relaxed weight x_3=1 tests the allowed larger moment lattice.
before={r:1 for r in range(1,9)}
after=step(before,0,8,0,3,1,8)
need(after=={r:0 for r in range(1,9)},'mask-drop gauge rescaling')
need(after!=before,'mask deletion is not an unscaled identity')

# Terminal Y_0 is not identified with the positive-index tail.
formal={1:1,2:-1}
need(moment(formal,1,0,2)==0 and moment(formal,1,1,2)==1,
     'terminal zero is not the tail')
support={3:8}
transported={m+1:(3*x)%8 for m,x in support.items()}
need(transported=={4:0} and len(transported)==1,'zero-residue semantic support')

# Odd division costs no extra binary bits.
odd_checks=0
for K in range(1,7):
    modulus=1<<K
    for d in range(1,20,2):
        for quotient in range(-20,21):
            numerator=d*quotient
            need((numerator%modulus)*pow(d,-1,modulus)%modulus
                 ==quotient%modulus,'odd-unit division')
            odd_checks+=1

# Independent unmerged Burnside calculation of the precision hostile.
def partitions(n,least=1):
    if n==0:
        yield ()
    for d in range(least,n+1,2):
        for tail in partitions(n-d,d):
            yield (d,)+tail

count=Fraction(0);partition_count=0
for p in partitions(50):
    multiplicity=Counter(p)
    z=prod(d**m*factorial(m) for d,m in multiplicity.items())
    exponent=sum((d-1)//2 for d in p)
    exponent+=sum(gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
    count+=Fraction(1<<exponent,z);partition_count+=1
need(count.denominator==1,'Burnside integral output')
a50=count.numerator
need(a50%(1<<25)==0,'order50 proved baseline')
need(a50%4==0 and (a50>>25)%4==2 and a50%(1<<27)==67108864,
     'absolute/normalized precision hostile')
print(json.dumps({'status':'PASS','constant_tail_basis_transfers':basis_checks,
    'sparse_dense_squares':commutations,'odd_division_checks':odd_checks,
    'order50_odd_partitions':partition_count,'absolute_two_bits':a50%4,
    'absolute_27_bits':a50%(1<<27),'normalized_two_bits':(a50>>25)%4,
    'mask_drop_nonzero_tail_to_zero':True,'product_lattice_hostile':True,
    'strict_larger_module_hostile':True,'terminal_Y0_hostile':True},sort_keys=True))

```

