# Independent arithmetic review: translations, exact domains and reset barrier

2026-10-25. **ACCEPTED / PROVED for the abstract arithmetic maps:** the
formulas and scoped claims in [README.md](README.md). **CITED:** the
native Antihydra parameter interface. **FINITE-EXACT:** the producer
and independent controls below. Neither the proof nor the controls
decide blank-tape halting, nonhalting from27, or Collatz convergence.

## Source and native-interface audit

The primary
[Antihydra section2.2](https://arxiv.org/html/2610.22603v2#S2.SS2)
was read independently. Its two parameter updates, initial A(0,4)
and halt predicate agree with the cited interface. The displayed halt
clock has the previously recorded negative boundary value; the present
arithmetic does not use it. No native transition system, Lean build
or clock proof was reproduced.

The primary
[Bernstein--Lagarias publication](https://doi.org/10.4153/CJM-1996-060-x)
is *Canadian Journal of Mathematics*48(6),1996,1154--1169. Its abstract
states the unaccelerated Collatz parity conjugacy. It is properly
credited as background; the accelerated maps and their precise
domains here receive their own proofs. The publication date is1996,
distinct from Cambridge's later online-publication date.

The earlier
[THM-2228, Mahler three-halves carry tails](../../../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md)
was read. Its ceiling map is exactly the translated map here. Its
ordinary-integer stabilization condition is not discarded by invoking
the formal binary shift.

## Translation, guards and odd-return clock

For b even and odd, direct substitution yields y=b+4 and
F(y)=floor(3y/2). Subtracting one gives the exact affine conjugacy
z=b+3 to ceil(3z/2). Since z and b have opposite parity, the credit
increment is3(z mod2)-1. Three allowed even-b updates from(0,4)
give(2,8),(4,14),(6,23), hence odd scalar27 with credit6.

For an odd scalar x, k=v2(3x-1) is the number of F macro steps until
the next odd scalar: one odd step followed by k-1 even steps. The
scalar factor and credit change in(A4) follow exactly. The smallest
credit inside that block is a-1, so the source guard a>=1 cannot be
replaced by a nonnegative net block change. The corrected finite
control is U_F(7)=15 with k=2; a=0 still halts before that block.
The positive-integer domain keeps every future odd return finite.

A5 follows by summing the guarded credit updates. Its interpretation
requires every earlier odd source to have positive credit; extending
the scalar map past a halt is a formal operation, not a native step.
No equality between macro-step counts and native-machine clocks is used.

## Shift coding and singular domains

Both inverse branches in(A6) map odd2-adic inputs to odd inputs and
have exactly the specified next valuation. Their difference scales
are2^k times odd units. Iterating the inverse branches gives(A7);
the omitted remainder tends to zero because the cumulative cost
A_j tends to infinity. This proves existence and surjectivity onto
points with infinitely many finite odd returns. Different first
valuation letters yield different2-adic distances, proving injectivity.

The excluded sets are exactly the two singular numerators and all
their finite inverse preimages. They are countable. On these points
the unaccelerated parity map continues through zero, but the claimed
infinite accelerated valuation word does not exist. The accelerated
intertwining statement is correctly restricted to the complements
of those sets.

A finite word of cost A gives the full odd residue cylinder modulo
2^(A+1) when only those finitely many returns are required. The image
of its infinite-word cylinder is that residue class intersected with
the appropriate singularity-excluding domain. This distinction was
clarified in the producer during review. The extra bit follows from
applying a scale2^A inverse prefix to the odd coset1+2Z_2.

At the first different letters k,l, the branch constants cancel and
the two residual summands have different valuations, giving min(k,l).
Multiplication by the common inverse prefix adds A_r. Thus(A8)
preserves actual low-bit common-prefix depth, not merely the number
of shared valuation letters.

For the constant letter2, the exact fixed points are1 on the Collatz
side and3/5 on the F side. Hence the conjugacy fails to preserve
positive integers. This is an explicit obstruction, not an absence
of an integer-preservation proof.

For an F period, its inverse affine composition has real slope
(2/3)^A<1 and a unique rational fixed point. The positive real
series gives that same rational, strictly above1/3 and at most1.
Every finite inverse prefix preserves the closed interval[1/3,1].
Consequently no ordinary integer x>1 has an eventually periodic
F valuation itinerary. The Collatz inverse series is not being
interpreted as a convergent positive real series in this argument.

## Fuel, resets and finite-prefix halting lifts

For odd x>1, k=1 iff x=1 modulo4. Along that branch
U_F(x)-1=3(x-1)/2, so d=v2(x-1) drops one per block until reaching1.
There are exactly d-1 such blocks. Credit also drops one per block.
The orbit halts before executing a reset exactly when d>a, including
the equality boundary d=a+1: after a one-blocks the credit is zero
at the reset's odd source. At d=a it reaches that source with
credit1 and can execute it.

For x=1 the repeated letter1 exhausts every finite credit. The
necessary-and-sufficient all-future condition retains the future
orbit quantifier; it is not an effective decision rule for27.

For R=5 modulo6, direct substitution in(A12) gives integrality,
d=1, k=2, next scalar1+2^R and new depth R. The recurrence
x_(R+6)=64x_R-49 and both small examples are correct.
At any fixed initial a>=1 the reset adds only one credit. If
R>=a+2, the next a+1 one-blocks reach a halt state. These are
different initial integers; their occurrence in the blank orbit
is neither asserted nor established.

For a surviving finite valuation prefix with remaining positive
credit c, append c letters1. Its finite cylinder has arbitrarily
large positive-integer lifts, including lifts with the native
scalar at least5. Those lifts reproduce every specified valuation
and first reach credit zero after the appended blocks. This refutes
a universal nonhalting certificate for all integer lifts of a
surviving cylinder. It does not alter the already fixed source27.
A prefix already ending with zero credit needs no extension.

## Haar probability and its exact limit argument

The finite-cylinder measure on odd Z_2 is2^-A, so the valuation
letters are independent with P(k=j)=2^-j; removing the countable
singular set does not change that measure.

For odd source credit a>=1, the first odd macro step is forced and
leaves credit a-1. The remaining binary parity word is a fair
independent sequence. The resulting walk has increments-1,+2.
Hitting-1 is equivalent to encountering the native halt predicate
before that forbidden odd step. It is not a permitted native
negative-credit configuration.

For rho=(sqrt5-1)/2, rho^3-2rho+1=0 makes
rho^(credit+1) harmonic. Stopping at ruin bounds it between zero
and one. The unrestricted walk has positive mean increment1/2,
so on nonruin paths its credit tends to infinity almost surely.
Dominated convergence then identifies the ruin probability as
rho^a. The displayed geometric valuation identity agrees.
This justifies the infinite probability formula; finite martingale
checks alone would not do so.

The probability is for Haar-random2-adic scalars at fixed credit.
Positive integers are a countable null set. Neither rho^6 nor the
finite4096-return blank prefix assigns a probability to the truth
of halting from the specified integer27.

## Reproduction

The producer passed normal and optimized Python with byte-identical
outputs, including3670 translation controls,4092 binary cylinders,
1360 finite valuation lifts,32640 weighted-LCA pairs,1499 fuel
sources,507 reset/credit cases,341 halting extensions and576 exact
finite-horizon martingale balances. Its4096 odd-return blank prefix
is finite data only.

The standalone code below imports no producer functions. It uses
exact rational arithmetic to test actual eventual poles, inverse
words, periodic valuations, interval bounds, credit guards and
integer halting lifts. Normal and optimized runs produce identical
PASS output. The checks are bounded corroboration of the proofs,
not a native-machine verification.

```python
from fractions import Fraction as Q
from itertools import product
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def val(x):
    x=Q(x)
    if not x:raise ValueError('infinite valuation')
    def order(n):
        n=abs(n);return (n&-n).bit_length()-1
    return order(x.numerator)-order(x.denominator)
def inv(mode,k,y):
    return (2**k*Q(y)-1)/3 if mode=='C' else Q(1,3)+Q(2,3)**k*Q(y)
def word(mode,w,y):
    for k in reversed(w):y=inv(mode,k,y)
    return Q(y)
def forward(mode,x):
    v=3*Q(x)+(1 if mode=='C' else -1)
    k=val(v)
    need(k>=1,'odd source')
    return (v/2**k)*(3**(k-1)if mode=='F'else 1),k
def residue(x,M):
    x=Q(x);return x.numerator*pow(x.denominator,-1,M)%M
def macro(a,b):
    if b%2:
        return None if a==0 else (a-1,(3*b+3)//2)
    return a+2,3*b//2+2

translations=returns=0
for a in range(8):
    for b in range(129):
        nxt=macro(a,b)
        if nxt is None:
            need(a==0 and b%2,'halt guard');continue
        ap,bp=nxt;y=b+4;z=b+3
        need(bp+4==3*y//2 and bp+3==(3*z+1)//2 and ap==a+3*(z%2)-1,
             'literal two translations')
        translations+=1
s=(0,4)
for _ in range(3):s=macro(*s)
need(s==(6,23),'first odd state is credit6, scalar27')
need(macro(0,3)is None and forward('F',7)==(15,2),'unguarded block hostile')
for x in range(5,512,2):
    target,k=forward('F',x)
    for a in range(1,9):
        credit=a;y=x;minimum=a
        for _ in range(k):
            need(y%2==0 or credit>=1,'block guard')
            credit+=-1 if y%2 else 2
            minimum=min(minimum,credit);y=3*y//2
        need((Q(y),credit,minimum)==(target,a+2*k-3,a-1),'exact odd return')
        returns+=1

prefixes=[w for L in range(4)for w in product(range(1,4),repeat=L)]
lca=poles=cylinders=0
for mode in ('C','F'):
    fixed2=Q(1)if mode=='C'else Q(3,5)
    fixed3=Q(1,5)if mode=='C'else Q(9,19)
    need(forward(mode,fixed2)==(fixed2,2)and forward(mode,fixed3)==(fixed3,3),
         'periodic tails')
    singular=Q(-1,3)if mode=='C'else Q(1,3)
    for w in prefixes:
        p=word(mode,w,singular)
        for k in w:
            p,observed=forward(mode,p);need(observed==k,'eventual pole prefix')
        need(p==singular,'exact eventual pole')
        try:forward(mode,p)
        except ValueError:pass
        else:raise RuntimeError('eventual pole has infinite itinerary')
        poles+=1
        for k,l in product(range(1,5),repeat=2):
            if k==l:continue
            left=word(mode,w,inv(mode,k,fixed2));right=word(mode,w,inv(mode,l,fixed3))
            need(val(left-right)==sum(w)+min(k,l),'weighted LCA depth')
            lca+=1
        M=2**(sum(w)+1);r=residue(word(mode,w,1),M)
        for j in (3,5):
            x=Q(r+j*M)
            for k in w:
                x,observed=forward(mode,x);need(observed==k,'finite exact cylinder')
            cylinders+=1
need(word('C',(2,),1)==1 and word('F',(2,),Q(3,5))==Q(3,5),
     'conjugacy is not integer preserving')

periods=0
for L in range(1,5):
    for w in product(range(1,5),repeat=L):
        beta=word('F',w,0);alpha=word('F',w,1)-beta;x=beta/(1-alpha)
        need(word('F',w,x)==x and Q(1,3)<x<=1 and val(x)==0,'periodic F interval')
        y=x
        for k in w:
            y,observed=forward('F',y);need(observed==k,'periodic exact valuations')
        need(y==x,'periodic forward cycle')
        A=0;partial=Q(0)
        for k in w:
            partial+=Q(1,3)*Q(2,3)**A;A+=k
        need(x==partial/(1-Q(2,3)**A),'periodic rational equals both sums')
        for prefix in ((1,),(2,3),(4,1,2)):
            y=word('F',prefix,x);need(Q(1,3)<=y<=1,'eventually periodic interval')
        periods+=1

fuel=reset_cases=extensions=0
for d in range(1,17):
    for odd in (1,3,5):
        x=1+2**d*odd;y=Q(x)
        for _ in range(d-1):
            y,k=forward('F',y);need(k==1,'one-run step')
        need(val(y-1)==1 and forward('F',y)[1]>=2,'exact one-run length')
        for a in range(d+3):
            credit=a
            for _ in range(d-1):
                if credit==0:break
                credit-=1
            need((credit==0)==(d>a),'marked ancestor barrier before reset')
        fuel+=1
for R in range(5,126,6):
    need((2**(R+2)+7)%9==0,'reset integrality')
    x=(2**(R+2)+7)//9;xp,k=forward('F',x)
    need(x%2 and val(x-1)==1 and k==2 and xp==1+2**R and val(xp-1)==R,
         'reset formula')
    need((2**(R+8)+7)//9==64*x-49,'reset recurrence')
    for a in range(1,13):
        if R<a+2:continue
        y=xp;credit=a+1
        for _ in range(a+1):
            need(credit>=1,'reset post-block guard')
            y,k=forward('F',y);need(k==1,'reset depletion');credit-=1
        need(credit==0,'reset halts after a+1 one-blocks');reset_cases+=1
for initial in range(1,7):
    for w in prefixes:
        credit=initial;valid=True
        for k in w:
            if credit<=0:valid=False;break
            credit+=2*k-3
        if not valid or credit<=0:continue
        extended=w+(1,)*credit;M=2**(sum(extended)+1)
        start=residue(word('F',extended,1),M)+M
        need(start>=5 and start%2,'positive native scalar lift')
        y=Q(start);a=initial
        for k in extended:
            need(a>0,'first halting boundary not bypassed')
            y,observed=forward('F',y);need(observed==k,'halting integer extension')
            a+=2*k-3
        need(a==0,'finite-prefix extension reaches halt');extensions+=1

# Exact quadratic-field harmonic identity; probability conclusion also
# needs the stopped-martingale and positive-drift argument in the proof.
def mul(x,y):
    a,b=x;c,d=y;return (a*c+b*d,a*d+b*c-b*d)
def power(k):
    x=(1,0)
    for _ in range(k):x=mul(x,(0,1))
    return x
need(power(3)==(-1,2),'rho^3-2rho+1=0')
for a in range(20):
    x=power(a);y=power(a+3);z=power(a+1)
    need((Q(x[0]+y[0],2),Q(x[1]+y[1],2))==z,'martingale harmonic equation')
need(power(6)==(5,-8),'rho^6')
print(json.dumps({'status':'PASS','translations':translations,'odd_return_controls':returns,
  'eventual_poles':poles,'weighted_LCA_pairs':lca,'finite_cylinder_lifts':cylinders,
  'periodic_F_words':periods,'one_run_sources':fuel,'reset_credit_cases':reset_cases,
  'positive_integer_halting_extensions':extensions,'harmonic_cases':20,'quadratic_root_and_power_checks':2,
  'native_clocks_checked':False,'infinite_blank_orbit_claim':False},sort_keys=True))

```
