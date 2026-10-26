# Congruence-stable unbounded resets: independent review

2026-10-25. **REVIEW PASS / PROVED elementary strengthening** of the reset
family in `README.md` section4. **FINITE-EXACT:** the independent control
below imports no producer and passes normally and with `python -O`.
The statement concerns the specified Mahler/Antihydra odd-return map,
not a direct theorem about the classical plus-one Collatz map or the
reachability of these sources from the blank configuration.

## Exact quantifiers

For odd positive integer x>1 use

    d(x)=v2(x-1),       k(x)=v2(3x-1),
    U_F(x)=3^(k(x)-1)*(3x-1)/2^k(x).

The guarded odd block requires credit a>=1 and updates the credit to
a+2k(x)-3. Fix any such a and any positive integer modulus

    M=2^m*q,       m>=0, q>=1 odd.

Choose an integer R0>=5 with R0=5 modulo6 and R0+2>=m. Put

    L=lcm(6,ord_(9q)(2)),
    R_t=R0+tL,       x_t=(2^(R_t+2)+7)/9,      t=0,1,2,... .

Then every x_t is an ordinary positive odd integer, all x_t have the
same residue modulo M, and

    d(x_t)=1,       k(x_t)=2,
    U_F(x_t)=1+2^R_t,       d(U_F(x_t))=R_t.

Thus the next depth is unbounded while a,d,k and the current residue
modulo M are fixed. The fixed residue is the one constructed by the
family; the statement does not assert this for every prescribed residue
class. For all t with R_t>=a+2, the reset is legal and the following
one-valuation run reaches credit0 and halts.

The qualifier a>=1 matters: at a=0 the native odd-source rule halts
before the reset, although the unguarded scalar formula still exists.

## Divisibility and the exact modulus

Because R_t=5 modulo6, 2^(R_t+2)=2 modulo9, proving divisibility of the
numerator by9. The quotient is odd. Direct substitution gives

    x_t-1=2*(2^(R_t+1)-1)/9,
    3*x_t-1=4*(2^R_t+1)/3.

The quotients multiplying2 and4 are odd integers, giving d=1,k=2 and
the displayed return value.

For t>=0,

    x_t-x_0 = 2^(R0+2) * (2^(tL)-1)/9.             (C1)

The period ensures 2^(tL)-1 is divisible by9q. After division by9,
the remaining factor is therefore divisible by q. The initial exponent
R0+2>=m supplies the power2^m. Since q is odd, this proves divisibility
of the difference by M.

For t>0, (2^(tL)-1)/9 is odd. Therefore the dyadic part is sharper:

    v2(x_t-x_0)=R0+2.                              (C2)

The chosen threshold is consequently necessary as well as sufficient
for these positive-period differences to vanish modulo2^m.

Since9 divides9q, the order of2 modulo9, which is6, divides
ord_(9q)(2). Thus the lcm in the definition is redundant but correctly
records the phase restriction. Euler's theorem supplies the explicit
alternative period L=phi(9q); it is also a multiple of6. Computing an
exact multiplicative order is not necessary for existence.

The factor9 in the odd modulus is essential to this argument. With
q=3, using only lcm(6,ord_3(2)) would choose6. Yet R=5,11 give
x=15,911, whose residues modulo3 are0 and2. The valid order is taken
modulo27. The same pair also shows the sharp dyadic boundary: their
difference is896, of valuation7, so they differ modulo256 when the
condition R0+2>=8 is violated.

## Halting after the reset and the permitted conclusion

At fixed a>=1 the k=2 reset raises credit only to a+1. For
0<=j<=a+1 and R>=a+2, the ensuing k=1 run is exactly

    credit_j=a+1-j,
    x_j=1+3^j*2^(R-j).

Before each of those a+1 steps, R-j>=2, so k=1 and the positive-credit
guard holds. At j=a+1 the credit is0, at an odd source of depth
R-a-1>=1. The next block is forbidden and the native rule halts. At
the boundary R=a+2, this is precisely the first possible new reset
source; it halts before that reset is taken.

For any finite list of fixed congruence moduli, take their least common
multiple as M. The same construction fixes every listed observation.
Consequently no finite-valued bound on next depth depending only on
(a,d,k) and those fixed current-integer residues can hold uniformly over
all legal positive-integer starting configurations.

This is not a statement about arbitrary additional features, adaptively
growing moduli, or an observation that retains the exact integer. Nor
does it assert that any x_t occurs on the specified blank orbit starting
at27. It strengthens the information-loss hostile for proposed uniform
finite-modulus reset bounds; it does not decide that fixed orbit.

## Independent finite control

The single executable block below checks all odd q<=51, m=0..16, and
t=0..4; both exact-order and Euler periods; the two failed shortcuts;
multiple simultaneous congruences; and the guarded halt for credits1..12
whenever the stated threshold holds. It calculates every valuation and
odd return directly and imports no earlier implementation.

From the repository root, extract this block from this file and execute
it with Python; the same code also runs under `python -O`.

```python
from math import lcm
import json

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def v2(n):
    need(n>0,'positive valuation argument')
    return (n & -n).bit_length()-1

def phi(n):
    answer=n; remaining=n; p=2
    while p*p<=remaining:
        if remaining%p==0:
            answer=answer//p*(p-1)
            while remaining%p==0: remaining//=p
        p+=1
    if remaining>1: answer=answer//remaining*(remaining-1)
    return answer

def order2(modulus):
    value=2%modulus; order=1; bound=phi(modulus)
    while value!=1:
        value=2*value%modulus; order+=1
        need(order<=bound,'multiplicative-order search terminates')
    return order

def x(R):
    value,remainder=divmod((1<<(R+2))+7,9)
    need(remainder==0 and value>1 and value%2,'positive odd integral source')
    return value

def odd_return(value):
    numerator=3*value-1; k=v2(numerator)
    return k,3**(k-1)*(numerator>>k)

families=values=halts=0
for q in range(1,52,2):
    order=order2(9*q); L=lcm(6,order); euler=phi(9*q)
    need(L==order and euler%L==0,'phase and Euler periods')
    for m in range(17):
        R0=5
        while R0+2<m: R0+=6
        modulus=(1<<m)*q; first=x(R0)
        need(x(R0+euler)%modulus==first%modulus,'Euler-period congruence')
        for t in range(5):
            R=R0+t*L; source=x(R)
            need(source%modulus==first%modulus,'same complete residue')
            need(v2(source-1)==1,'old depth fixed')
            k,returned=odd_return(source)
            need(k==2 and returned==(1<<R)+1 and v2(returned-1)==R,
                 'exact reset and new depth')
            if t:
                need(v2(source-first)==R0+2,'sharp dyadic threshold')
            values+=1
            for a in range(1,13):
                if R<a+2: continue
                credit=a+1; at=returned
                for j in range(a+1):
                    need(credit>0 and v2(at-1)>=2,'ones-block source guard')
                    next_k,at=odd_return(at)
                    need(next_k==1,'actual following valuation is one')
                    credit+=2*next_k-3
                need(credit==0 and at==1+3**(a+1)*(1<<(R-a-1)),
                     'exact guarded halt state')
                halts+=1
        families+=1

need(x(5)==15 and x(11)==911,'small family anchors')
need(lcm(6,order2(3))==6 and x(5)%3!=x(11)%3,
     'reducing before division by nine loses the odd residue')
need(x(5)%256!=x(11)%256 and v2(x(11)-x(5))==7,
     'insufficient dyadic stabilization really fails')

for observations in ((4,9,5),(8,3,25,7)):
    M=lcm(*observations); m=v2(M); q=M>>m
    R0=5
    while R0+2<m: R0+=6
    L=lcm(6,order2(9*q)); residues=tuple(x(R0)%d for d in observations)
    need(all(tuple(x(R0+t*L)%d for d in observations)==residues for t in range(5)),
         'simultaneous fixed congruences')

print(json.dumps({'status':'PASS','odd_q_range':[1,51],'dyadic_exponents':[0,16],
                  'period_indices':[0,4],'congruence_families':families,
                  'exact_reset_values':values,'guarded_halt_checks':halts,
                  'credit_range':[1,12],'simultaneous_modulus_lists':2,
                  'hostiles':['q instead of 9q','R0+2 below m'],
                  'scope':'fixed modular observations across different starting integers; no blank-orbit claim'},indent=2))
```
