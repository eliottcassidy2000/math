# Prime-power lift delay, full quadratic rings, and invariant tournaments

2026-10-15. **PROVED:** all-depth orbit counts for odd-order norm-one
metallic units, in both split and inert unramified quadratic rings; the
split full-ring remark of THM-4656 is thereby proved. **FINITE-EXACT:**
61 complete ring censuses, 2,442,949 elements, with direct orbit and
independent Burnside counts. This is an extension of the existing orbit
dictionary, not a new claim of discovering classical order lifting.

## 1. The coordinate forgotten by reduction modulo p

Use the setting of
[THM-4656, wall-sun-sun-golden-orbit-explosion](../../01-canon/theorems/THM-4656-wall-sun-sun-golden-orbit-explosion.md):
mu^2-a mu-1=0, a>=1, D=a^2+4, p an odd prime not dividing D.
Let R_k=Z[mu]/p^k, pi=ord_(R_1)(mu), m the odd part of pi,
and g=mu^(pi/m). Since N(mu)=-1, pi is even and N(g)=1.
Set

    h=v_p(g^m-1)=v_p(mu^pi-1),

where valuation means the minimum valuation of the two integral-basis
coefficients. It is finite: the real value mu^pi is greater than1.

**Theorem U1 (PROVED).** For every k>=1:

    ord_(R_k)(g)=m p^max(0,k-h),                               (1)

and every nonzero element of content p^(k-l), 1<=l<=k, has precisely
`m p^max(0,l-h)` elements in its g-orbit. Consequently the number of
nonzero orbits is

    N_k=sum_(l=1)^k (p^(2l)-p^(2l-2))/(m p^max(0,l-h)).         (2)

Negation pairs these orbits without a fixed orbit. The number of
g-invariant Cayley tournaments on the additive group of R_k is

    2^(N_k/2).                                                (3)

The vertices are elements of R_k. The binary relation is x->y iff
y-x belongs to a chosen union S of g-orbits. There are no ties for
distinct vertices because exactly one of each paired orbit O,-O is
chosen. The count is of labelled Cayley tournaments with the specified
g symmetry, not of isomorphism classes or all tournaments.

**Proof.** In the p-adic completion, g^m=1+p^h u. In the inert case
the algebra is the unramified quadratic DVR and u is a unit. In the
split case it is Z_p x Z_p and g=(v,v^-1). The valuations of v^d-1
and v^-d-1 are equal for every d, since they differ by a unit factor.
Thus u is a unit in both coordinates. This is the missing argument
for zero divisors in the earlier full-ring remark.

For odd p, binomial expansion gives
`v_p((1+p^h u)^e-1)=h+v_p(e)`: prove first for e coprime to p,
then for e=p, and iterate. For d not divisible by m, g^d-1 is a
unit already modp: in the split case one component being1 would
force the reciprocal component also to be1. Therefore
g^d-1 vanishes to depth l exactly when m divides d and
h+v_p(d/m)>=l. This proves (1).

Write a nonzero x=p^(k-l)y, with y primitive. Inert primitive y is
a unit. Split primitive y has at least one unit coordinate, and
g^d-1 has the same valuation in both coordinates. Hence
`g^d x=x` iff g^d=1 modp^l in either case. The shell contains
p^(2l)-p^(2l-2) elements, proving (2). If g^j x=-x, the orbit's
odd order implies g^j acts identically on x after squaring and
then inverting2 modulo that order; hence x=-x, impossible for
nonzero x at odd p. This proves (3). ∎

At k=2, h=1 gives N_2=(p^2-1)(p+1)/m; h>=2 gives
N_2=(p^2-1)(p^2+1)/m. These are exactly the previous inert
formulas and the previously finite-checked split full-ring formulas.
The same leading modp orbit can persist for several levels before
growth starts; the integer h records that delay. Reduction modp
forgets h, just as projectivizing a clock matrix forgets its scale.
This is a precise information-loss analogy, not a conjugacy with Collatz.

### Why norm one matters in the split case

For a general unit (u,v) in Z/p^k x Z/p^k, the component orders and
lift delays may differ. For example at p=7, k=1, the odd-order unit
(2,1) has orbit length3 on (1,0) and length1 on (0,1), although both
elements are primitive. The content-only shell formula fails. The
reciprocal condition (u,u^-1) is the additional coordinate that makes
U1 valid. Ramified primes and p=2 are explicitly outside the theorem.

## 2. Non-generous and Wagstaff controls

For a scalar odd-order g modulo an odd prime p, with order m and
h=v_p(g^m-1), the same proof gives

    N_k^scalar=sum_(l=1)^k (p^l-p^(l-1))/(m p^max(0,l-h)).        (4)

The scalar invariant tournament count is 2^(N_k^scalar/2).
A non-generous prime is one whose least positive primitive root modp
is not primitive modp^2; see the defining
[OEIS entry A055578](https://oeis.org/A055578/internal) and
[Kerr--McGown--Trudgian](https://arxiv.org/abs/1908.11497).
For odd p this is exactly a lift delay h>=2 for that selected integer
base. The adjective depends on the least representative; it is not
an obstruction to choosing some primitive root modp^2.

The audit verifies p=40487 is prime, its least primitive root is5,
and 5^(p-1)=1 modp^2 but not modp^3. Its delay is exactly2.
Taking g=25 removes the 2-part and puts this example directly in (4).
No density or completeness assertion about non-generous primes is used.

A Wagstaff prime is W_q=(2^q+1)/3 for odd prime q; definition:
[A000979](https://oeis.org/A000979). For q>3,

    ord_(W_q)(2)=2q,  ord_(W_q)(4)=q,
    #4-invariant circulant tournaments = 2^((W_q-1)/(2q)).     (5)

Indeed 2^q=-1 modW_q, so the order divides2q but notq; order2
would force W_q=3, excluded. This proves (5) by the scalar orbit
argument. At q=5, W=11 and the exponent is1: the two choices are
Paley and its reverse. At q=7, W=43 and the exponent is3: there
are eight choices. Thus the 11 example does not justify identifying
every Wagstaff symmetry with the full Paley symmetry. The exceptional
q=3 gives W=3 and ord_3(2)=2, not6.

## 3. Verification and extension boundaries

The checker enumerates every element and orbit for a=1..4,
p in {3,5,7,11,13,17,19}, excluding p|D, at k=1,2;
also k=3 for p<=7 and the split silver exception a=2,p=31,k=2.
It independently checks the Burnside fixed-point count and confirms
no orbit is fixed by negation. Exceptional examples:

| a | p | k | nonzero full-ring orbits |
|---|---|---|---:|
| 2 | 13 | 2 | 4080 |
| 2 | 31 | 2 | 61568 |
| 4 | 3 | 2 | 80 |
| 4 | 3 | 3 | 296 |

The earlier count64 at a=2,p=31 is on the scalar split factor;
61568 is on the full rank-two ring. They count different vertex sets.

```sh
python3 04-computation/experiments/unit_lift_clock_20261015a.py 05-knowledge/results/unit_lift_clock_20261015a.json
```

[Exact output](unit_lift_clock_20261015a.json).
The proof extends to any norm-one element of an unramified quadratic
Z_p-algebra having odd residue order and finite h. It does not require
a=1 or2. A broader split theorem needs both component valuations and
both element-content coordinates; a single content shell is insufficient.
