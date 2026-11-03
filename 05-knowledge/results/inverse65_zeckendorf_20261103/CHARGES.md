# Signed charges, quadratic recurrences, and a norm-13 clock

Status: **PROVED** identities for actual source words; **FINITE-EXACT** bank
classification and listed primality checks. Infinite prime occurrence and
complete prime/7p coverage remain **OPEN**.

## C1. Classify the single-symbol conserved charges

For each of the65 matrices solve `P ell = sigma ell`, sigma=+1 or-1, over Q.
Exactly seven letters have a nonzero solution. Each solution space has
dimension1. Normalize ell to primitive integers with first nonzero entry
positive. Let J=v ell at the initial state. The full classification is:

| Letter | sigma | ell | t | d | kappa | Quadratic field |
|---|---:|---|---:|---:|---:|---|
| b |-1|(1,0,-1)|8|-1|8|Q(sqrt17)|
| q1 |-1|(3,0,-4)|18|-9|6|Q(sqrt10)|
| q2 |1|(1,0,-1)|18|5|-12|Q(sqrt19)|
| q4 |-1|(3,2,-4)|18|-9|6|Q(sqrt10)|
| q5 |-1|(1,1,-2)|16|-21|4|Q(sqrt85)|
| r21,r38 |1|(3,-2,-2)|40|13|-6|Q(sqrt43)|

In each case `char_P(X)=(X-sigma)(X^2-tX+d)` and

    H_(k+2) = t H_(k+1) - d H_k + kappa J sigma^k.       (C1)

Proof: `((P^2-tP+dI)e_h)=kappa ell`, checked exactly; multiply by `vP^k`.
The shifted sequence `(1-t*sigma+d)H_k-kappa J sigma^k` obeys the homogeneous
quadratic recurrence. The sign phase must be retained when sigma=-1.
The characteristic factors and nullspaces are independently reproducible
with rational linear algebra in `charges.py`.

This classifies single-symbol rational charges, not all charges preserved by
arbitrary longer products. The maximal listed subalphabets sharing one charge
line are {b,q2}, {q1}, {q4}, {q5}, and {r21,r38}. In the first, b reverses
the charge and q2 preserves it. The other58 symbols have no such charge.

Thus17 and19 here are discriminant fields of actual transition matrices,
not an identification of a prime list with an ear spectrum. The quadratic
recurrence together with its charge/phase is a real map between the objects.

## C2. An exact sevenfold attachment

Use the opposite-sign charge `J=2x+2y-3h` for convenience. Both r21 and r38
have the identical matrix

    [[13,8,10], [20,14,17], [16,11,14]]

and preserve J. Their literal labelled charts differ: the choice at every
occurrence stays in the source word even though the three-coordinate observer
cannot see it. The terminal triangle b gives

    H_b=h+4(x+y)=7h+2J.

If J=7r, put p=h+2r; the actual terminal attachment has exactly7p paths.
Two useful seeds are:

| Seed | Initial (h,x,y) | J | r | First p values |
|---|---|---:|---:|---|
| r29 |(15,8,11)|-7|-1|13,529,21001,833173,33053917,1311325441|
| r43 |(11,9,11)|7|1|13,501,19861,787917,31258477,1240096149|

The processing word is `seed (r21 or r38)^k b`, on exactly `8+5k` vertices.
Source order is the reversed word. All literal graphs are Boolean; the
multi-component source gives a strong target: the preserved spine and the
reversed nonadjacent last-to-first cross edge form a directed Hamilton cycle
(all these targets have at least8 vertices). Choosing r21 throughout fixes a
canonical representative.

These families satisfy

    h_(k+2)=40h_(k+1)-13h_k+6J,
    p_(k+2)=40p_(k+1)-13p_k-10r.                         (C2)

Their ordinary generating functions are

    sum p_k z^k = (13-4z+z^2)/((1-z)(1-40z+13z^2))     [r=-1],
    sum p_k z^k = (13-32z+9z^2)/((1-z)(1-40z+13z^2))   [r=+1].

The H-weighted vertex-clock series is `7u^8 P(u^5)`. It sums Hamilton counts;
it does not enumerate isomorphism classes or certify that its coefficients
exhaust all target values.

Both p sequences increase strictly. Their successive differences satisfy
`D_(k+1)=40D_k-13D_(k-1)` with positive first differences and D1>D0; induction
gives continued growth. Consequently sequence membership has a uniform inverse:
iterate until equality or overshoot, then output the seed, k copies of r21,
and b. The number of steps is O(log p). This inverse is restricted to these
two explicit lanes; the general bank inverse is in INVERSE.md.

Trial division certifies13,499,1009,19861,21001 prime. Subset DP independently counts
the actual graphs with H=91,499,147007=7*21001,139027=7*19861, the latter two
at18 vertices. As a hostile test,529=23^2. For the r=+1 family, every odd k
has p_k divisible by3 and p_k>3, by (C2) modulo3. No prime-producing claim
follows from the recurrence alone.

## C3. The quadratic norm carries a clock, not an orbit certificate

Set `U_k=13h_k+3J=13p_k-5r`. Then

    U_(k+2)=40U_(k+1)-13U_k,
    Q(U,V)=V^2-40UV+13U^2,
    Q(U_(k+1),U_(k+2))=13 Q(U_k,U_(k+1)).               (C3)

The roots are `20 +/- 3sqrt43`, of norm13. A backward clock step is
`(U,V)->((40U-V)/13,U)`, requiring integrality and eventual arrival at the
chosen seed. Merely satisfying a norm equation is not sufficient orbit
membership. The retained consecutive pair, divisibility, and seed are the
arithmetic counterparts of ports, legal predecessors, and actual sources.

## C4. Every modulus, including the singular prime13

For modulus m the p-pair evolves by the affine map

    (a,b) -> (b,40b-13a-10r) mod m.                     (C4)

If gcd(m,13)=1 it is a permutation, so every initial pair is periodic with
no transient. For each of these two seeds, modulo `13^e` (e>=1), the transient
has **exactly e steps**, followed by period **exactly13^e**. In that periodic
part the first coordinate visits every residue once.

Proof. Work modulo `13^(e+1)` with U=13p-5r. The roots alpha,beta of
`X^2-40X+13` lift uniquely from1,0 modulo13, since the derivatives are units.
The first lifts are alpha=27 and beta=13 modulo169. Thus

    v13(alpha-1)=1, v13(beta)=1, alpha-beta is a unit.

Write `U_k=A alpha^k+B beta^k`. Here A is a unit and v13(B)=1 for both seeds.
Indeed their `(U0,U1)` are `(174,6882)` and `(164,6508)`; the quantities
`(27U0-U1)/13` are1 and-4 modulo13, respectively, and U0 is a unit.
The change from consecutive U coordinates to the two eigencoordinates has
unit determinant. Equality of consecutive p pairs modulo13^e is equivalent
to equality of consecutive U pairs modulo13^(e+1).

The beta coordinate can repeat only after it vanishes: since beta^L-1 is a
unit for L>0, this first happens at k=e. The alpha coordinate has period13^e:
binomial lifting gives `v13(alpha^L-1)=1+v13(L)`. This proves the exact tail
and pair period. Once k>=e, the same valuation, after division by13, shows
that `p_(k+L)-p_k` has valuation v13(L), proving the complete residue visit.

For arbitrary `m=13^e n`, gcd(n,13)=1, CRT now gives tail e and period
`lcm(13^e, period modulo n)`, with the convention e=0 gives no transient.
This is a uniform compatible-clock theorem; the coprime period is computed
by finite pair iteration, not claimed to have a single simple formula.

| Modulus |2|3|5|7|11|13|17|19|23|27|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| r29 period |1|1|8|6|120|13|16|18|33|6|
| r43 period |1|2|8|6|120|13|16|18|33|18|

Only the13 column has a transient (length1). These exact clocks expose the
prime3 obstruction and preserve the distinct phases at27. Exhaustive checks
through13^3 and mixed moduli13*11,169*3,169*19 agree with the all-e proof.
