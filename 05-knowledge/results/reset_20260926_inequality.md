# Split-product contraction, reset debt, and a growing first-collision family

**PROVED elementary inequalities and obstructions / FINITE-EXACT controls /
OPEN global reset inequality.** No Collatz proof or novelty claim.

The positive result is a sharp reset inequality. If the two integer
summands at a precision collision each carry between one quarter and
three quarters of the actual integer, immediately splitting the collision
endpoint contracts their product by a factor strictly below3/4. The
constant3/4 is best possible. However, an exact family of canonical
episodes reaches its first collision with an unbounded net increase in
the actual integer. Thus local collision contraction does not supply the
missing amortized bound for a whole episode.

## 1. Inheritance and the clock

Closest proved mechanism: the three-case exact controller in
[the precision note](nextforest_20260926_precision.md), with
`y=q*2^r+u`, q,u positive odd and r>=1. Canonical hostile: the original
core can return in a larger state, as27 does inside91. Corrected near
miss: subtracting precision exponents is not a decreasing resource.
Least-used sidecar: the *two labelled integer summands* before a swap.

Anchor / niche / wildcard: an actual amortized reset inequality;
role-symmetric product geometry; and rigidity of rejoining from
[the seven3 itinerary note, Lemma R](procgen_seven3_20260926_itinerary_strategies.md).
The six live concepts are source value, split product, projective balance,
collision clock, valuation-one growth, and reset imbalance debt.

This note adopts an explicit **immediate-reset convention**: after a
collision gives x=U(y), split x itself into a new canonical triple, unless
x=1. The split costs no orbit step. This differs from(P5) of the inherited
precision note, which first takes an additional step x->U(x). That
additional step must be separately charged; our3/4 theorem cannot be
applied to its two-step reset without doing so.

## 2. A role-symmetric potential with an exact update

Put X=q*2^r and Y=u, so y=X+Y, X is positive even and Y positive odd.
Define

    K(q,r,u)=XY=q*u*2^r,
    p=u/y,              K/y²=p(1-p).

Write3u+1=2^a v, v odd. Away from a collision a=r, the actual odd
division exponent is m=min(a,r). The two updated summands, before
possibly exchanging their roles, are

    3X/2^m,              (3Y+1)/2^m.

Consequently BOTH consume and swap branches obey

    K'/K=(9+3/u)/4^m.                                    (1)

In particular, m>=2 gives K'/K<=3/4, and m=1 gives K'/K<=3.
This inequality treats a swap as an exchange of real summands, avoiding
the false assumption that a newly large precision exponent is free fuel.
It does not make every transition contracting.

The corresponding projective update is exact:

    p'=(3u+1)/(3y+1)        on consume,
    p'=1-(3u+1)/(3y+1)      on swap.                       (2)

Thus the unordered imbalance d=|2p-1| satisfies

    |d'-d| <= 2/(3y+1).                                  (3)

The large exchange of coefficient/core roles is not itself a large
change of the unordered split geometry. In contrast, a collision
destroys the two-integer decomposition and a fresh split may change
that geometry substantially.

## 3. A sharp contraction theorem at balanced collisions

At a collision a=r, let b=v2(3q+v). Both3q and v are odd, so b>=1 and

    x=U(y)=(3y+1)/2^(r+b),       r+b>=2.

Therefore every collision strictly decreases the immediate actual
integer y>1. In particular y=1 mod4.

Call the input split balanced when y/4<=u<=3y/4. As y=1 mod4 and both
summands are integers, their smaller member is at least(y+3)/4. Hence

    K_old >= 3(y-1)(y+3)/16 =: K_*(y).                    (4)

Reset x into any positive integer even/odd split, not necessarily the
canonical one; define K_new=0 at terminal1. An odd integer x has split
product at most(x²-1)/4, while x<=(3y+1)/4. Thus

    K_new <= (9y²+6y-15)/64,

and division by(4) gives

    K_new/K_old <= (3y+5)/(4(y+3)) < 3/4.                 (5)

**Sharpness.** For k>=1 set

    y=(2^(2k+3)-5)/3,
    u=(y+3)/4,       r=1,       q=3(y-1)/8.

These are positive integers with q,u odd, the input is balanced, and
v2(3u+1)=1=r. The actual collision has valuation2 and gives

    x=2^(2k+1)-1.

Its canonical split is2^(2k)+(2^(2k)-1), whose product attains the odd
integer maximum. Equality holds in the first inequality of(5), and the
ratio tends to3/4. Thus a uniform smaller contraction constant is false.

Without balance, the same computation gives the explicit debt bound

    K_new/K_old <= rho(y)*K_*(y)/K_old,
    rho(y)=(3y+5)/(4(y+3)).                               (6)

For instance, putting D=max(1,K_*/K_old) gives a factor at most(3/4)D.
This isolates the price of an unbalanced reset rather than silently
discarding it.

## 4. That debt is genuinely unbounded

For ell>=3 take

    q=3*2^ell+1,       (q,r,u)=(q,2,1).

This is a legal, generally noncanonical representation. It collides and

    x=U(q)=9*2^(ell-2)+1,
    K_old=4q,
    K_new=2^(ell+1)*(2^(ell-2)+1)

under immediate canonical reset. Therefore K_new/K_old tends to
infinity even though x<y. The input has a tiny odd component, and reset
converts that hidden imbalance into product size.

This witness concerns the class of all legal representations. It is not
by itself a theorem about which representations occur under every fixed
canonical-reset policy. Section5 supplies a separate hostile entirely
inside the actual canonical policy.

## 5. An actual first-collision episode can grow without bound

Start at n_H=2^H-1, H even and H>=2, using its canonical split

    (q,r,u)=(1,H-1,2^(H-1)-1).

After j consume steps,0<=j<=H-2, its exact state is

    (3^j, H-1-j, 3^j*2^(H-1-j)-1).                     (7)

Whenever r>=2, the core valuation is1, so the first H-2 transitions
really are consumes. At r=1, the core is2*3^(H-2)-1. Since H-1 is odd,
its next valuation is2, giving the swap

    ((3^(H-1)-1)/2, 1, 3^(H-1)).                       (8)

The new core valuation is1 because H is even. Hence the next transition
is the first collision, whose endpoint is

    F(n_H)=oddpart(3^H-1)
          =(3^H-1)/2^(2+v2(H)).                         (9)

The last equality follows by factoring3^H-1: for H=2^s h with h odd,
the odd h factor contributes no further2, and successive squarings
give v2(3^H-1)=s+2. The episode has exactly H odd steps: H-2 consumes,
one swap, one collision.

For H=4k+2,

    F(2^H-1)=(3^H-1)/8,
    F(2^H-1)/(2^H-1) -> infinity.                       (10)

The input split at the final collision is asymptotically perfectly
balanced: its summands are3^(H-1)-1 and3^(H-1). It is therefore covered
by the locally contracting theorem(5). Nonetheless the complete episode
has unbounded growth. This is the first failed implication in the
proposal that every sufficiently balanced collision pays back its
preceding precision consumption.

Consequently no positive multiple of log n plus a globally bounded
correction can be nonincreasing at **every first-collision episode** of
this immediate canonical policy, even after excluding a finite core.
This follows from one episode with arbitrarily large endpoint ratio;
it does not exhibit an infinite divergent orbit. A variable collection
of collision episodes is outside this obstruction.

This is distinct from a collision being locally descending. For example,
63 reaches its first collision endpoint91, and1023 reaches7381. On the
finite control universe of2047 sources3..4095,329 first-collision
episodes grow; no finite frequency is used in the proof.

## 6. What a coarse amortized account does prove

On any finite segment with immediate resets, let B count noncollision
steps of actual valuation1, and let G count all remaining steps.
For each collision let D be the factor defined after(6). Multiplying
the proved bounds yields

    K_end <= K_start * 3^B * (3/4)^G * product D.          (11)

Equation(11) is an actual inequality, and retains the unbalanced reset
charges. It can certify particular segments when its right-hand side
is small enough. It is not a global estimate of B,G or the reset debt.
It deliberately discards the extra contraction from valuations>2;
the exact noncollision expression(1) should be used when available.

For translation to integer descent, retain endpoint balance. If
B_i=4K_i/y_i², then

    y_end<y_start iff K_end/K_start < B_end/B_start.

For example, when the terminal split is balanced, B_end>=3/4 and
B_start<=1, so K_end/K_start<3/4 is sufficient for actual descent.
Without this sidecar, K reduction alone need not decrease y. A terminal1
is handled directly rather than assigned a fictitious positive product.

The all-ones family shows precisely what remains unpaid in(11): long
valuation-one runs before a strong collision. The unbalanced family in
section4 shows a second independent debt. An improvement must control
both on a source-preserving sequence of selected episodes.

## 7. A simple logarithmic role rank also fails

Consider on the class of **all legal triples**

    V(q,r,u)=alpha log q+beta log u+gamma r log2+h(q,r,u),

where h is bounded and alpha,beta,gamma are nonnegative. There is no
nonzero coefficient triple for which V is nonincreasing on every
noncollision controller step. This scope does not presume that every
legal triple is reached by one chosen reset policy.

First the swaps(q,1,1)->(1,1,3q), with q arbitrarily large, require
beta<=alpha. Next the swaps

    (1,1,(4^k-1)/3) -> (1,2k-1,3)

require gamma<=beta. Finally start at(1,2m,2^m-1) and take m-1
consecutive consume steps. The net change per step tends to

    alpha log3+beta log(3/2)-gamma log2.

Using gamma<=beta<=alpha, this is at least alpha log(9/4), strictly
positive unless alpha=0, which then forces every coefficient to vanish.
The bounded correction cannot absorb the changes diverging in q,k,m.

Nonnegative coefficients are exactly what global boundedness below
requires for this logarithmic form up to bounded h: q,u,r can each tend
to infinity with the other two fixed. Thus simply pricing the three
stored magnitudes separately does not solve their coupled reset cost.
This does not exclude a nonlinear interaction, an unbounded correction,
or an inequality restricted to an appropriately selected reachable set.

## 8. The rejoining analogy has a decisive arithmetic boundary

The seven3 note's Lemma R compares two affine branch words. If their
slopes differ, their equality has at most one rational solution. The
same elementary argument works for3n+1: equality on infinitely many
distinct source values forces identical affine maps, hence equal odd
counts and division counts by unique factorization.

That is useful for checking proposed **all-height families** of rejoined
paths. It is not a theorem that integer paths must pay equal counts on
rejoining. Positive integers are rational, so the lemma gives no guarantee
that a given integer avoids its exceptional set. Indeed, for any known convergent
integer, its path to1 and that path followed by1's cycle are two different
affine words agreeing at that source. The exceptional arithmetic is
exactly where a fixed-integer proof must operate.

The three role cases likewise do not create extra arithmetic choices:
consume, swap and collision are forced by comparing a with r. They do
retain information that a parity-only quotient loses, but the full
integers q,u and the unbounded precision r remain necessary sidecars.

## 9. Reproduction and strongest next target

Run `python -B 04-computation/experiments/reset_20260926_inequality.py`
and repeat with `-O`. The script checks20,480 triples with q odd1..63,
u odd1..127,r1..10;1078 balanced collisions;80 sharp examples;100
all-ones episodes through H=200; and98 unbalanced collision examples.
It separately records all first-collision episodes for odd sources
3..4095. Every test uses explicit exceptions, not optimizable assertions.

[Script](../../04-computation/experiments/reset_20260926_inequality.py)
and [output](reset_20260926_inequality.out).

The strongest surviving question is whether a computable, variable-depth
selection of **several** collision episodes can bound the accumulated
valuation-one growth and the imbalance factors D by their actual enclosing
boundary, while retaining the endpoint balance needed for integer
descent. A first-collision horizon or an uncharged ternary role label
does not discharge that question. The displayed inequalities and two
independent hostile families make its obligations explicit.
