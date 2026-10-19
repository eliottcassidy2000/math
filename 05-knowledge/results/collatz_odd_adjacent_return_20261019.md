# Odd-lower expanding adjacent returns and bounded-head reset cylinders

**Status, 2026-10-19:** PROVED an infinite family of expanding, oriented
adjacent returns whose odd lower source has no earlier descending prefix.
This disproves a proposed structural screen, not a theorem in the repository.
PROVED a general bounded-head reset lemma and two all-L families of actual
adjacent-return descent receipts. Universal entry into these families is OPEN.

## 1. A genuinely first-descent-free adjacent return

Use the Terras map T(x)=x/2 for even x and T(x)=(3x+1)/2 for odd x. Put

    P=1162261467=3^19, Q=1073741824=2^30,
    C=4912087879, D=P-Q=88519643,
    r=815464443, c=C/D.

Consider the two length-30 parity words

    lower: 110110110110110101011011101010
    upper: 001111111111111101010110000100.

Both have 19 odd steps. Their respective affine carries are

    C=4912087879,
    C_upper=4823568236=C+Q-P.

Their exact lawful cylinders are respectively r modulo Q and r+1 modulo Q.
Consequently every positive source x=r+tQ, t>=0, satisfies

    T^30(x)=K(x)=(Px+C)/Q,
    T^30(x+1)=K(x)+1.                             (1.1)

This is an **oriented adjacent return**: the two independently guarded
actual trajectories remain ordered at the return. The lower source is odd,
and K(x)>x because P>Q and C>0. For example,

    (815464443, 815464444)
        -> in 30 steps -> (882691615, 882691616).  (1.2)

The upper source is even and descends immediately. The next result concerns
the odd lower source, which is the endpoint relevant to a least-counterexample
screen; no first-descent-free assertion is made about the upper source.

### Every lower prefix expands

For a parity prefix of length t and odd count a, its affine map is
(3^a x+C_t)/2^t, where C_t>=0. After an initial run of odd steps, minima of
its multiplier occur at the ends of runs of even steps. For the lower word,
these pairs (t,a) are

    (3,2), (6,4), (9,6), (12,8), (15,10), (17,11),
    (19,12), (22,14), (26,17), (28,18), (30,19).

Every pair satisfies 3^a>2^t. The smallest multiplier over all 30 nonempty
prefixes is 3^12/2^19=531441/524288>1. Thus every nonempty actual prefix of K
is strictly above its positive starting source. A complete K block also
increases its source. Induction proves the same assertion for every
nonempty prefix of every lawful K^L.

In particular the following proposed screen is **false**:

> Every oriented adjacent return with odd lower source either contracts
> that source or has an earlier descending prefix on that source.

Equation (1.1), with the prefix inequalities above, supplies a whole infinite
cylinder of counterexamples, and lawful repetitions supply arbitrarily long
ones. This is not a claim that these sources never subsequently descend.
No minimal-length claim for this counterexample is made.

For comparison, the earlier
[switching-height packet](collatz_switching_height_certificate_20261019.md)
used an expanding adjacent return with even lower source, so it was already
screened by first descent. The
[maximal-run reset packet](collatz_maximal_run_reset_20261019.md) supplied a
first-descent-free odd driver returning the relation (3u+2,u). K supplies
both features together: odd-driver survival and an actual adjacent return.

## 2. Exact compression of a maximal expanding K run

For any positive integer x in K's cylinder, write

    Dx+C=2^(30L+s) m,
    m positive odd, L>=1, 0<=s<30.                (2.1)

Then the maximum number of consecutive complete K blocks is exactly

    L=floor(v_2(Dx+C)/30),

and their actual exit is

    y=K^L(x)=[P^L 2^s m-C]/D.                    (2.2)

Its remaining K rank is exactly s. In particular another full K block is
not lawful. The companion satisfies T^(30L)(x+1)=y+1 throughout this run.

**Proof.** Let

    C_L=C(P^L-Q^L)/(P-Q).

The word K^L has affine numerator P^L x+C_L and modulus Q^L. Since D is odd,
its lawful cylinder is equivalently Dx+C=0 modulo Q^L. This follows either
from the parity-word recurrence or from

    D(P^L x+C_L)=P^L(Dx+C)-C Q^L.

Hence its exact repetition budget is (2.1). The centered identity
DK(x)+C=(P/Q)(Dx+C) gives (2.2) and the residual rank. Integrality follows
also from P=Q modulo D. The separate upper parity-word guard in (1.1)
propagates by induction. ∎

The pair (y,s) stays unchanged as the run is traversed; only the number of
remaining blocks drops. Its computation uses integers of bit length
O(log(x+1)), since 30L<=log_2(Dx+C). This is an exact maximal-run compression,
not a freely imposed switching policy.

## 3. General bounded-head reset lemma

The following statement is independent of the particular K word.

**Theorem.** Let U and W be lawful affine branches of finite Terras parity
words, with

    U(x)=(Px+C)/Q,   P=3^a>Q=2^n, C>0,
    W(x)=(Ax+B)/q,   A=3^b<q=2^ell, B>=0.

Write their centers as -c and p, where

    c=C/(P-Q)>0,   p=B/(q-A)>=0.

Fix any L>=1 and choose

    j_0=min {j>=1 : A^j P^L < q^j Q^L}.          (3.1)

For j=j_0, set

    C_L=C(P^L-Q^L)/(P-Q),
    B_j=B(q^j-A^j)/(q-A),
    A'=A^j P^L,   M=q^j Q^L,
    C'=A^j C_L+Q^L B_j,   Delta=M-A'>0,
    rho=-C'(A')^(-1) modulo M,   0<rho<M.        (3.2)

The full lawful cylinder of the word U^L W^j is x=rho modulo M. On its
positive lifts x=rho+tM, t>=0, the endpoint is smaller than x if and only if

    Delta x>C'.                                 (3.3)

Its exact first successful lift is

    t_min=max(0, floor((C'-Delta rho)/(Delta M))+1),

and, uniformly in L,

    t_min <= ceil(c+p).                           (3.4)

Thus **at most ceil(c+p) initial lifts are excluded**. Every later lift has
a lawful, explicit smaller-source receipt. The same bound holds for any
j>=j_0, not just the minimal choice.

**Proof.** The usual parity-word carry recurrence gives (3.2) and its unique
source cylinder. For completeness, the full-word congruence also implies
all its prefix congruences: reducing the final carry modulo a prefix's
power of two leaves the prefix carry multiplied by a power of three,
which is invertible modulo that power of two. The parity-word induction
then gives the actual branch at every step. Since U contains odd steps,
the complete word is not all even and its residue is nonzero.

Its endpoint is (A'x+C')/M, proving (3.3) and the exact strict cutoff.
The centered expression for the composite carry is

    C'=c A'+p M-(c+p) A^j Q^L < (c+p)M.          (3.5)

It is positive, and Delta is a positive integer. If t>=ceil(c+p), then
x=rho+tM>=(c+p)M>C'/Delta. This proves (3.4), without a power-gap estimate
for 2^u-3^v. Fewer than j_0 resets cannot decrease any positive source:
their multiplier is at least one and their carry is positive. ∎

**Sharp ceiling wording.** Take U(x)=(3x+1)/2, with c=1, and W(x)=x/2,
with p=0, at L=j_0=1. The composite is (3x+1)/4 on x=1 modulo 4. Its first
lift x=1 gives equality; all lifts x=1+4t with t>=1 decrease. Therefore
exactly ceil(c+p)=1 initial lift is excluded. "Fewer than ceil(c+p)" would
be false.

**What the theorem does and does not assert.** The cylinder is an explicit
sufficient descent class, and the bound is on its initial lifts, not on all
positive integers. The theorem itself makes no first-stopping-time claim;
that requires additional control of the prefixes. It also does not assert
that an arbitrary U exit enters W's cylinder.

## 4. All-L resets that preserve the adjacent relation

Use the following two genuine oriented adjacent-return branches:

    F(z)=(3z+1)/4,       lower word 10, upper word 01,
       lower cylinder z=1 modulo 4, center p_F=1;

    H(z)=(81z+101)/512,  lower word 110101000,
       upper word 001101001,
       lower cylinder z=43 modulo 512, center p_H=101/431.

Their upper cylinders are respectively 2 modulo 4 and 44 modulo 512.
The upper carries are 2 and 532, namely the respective lower carry plus
modulus minus multiplier numerator. Thus, at both endpoints with their
separate guards checked,

    T^2(z+1)=F(z)+1,
    T^9(z+1)=H(z)+1.                              (4.1)

Apply section 3 to U=K and W in {F,H}. Denote the reset block length by
ell_F=2 and ell_H=9, and write j_0(L,W) as in (3.1). For every L>=1 this
constructs an exact lower cylinder rho modulo M and a corresponding upper
cylinder rho+1 modulo M. Every source in the certified tail satisfies

    (x,x+1) -> after 30L+ell_W j_0 steps -> (z,z+1),
    0<z<x.                                      (4.2)

Both endpoints therefore finish below their own starting source, and the
entire composite remains an oriented adjacent return. Its lower source has
no descending prefix during any of the first 30L steps.

The general bound specializes to

    t_min<=57 for K^L F^j_0,
    t_min<=56 for K^L H^j_0,                     (4.3)

because 55<c<56, 56<c+1<57, and 55<c+101/431<56.
There is no finite exception census hidden in (4.3): it follows uniformly
from the strict carry inequality and the positive integer Delta.

### These cylinders have the exact maximal K-run length L

The relevant center separations are

    v_2(c+1)=1,
    v_2(c+101/431)=4.                             (4.4)

Indeed

    C+D=2*2500303761,
    431C+101D=16*132878147487,

with both displayed cofactors odd. A lawful F block has
v_2(y-1)>=2>1; a lawful H block has v_2(y-101/431)>=9>4. The unequal-valuation
rule therefore gives v_2(y+c)=1 or 4 respectively. By section 2, the source
rank is exactly 30L+1 or 30L+4. The initial K run is maximal, not merely an
arbitrarily truncated expanding segment. Different L values give disjoint
source cylinders, as do the F and H families.

### Exact lower first stopping time for the F tails

For every certified lift in K^L F^j_0,

    first stopping time = 30L+2j_0.              (4.5)

Every K prefix expands. At an earlier completed F block the composite
multiplier is at least one, by minimality of j_0, and its positive carry
rules out descent. Within an F block its odd first step increases the
current input, so no intermediate descent is missed. The last complete
block is the first descent by (3.3).

For H, (4.2) is a valid receipt and an upper bound on first stopping time,
which lies strictly after 30L. No claim that the final H block boundary is
the first descent is made. The upper source x+1 already descends at step 1
in both families; the two endpoint claims must not be conflated.

### Removing the exceptional head with at most one extra reset block

Let r_L be the least positive source in K^L's cylinder and y_L=K^L(r_L).
Let j_all be the least j>=1 for which W^j(y_L)<r_L as a real-affine
comparison. Then

    j_0<=j_all<=j_0+1.                            (4.6)

Every positive source in the full composite cylinder K^L W^j_all descends;
there are no excluded positive lifts.

To prove this, write lambda=P/Q and mu=A/q. Since r_L>=r and
r>3c+4p_W, one has

    (y_L-p_W)/(r_L-p_W)
        < lambda^L (r_L+c)/(r_L-p_W)
        < (4/3)lambda^L.

Both resets have mu<=3/4, so j_0+1 blocks make this ratio strictly less than
one. Fewer than j_0 cannot suffice. Once the affine composite is below the
identity at r_L, its contracting slope puts it below the identity for all
x>=r_L, hence for every source of its lawful composite cylinder. The
comparison at r_L is not an assertion that W is lawful at y_L; lawfulness
is imposed by the complete composite cylinder. No first-stopping-time
identity is asserted for this optional longer word.

## 5. A source-sensitive decoder, not a forced reset policy

For an arbitrary maximal K exit (2.2), the exact available F/H capacities
can be read from the exit rank and odd part m:

- If s=1, then y=1 modulo 4 and

      v_2(y-1)=1+v_2(P^L m-2500303761).

  Exactly floor(v_2(y-1)/2) consecutive F blocks are available.
- If s=4, then

      v_2(y-101/431)=4+v_2(431P^L m-132878147487).

  Exactly floor(v_2(y-101/431)/9) consecutive H blocks are available.
  This number can be zero: shell 4 alone does not force a full H block.
- At every other residual rank, neither F nor H is lawful as the next
  complete reset word, by the center-separation test (4.4).

The decoder compares these available blocks with the least j satisfying

    mu^j(y-p_W)<x-p_W.                            (5.1)

All quantities are exact rationals, and x>p_W. It returns the actual paired
receipt only when capacity suffices. Otherwise it reports an uncertified
exit. Thus the arithmetic controls rank recharge for these explicit source
classes; it does not assume a graph transition which the source might lack.

For example, sources with exactly one F block available can succeed at
L=1 but fail at L=4, where even the composite multiplier after one F already
exceeds one. Sources can likewise have exactly one H block available when
L=24, which cannot suffice. The checker supplies source-derived instances
of these failures, as well as exits with no full adjacent reset available.

## 6. Concrete receipts and replay

A smallest lift of the L=1 F composite is

    lower: 2962948091 -> K 3207214549 -> F 2405410912,
    upper: 2962948092 -> K companion 3207214550
                                  -> F companion 2405410913.

The lower first stopping time is exactly 32. The full lower source cylinder
has modulus 2^32 and residue 2962948091. A smallest lift for L=1,H is

    520506507259 -> K 563417241643 -> H 89134368307,

with upper endpoints exactly one higher. Its receipt length is 39, while
its actual first stopping time is 35, illustrating the necessary distinction
for H.

Replay:

    python3 05-knowledge/experiments/collatz_odd_adjacent_return_20261019.py
    python3 -O collatz_odd_adjacent_return_20261019.py

The companion script uses explicit exception-raising checks, not Python
assert statements. It checks both parity guards independently, all 30 lower
prefix inequalities, the sharp general ceiling example, 22 parameterized
F/H cylinder controls with L up to 128, 66 literal paired receipts, the
maximal-run invariants, exact F first stopping times, the all-positive-lift
refinement, and successful and unsuccessful source-specific exit controls.
There are 219637 explicit checks; normal and optimized JSON receipts match
byte for byte. These finite controls audit the formulas. The all-L proofs
are the arguments above, not an extrapolation from those controls.

## 7. What is still missing

The proposed odd-lower adjacent-return screen is disproved, and the
expanding-to-contracting mechanism now has genuine first-descent-free
adjacent examples with exact guards. Every displayed reset family gives
actual smaller-source receipts, not an assumption about an uncomputed tail.

The remaining coverage problem is explicit: an arbitrary positive source
need not enter K, and an arbitrary maximal K exit need not have enough
lawful F/H capacity. Nothing here proves that all surviving odd sources
enter these cylinders, or that all uncertified exits later descend. The
Collatz conjecture remains open.
