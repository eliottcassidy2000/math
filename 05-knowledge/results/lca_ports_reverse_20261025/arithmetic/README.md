# Antihydra, the Mahler map, and a marked-ancestor barrier

2026-10-25. **PROVED elementary derivations for the stated arithmetic
maps; CITED native macro-rule interface; FINITE-EXACT controls.** No
Antihydra nonhalting or Collatz convergence claim is made. The symbolic
conjugacy preserves a valuation itinerary, not positive-integer realization.

## 1. Inheritance and source boundary

The closest earlier mechanism is July's
[THM-2228, Mahler carry tails and integral stabilization](../../../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md).
It already proves the full binary-shift coding of ceil(3z/2), and separates
formal 2-adic coding from ordinary positive-integer stabilization. The
[Session56 valuation/LCA carrier](../../session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md)
turns v_2(x-y) into low-bit-first common-prefix depth and recovers the
finite repetition fuel of27. The
[guarded debt monoid](../../session56_rmq_lca_connections_20261025/debt-range/RESULT.md)
retains the minimum over *all* internal prefixes. Its operation is
associative but not idempotent: overlapping blocks cannot be counted twice.

The less-used older sidecar is the labelled split, rather than bare tree
shape, in
[THM-3227, residue contact tries and delayed resplitting](../../../../01-canon/theorems/THM-3227-selected-root-residue-contact-trie-primitive-carry-and-delayed-resplitting.md).
Its scope is selected polynomial germs; no germ-to-Collatz map is assumed.
Here the transferred mechanism is specifically prefix depth plus a marked
root and exact affine carry. The hostile is a reset with bounded old
depth and arbitrarily large new depth.

The Antihydra source is
[At the Limits of Computability, version2, section2.2](https://arxiv.org/html/2610.22603v2#S2.SS2),
checked2026-10-25. We use its parameter updates, starting configuration
and halting predicate as **CITED** native inputs. We do not use its
printed negative clock at the b=0 halt boundary; the earlier
[source audit](../../normseven_boundary_20261025/PAPERS_AND_TRANSPORT.md)
records that discrepancy. No new native-machine or Lean proof is claimed.
The classical Collatz parity conjugacy is credited to
[Bernstein--Lagarias, The 3x+1 Conjugacy Map](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/3x-1-conjugacy-map/6975BB4A8C46CF6842217043AAF9EC13).
The formulas below have independent elementary proofs.

Board: labelled prefix depth / affine inverse branches / credit guard /
integer stabilization / restored precision / finite reverse potential.
Anchor: exact fixed-source barriers. Niche: the recovered Mahler translation.
Wildcard: the golden-ratio Haar survival law, separated from one integer.

## 2. A literal integer translation to the old Mahler map

Write the abstract native state as(a,b). Its nonhalting rules are

    b even: (a,b)->(a+2, 3b/2+2),
    b odd, a>=1: (a,b)->(a-1, (3b+3)/2).

At b odd,a=0 the machine halts. Set y=b+4. Both scalar branches become

    F(y)=(3y-(y mod2))/2=floor(3y/2) on integers.     (A1)

More strongly, z=y-1=b+3 satisfies

    z'=ceil(3z/2).                                  (A2)

Indeed F(y)-1=(3(y-1)+((y-1) mod2))/2. This is an
integer-preserving affine conjugacy to THM-2228's map, not merely a
shared symbolic picture. If c=z mod2, the other coordinate obeys

    a'=a+3c-1.                                      (A3)

Starting from(a,b)=(0,4), three even-b steps give(6,23), so the first
odd-y state is(a,x)=(6,27). Native step clocks are not identified with
these three macro steps.

At an odd x, put k=v_2(3x-1). The first F step is odd, followed by k-1
even F steps, and the next odd point is

    U_F(x)=3^(k-1)*(3x-1)/2^k,
    a'=a+2k-3.                                      (A4)

The guard a>=1 is required BEFORE this block. Its minimum internal
credit is a-1, even when its net change is positive. For example a=0,
x=7 would have net change+1 under k=2, but the native rule halts
immediately. The scalar odd-return formula is an unguarded extension.

If k_0,...,k_(j-1) are the actual valuation letters from27 and
A_j=sum_(i<j)k_i, then formal odd-return credits are

    a_j=6+2A_j-3j.                                  (A5)

The native orbit survives exactly while these credits remain positive
at odd sources. A zero credit at an odd source is already a halt state.

## 3. Exact shift coding, with the domain stated

Let U_C(x)=(3x+1)/2^v_2(3x+1) on odd2-adic inputs for which the
numerator is nonzero. The inverse branches on odd2-adic y are

    C_k(y)=(2^k*y-1)/3,
    F_k(y)=1/3+(2/3)^k*y,              k>=1.          (A6)

Each is a strict 2-adic contraction of scale2^-k, maps odd inputs to
odd inputs, and has exact next valuation k. For an infinite positive
valuation sequence k and A_0=0,A_j=sum_(i<j)k_i, contraction gives

    Psi_C(k)=-sum_(j>=0) 2^A_j/3^(j+1),
    Psi_F(k)= (1/3)*sum_(j>=0) (2/3)^A_j.             (A7)

Both sums here are 2-adic. They biject the full one-sided shift on
positive integer letters with the corresponding odd points having
infinitely many odd returns. The excluded points are those eventually
reaching the singular numerator: -1/3 on the Collatz side,1/3 on the
F side, and their finite inverse preimages. The accelerated conjugacy
is not asserted on all odd Z_2. The unaccelerated parity systems do
have their full binary-shift conjugacies on all Z_2.

The map H=Psi_F composed Psi_C^(-1) preserves the valuation sequence
and intertwines the two accelerated maps. For a finite word of cost A,
its exact finite-return cylinder is a full residue class modulo2^(A+1),
not merely2^A. The image of the corresponding infinite-word cylinder
is that residue class intersected with the infinite-return domain:
eventual poles must still be excluded. If two infinite
words first differ after r common letters, at next letters k,l, then

    v_2(Psi_*(word1)-Psi_*(word2))
       = A_r+min(k,l),             *=C or F.         (A8)

Proof: the common inverse prefix multiplies differences by2^A_r times
an odd unit. At the first different branch, the two numerators have
different valuations k,l, so their difference has valuation min(k,l).
Thus H preserves the actual2-adic distance. This is the weighted LCA
metric of unary-coded valuation words; the number of common valuation
letters alone is not the digit depth.

The decisive integer-realization hostile is the constant word k=2:

    Psi_C(2,2,...)=1,       Psi_F(2,2,...)=3/5.       (A9)

The conjugacy therefore does not preserve positive integers. A periodic
F valuation word has a rational decoded value in[1/3,1], as follows by
summing its positive real geometric series in(A7); the 2-adic periodic
sum gives the same rational. A finite inverse prefix preserves that
interval. Consequently no ordinary integer x>1 has an eventually
periodic F valuation itinerary. This is also the translated positive-
integer nonperiodicity mechanism already present in THM-2228.

Real growth differs from symbolic recurrence: for x>=3,
U_F(x)=(3/2)^k*(x-1/3)>x. For U_C the multiplier is3/2^k with an
additional positive carry. The Haar valuation mean2 gives contrasting
mean logarithmic multipliers, but is not a law for a specified integer.

## 4. A marked-LCA barrier and an explicit precision-reset family

For odd integer x>1 let d=v_2(x-1). The marked low-bit-first LCA with
the fixed point1 has depth d. The exact number of consecutive k=1
steps before the next k>=2 step is

    L(x)=d-1.                                       (A10)

Indeed k=1 iff x=1 mod4. Along that branch
U_F(x)-1=3(x-1)/2, so its valuation drops exactly one until it reaches1.
During the same run a drops exactly one each step. Thus a-d is
constant during the run. The block hits the halt barrier before a
reset exactly when

    d>a.                                           (A11)

For x=1 set d=infinity; every finite credit is then exhausted.
At every actual odd source, d<=a is necessary for infinite survival.
If the inequality holds at every odd source it is also sufficient:
all credits are then positive. This remains a fixed-source orbit
condition, not an eliminated quantifier over all future states.

The difficulty is regenerated depth at a k>=2 reset. For every
R=5 mod6, R>=5, define

    x_R=(2^(R+2)+7)/9.                              (A12)

This is an odd positive integer with

    v_2(x_R-1)=1,  v_2(3x_R-1)=2,
    U_F(x_R)=1+2^R,  v_2(U_F(x_R)-1)=R.             (A13)

For example x_5=15, x_11=911, and x_(R+6)=64*x_R-49.
The formulas follow by direct substitution; divisibility by9 uses
2^6=1 mod9. At any fixed initial credit a>=1 the reset increases
credit only to a+1, while R is unbounded. For R>=a+2, the following
ones-run exhausts credit and halts. At(a,x)=(6,911), the reset lands
at(7,2049) and halts after seven subsequent k=1 blocks.

Thus there is no bound on regenerated LCA depth from(a,d,k) alone,
even with d=1,k=2 fixed and every scalar source an ordinary integer.
These are different starting configurations, not points asserted to
occur in the blank orbit. They test the information lost by a proposed
finite-feature rank. The blank source27 still requires its own proof.

The failure survives any fixed finite collection of congruence
observations. Combine their moduli into M=2^m*q with q odd. Choose
R0=5 mod6, R0>=5, R0+2>=m, and let

    L=lcm(6, ord_(9q)(2)),       R=R0+t*L, t>=0.    (A13b)

Then every x_R has the same residue modulo M, the same d=1 and k=2,
and arbitrarily large new depth R. This constructs one residue class
for each M; it does not prescribe arbitrary residue values in advance.
Indeed

    x_R-x_R0 = 2^(R0+2)*(2^(t*L)-1)/9.

The numerator's odd factor is divisible by9q, and its power of two
is at least2^m. The factor9 in the order modulus is necessary because
the formula divides by9. At any fixed a>=1 all sufficiently large
members exhaust credit after the reset. This does not address arbitrary
finite-valued features: an observation could be defined using the very
future property at issue. The scope is fixed finite congruence data,
with the exact current a,d,k. The
[independent congruence review](CONGRUENCE_REVIEW.md) gives both omitted-
hypothesis hostiles and an independent finite control.

Every finite surviving valuation prefix also has a positive-integer
extension that halts: append as many k=1 letters as its remaining
credit. The resulting exact finite cylinder contains positive integer
lifts. This shows that no nonempty surviving cylinder can certify
nonhalting for all its integer lifts. It does not exclude a proof that
retains the exact integer and uses an additional global invariant.

## 5. A golden-ratio probability, not a fixed-source decision

Normalized Haar measure on odd Z_2 makes the valuation letters
independent with probability P(k=j)=2^-j, because the cylinder of
cost A has conditional odd measure2^-A. Let

    rho=(sqrt(5)-1)/2.

For an odd source with credit a>=1, the exact Haar probability of
eventually hitting the credit barrier is

    P_Haar(halt | a)=rho^a.                         (A14)

One proof uses full binary parity steps. After the first forced odd
step, the credit is a-1 and the future parity bits are fair independent
bits. Credit increments are+2 or-1 with equal probability, and ruin is
at-1. For the unrestricted walk the bounded stopped martingale
rho^(credit+1) follows from rho^3-2rho+1=0. The strong law gives credit
tending to infinity on nonruin paths, so dominated convergence gives
ruin probability rho^(initial credit+1), hence(A14). Equivalently,
sum_(k>=1)2^-k*rho^(2k-3)=1.

This computes a measure of a shift-language event. It assigns no
nontrivial probability to the truth of nonhalting from the already
fixed integer27. Positive integers form a countable Haar-null set.
The explicit family(A12) also shows why positive drift and small old
LCA depth do not provide a pointwise safety certificate.

## 6. What transfers to the tournament work

The recovered LCA/debt apparatus transports exact fixed-word data:
prefix ownership, a declared minimum, and affine continuation guards.
It does not replace the prime quotient orientations in a tournament's
module hierarchy. The six repaired tournament failures retain those
literal arc ports; their automorphisms force the new ledger guards.

The reverse valuation algorithm in[the counting lane](../counting/README.md)
computes a suffix minimum on a finite acyclic graph whose stage always
decreases. Antihydra and Collatz have no supplied finite terminating
graph of that kind. Taking a minimum over an unproved infinite
continuation is not an implementation of that algorithm. The common
method is to preserve the complete boundary state before optimizing;
the missing unbounded fixed-source certificate stays explicit.

Run check.py under normal and optimized Python. The controls verify
translations, exact finite valuation cylinders, weighted LCA depths,
the reset family and integer halt extensions, and exact quadratic-field
martingale identities. They are not native-machine proofs or an
extrapolation to the blank orbit's infinite future.

The [independent arithmetic review](INDEPENDENT_REVIEW.md) checks the
proofs and supplies a separate Fraction-based implementation, including
actual eventual poles and periodic rational points.
