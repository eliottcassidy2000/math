# A reindexed single-module clock faces factorial and phase gates

**PROVED necessary gates; FINITE-EXACT inherited-fibre comparisons with
independent literal controls.** This addresses the reindexing boundary
left open by the earlier one-TT-clock result. It does not exclude all
delayed clocks, two changing modules, selective incidences, or growing
cores. Chan--Pak and the unit-ready405 return remain OPEN.

## 1. The exact clock language and its inherited coefficients

Fix a tournament context Q and replace one vertex uniformly by TT_m.
Keep the marked ready root outside that module. Suppose, for some
integer d>=1 and integer c, its full count is

    H(Q[TT_m])=405*2^(dm+c)+71                           (1)

for every sufficiently large positive m, and the marked root is unit-ready
there. The offset c is unrestricted: it cannot be reduced modulo d
without proving a physical backward deflation of the context.

The inherited module theorem
[M1--M3](../../mrb_seven_cutoffs_20261113/module/PROOF.md)
gives nonnegative integer skeleton coefficients A_j with

    H(Q[M])=sum_j A_j K_j(M),
    K_j(TT_m)=j! S(m,j),
    (j-1)! divides A_j.                                (2)

The factorial is already a divisor of A_j; it must not be counted a
second time in K_j. The coefficient normalization is essential here.

Put q=2^d and alpha=405*2^c. Since

    q^m=sum_(j=1)^q binom(q,j)K_j(TT_m),

finite exponential independence forces

    A1=alpha*q+71,
    A_j=alpha*binom(q,j),  2<=j<=q,
    A_j=0,               j>q.                          (3)

These coefficients follow even from an eventual identity. In particular
the positive top coefficient is A_q=alpha.

## 2. All large exponential rates are excluded without a search

Since A_q is an integer, c>=0. Divisibility in(2) requires

    (q-1)! divides 405*2^c.                            (4)

If d>=3 then q>=8 and the left side has a factor7; the right side
does not. Thus **no d>=3 is possible at any integer offset c**, even
without a readiness requirement. If d2, divisibility by3!=6 further
requires c>=1. The only remaining rates are

    d1,c>=0,       or       d2,c>=1.                   (5)

This is a direct application of inherited excursion-permutation
divisibility, not a new positivity theorem or a numerical rate fit.

## 3. The complete known ready fibres exclude ten early clock phases

For the remaining rates, eventual readiness extends to m1. The module
cannot contain the external ready root or its sole predecessor. Its
root-deleted starting-path observable is another nonnegative K_j
expansion. Being constantly1 eventually forces coefficient K1=1
and all later coefficients0. Strongness of a fixed uniform quotient
does not depend on the positive module size. The count identity also
extends to m1 and m2 by(3).

At m1 the source count is N_(d+c), where N_k=405*2^k+71. At m2 a
single transitive twin must give N_(2d+c). The inherited complete
ready fibres for k1..6 and their exact all-vertex TT2 counts are in
[inflation_descent.json.gz](../../multicore_metagraph_20261113/inflation_descent.json.gz).
They were proved complete using the unique-start inverse and the
published strong-order bound; this lane does not regenerate that census.

Comparing every recorded vertex gives no match in the following phases:

|rate d|offset c|source exponent d+c|required next exponent 2d+c|
|--:|---|---|---|
|1|0,1,2,3,4,5|1,2,3,4,5,6|2,3,4,5,6,7|
|2|1,2,3,4|3,4,5,6|5,6,7,8|

For example, the only canonical odd-phase quartic clock would require
3311 ->13031 by one TT2 inflation. All three ready3311 cores miss that
count at every vertex. Their strong orders are11,13,13; the complete
ready13031 cores have strong order15, also exposing the incompatible
one-vertex growth in this case.

Weak ready sources consist of one strong ready component followed by
a transitive tail. Inflating a tail vertex leaves H unchanged; inflating
the strong component uses the same table. Thus weak padding does not
repair any displayed miss.

Combining the arithmetic and finite gates, a surviving clock of the
exact form(1) would have to lie in one of the still-open regions

    d1,c>=6,       or       d2,c>=5.                   (6)

There is no claim that these remaining charts exist. Conversely, the
finite head does not exclude them. A delayed clock is not equivalent
to an early clock without an actual smaller-context construction.

## 4. Verification and what this changes

Run `python -B 05-knowledge/results/growing_t_returns_20261113/inverse/reindexed_clocks.py`
and repeat with `-O -B`. The program pins the inherited complete table,
compares all its entries in the ten stated phases, and independently
counts the source graph and the closest TT2 candidate in each source
row:54 full literal Hamilton counts. It also checks the Stirling-basis
normalization independently for q2,4,8 and lengths1..9. The receipt
does not claim a fresh census or fresh recount of every inherited clone.

The resulting map is a proposed reindexed geometric growth law -> its
required positive module coefficients -> a factorial obstruction or a
specific finite ready-fibre test. It preserves full H and explicitly
retains the ready root's location. It loses selective-port information
because uniform module incidences are a hypothesis of the map.
The practical boundary is now narrower: rapid reindexing cannot rescue
one fixed module, while delayed dyadic/quartic charts and genuinely
changing cores remain legitimate targets.
