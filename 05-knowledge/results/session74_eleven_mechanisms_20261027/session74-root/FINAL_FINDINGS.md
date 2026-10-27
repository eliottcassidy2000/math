# Session 74 — Leveraging eleven through operations and retained information

## Executive findings

This session began with the actual inherited `{2,3,11}` theorem, not a search
for numerological matches. The resulting bridges use explicit maps and expose
the information needed to make them compositional. Results below are proved
in the accompanying packets and checked independently. Historical novelty is
not claimed; the golden coding has a classical arithmetic-coding precedent.

### 1. Paley11 has exact complement distances six and fifteen

Let L(T) be half the sorted outdegree L1 distance from `(0,...,n-1)`, and
let beta(T) be the minimum cycle rank of the backward-edge graph over all
vertex orders. For every n>=2,

    beta(T) >= max(0, L(T)-n+2).

In particular every backward-forest frame has L<=n-2. The inherited legal
Hamiltonian-path complement changes L by at most one; full converse and
relabeling preserve it. Paley11 is regular with L=15, so reaching a forest
requires at least six moves and reaching a transitive tournament at least
fifteen. Explicit independently verified paths attain both bounds:

    P11 -- six complements --> forest -- nine more --> transitive.

Every prefix and suffix of this fifteen-move word is geodesic to its stated
endpoint. Along the first six moves, intrinsic cycle debts are exactly
5,4,3,2,1,0. The inherited initial value is beta(P11)=10, so beta itself is
NOT a one-step Lipschitz potential. The score-certified lower envelope is the
valid clock bound. The final forest has ten edges; L<=n-2 does not reduce the
forest feedback-edge ceiling from n-1 to n-2.

Scope: fixed order11, actual named legal paths, distance to the forest or
transitive CLASS. These are not all-order accessibility results and do not
specify distance to every separately prescribed labelled transitive tournament.
The discarded L=1 trap has243 Hamiltonian paths and no L-decreasing move,
showing that descending score greedily need not finish.

### 2. The golden eleven model now has an all-period proof

For rooted cyclic binary words with no cyclic `11`, golden reading into
R_n=Z[phi]/(phi^n-1) is bijective for odd n. For even n it is surjective with
exactly the three-way zero fibre consisting of the zero word and the two
alternating words. The earlier repo statement was bounded by computation.
A complete carry bound and finite transition classification now prove the
all-period result; a separate index argument proves surjectivity.

This is also a coding of rational 2-adic standard-Collatz periodic sources,
not a theorem of convergence for positive integers. Here a period-n source
means a point fixed by n iterations, not necessarily of least period n.
Standard and shortcut
clocks remain different: the inherited -17 orbit has standard period18 and
shortcut period11. The period-five golden eleven is instead Lucas L_5=11.

Period extension has an exact additive lift J_(m,r) obtained by multiplying
by `1+phi^m+...+phi^((r-1)m)`. It is injective, but return projection obeys
`rho J = multiplication by r`, not identity. Smith invariant factors give
the exact size of the return kernel. The chosen word length changes from m
to rm, even though the underlying periodic source is unchanged.

At m=5, the three repeat counts2,3,11 give useful distinct return operators:
2 swaps Paley colours, 3 preserves them and acts as a five-cycle rotation,
and11 kills every projected label. In the last case all eleven lifted
labels remain distinct and form a square-zero ideal. The scalar examples
instantiate residue/nonresidue/zero types; they do not explain the separate
bracket selection of those primes.

### 3. Odd-square brackets have a weighted residue formula

The inherited doubling-stable integers are `{2,3,4,10,11,12}` and the prime
filter gives `{2,3,11}`. For the relaxed multiplier q/(q-1), q>=2, exactly
brackets1 through2q-2 contribute. Their total integer population is

    N_q = (16q^2-24q+8)/3 - 2 T_q/q,
    T_q = sum_{m=0}^{q-1} [-(2m+1)^2]_q.

For an odd prime p this becomes a polynomial plus `2 S_p/p`, where
`S_p=sum a chi_p(a)`. This weighted quadratic-character sum is a marked
Paley-row moment when p=3 mod4. At11, S_11=-11 and N_11=548.
The arithmetic packet proves an exact prime-power lift, including nonunit
and zero-root strata. The integer count is NOT a prime count.

The integer ruler is essential: moving the Paley origin while holding ruler
weights fixed changes the eleven row moments through
`-11,0,-11,0,11,22,11,0,-11,0,-11`.
Transporting both ruler and origin restores covariance; forgetting the ruler
does not. Likewise arbitrary source mass must be pushed through square
fibres, not replaced by uniform residue weights.

### 4. A concrete golden-to-341 connection, with a source warning

The CRT supplement supplies the exact ring decomposition

    R_15 = F_4 x Z/341 = F_4 x F_11 x F_31,
    phi -> (alpha,81) -> (alpha,4,19).

This makes the old joint-clock subgroup `<4,81>=QR_11 x QR_31` a literal
quotient of the golden arithmetic, not merely another group of order75.
Scalar4 is a unit clock only AFTER discarding the F_4 component, not on
all of R_15. The quotient has four-element fibres. Its four-cell kernel is
exactly the lifted period-dividing-three source set `{0,1,2,4}`: the ordinary positive
cycle's phases and zero disappear under this quotient. This does not identify
other fibres with ordinary rational-additive translates of those sources.

Crucially, extending the eleven period-five sources to period15 does not
put them in that unit subgroup. Their images are the conductor11 axis,
`124a mod341`, retaining the return phase3. A units-only representation
therefore drops every one of those sources. For the bracket statistic the
two nonunit axes contribute an exact missing correction36; the11-axis
alone contributes12. Even mass normalization does not repair the wrong
source distribution: its ceiling-tax mean is186 versus174 for ambient
uniform roots.

## What unifies the results

The useful common mechanism is not the numeral11. It is the discipline of
transporting a source through an operation while retaining the coordinate
that its next observer needs:

- score potential bounds actual complement clocks, but not each raw cycle-debt change;
- golden ideal lifts preserve source labels, while return projections can annihilate them;
- bracket counts need height, residue multiplicities and the integer lift;
- the341 bridge needs conductor strata, not just the attractive unit subgroup;
- square-sum Q25 uses the named11,24 ports and a two-component forest response,
  an inherited attachment mechanism that endpoint-only summaries miss.

The future questions remain typed: larger tournament accessibility and prime
atlases, positive-integer Collatz dynamics, prescribed-endpoint or simultaneous
square-sum edge constraints, and prime distribution in brackets are not closed
by this session. The general finite-valuation potential obstruction was already
proved and is not rebranded as a new result.
