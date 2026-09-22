# Glued affine blueprint: repairs and a stronger discrepancy dichotomy

**2026-09-21. Status: PROVED elementary statements; CITED classical
orbit-sparsity input with explicitly derived corollaries; FINITE-EXACT
controls. Global Collatz convergence remains OPEN.** No novelty claim.
The attached blueprint is preserved as an
[unaudited source](../reference/GLUED-AFFINE-BLUEPRINT-2026-09-21-SOURCE.md).
Its formalization claims must be checked against actual declarations.

## Inheritance and live concepts

The closest proved mechanism is the product/carry identity in
[the discrepancy note, D1](collatz_guards_20260921_discrepancy.md).
Its bounded-strip exclusion used a source-density argument; the elementary
argument below removes the need for that input and strengthens the
conclusion to a one-sided statement. The classical sparsity input then
strengthens it further. The known negative seven-cycle is the canonical
hostile to a unit-gap restriction; finite positive word realization is the
hostile to treating a long finite word as a divergent integer orbit.
The corrected near miss is deleting either the carry or the no-repeat
hypothesis. The least-used sidecar is the spacing of distinct orbit values.

Anchor: actual signed dynamics and the direction of clock discrepancy.
Niche: reciprocal sums of distinct states. Wildcard: compare additive and
multiplicative budgets on the two sign sheets.

| Live object | Invariant / operation | Missing coordinate / cheap test |
|---|---|---|
| Actual halving word | q=2^K/3^L | Source, carry and valuation guards |
| Distinct orbit values | Rearrangement of reciprocal sum | Repetitions; fixed point 1 is hostile |
| Orbit sparsity | Dyadic shell summation | A power saving, not merely density zero |
| Signed parameter | q*n=n0+(b/3)sum q | b changes the budget's sign |
| Inverse fibre | Consecutive gcd=1 | These are siblings, not successive orbit states |
| Formal package | Public root and actual theorem text | A theorem count is not a theorem identifier |

## 1. Audit of the supplied claims

| Source section | Status and repair |
|---|---|
| 1.1 transported operations | **PROVED inherited ring isomorphism** via h(x)=x+1. It is a ring, not a field; the additive identity is -1 and multiplicative identity is 0. |
| 1.2 triangular identity | **PROVED inherited** with index N-2; adding the root star restores N-1 edges. Combinatorial nonnegative-edge interpretations need the stated nonnegative parameter range. |
| 2.1 affine guard | **PROVED** for actual odd-orbit steps. The source reuses U_b for the raw affine map; this differs from the accelerated U_b in our notes. |
| 2.2 unit-gap cycles | The three clocks do **not** all give the positive trivial cycle: they give (-1), (1), and (-5,-7) at b=1. |
| 2.2 rational cancellation | Cancellation is needed for **integer**, not rational, realization. Every positive exponent word gives its rational fixed cycle. For example (1,3) gives 5/7 <-> 11/7 with gap 7 and uncancelled denominator. |
| 2.2 q invariant | **PROVED** as an integrality/minimal-parameter invariant under repetition; no analytic stability claim follows. |
| 2.3 coprime fibres | **PROVED** for adjacent siblings. With seed 1 the fibre is 1,5,21,85,...; 5 and 85 share 5. All map to 1, so this is not a trajectory of independently refreshed prime factors. |
| 3 divisor sign law | **PROVED inherited**, with uniform integer sampling and proper nontrivial divisors. It does not define an orbit measure. |
| 4.1 Hamiltonian count | **PROVED inherited** for n>=3. At n=2 the reversed single edge has one path, not the formula's two. |
| 4.2 digit indexing | With the displayed formula R_1=31, not 331. The sharp pair is R_1,R_16 with gcd31; the fifteen-term theorem survives. |
| 4.2 Fermat-prime clocks | **PROVED inherited** with the original parameter scope and first-hit indices. No prime-infinitude consequence is added. |
| 5.1 matrix action | **PROVED** order-twelve vector action C2 x S3. Projectivization kills {I,-I}, leaving S3; only the C3 subgroup preserves directed cycle arrows. No fibration or monodromy representation is defined here. |
| 5.2 fruit correspondence | **PROVED inherited** after the decimal repair. The sqrt(5) statement concerns all positive triples with varying fraction sum; the fixed sum-4 curve is one member. |
| 5.3 square-value model | The curve equation is an **additional question**, not a consequence of s+t=(x+1)^2. Its nonisogeny statement is correct. |
| 6.1 full descent margin | **PROVED inherited equivalence** for actual counters; the universal margin remains unproved. |
| 6.2 Lean attribution | **REFUTED attribution:** the current CollatzBlueprintAudit root contains no logarithmic discrepancy theorem. Its 32nd listed certificate is the finite `twentyThree_tally_control`, not bounded-strip exclusion. |
| 6.2 divergent D -> +infinity | **REFUTED proposed behavior for actual positive integer orbits**, by the elementary theorem below; a classical input gives D -> -infinity instead. |

The exact inherited repairs route through the
[Catalan note](catalan_elliptic_20260921_catalan.md),
[prime note](arithmetic_seams_20260921_primes.md),
[vector-lift note](arithmetic_seams_20260921_dynamics.md),
[elliptic note](catalan_elliptic_20260921_elliptic.md), and
[formal package README](../../04-computation/lean/CollatzBlueprintAudit/README.md).
The labels "fluid", "manifold", "geodesic", and "monodromy" do not supply
a defined transport law, measure, connection, or fibration. No assertion
about such structures is used as a mathematical dependency.

## 2. The exact product identity

Let n_0,n_1,... be a positive odd orbit of
`U_b(n)=(3n+b)/2^v2(3n+b)`, with b=+1 or -1. Let k_i be the actual
valuation at the step from n_(i-1), and define

```text
K_L=sum_(i=1)^L k_i,
D_L=K_L-L log_2(3),
q_L=2^K_L/3^L=2^D_L,       q_0=1.
```

Here D is a **clock discrepancy**, not the divisor defect in another lane.
Telescoping the actual equations gives the two exact finite identities

```text
n_L q_L = n_0 product_(i=0)^(L-1) (1+b/(3n_i)),       (B1)
n_L q_L = n_0+(b/3)sum_(i=0)^(L-1)q_i.              (B2)
```

All factors are positive for either sign. No stochastic model is used.

## 3. An elementary one-sided theorem from distinct-state capacity

**PROVED without a literature input.** If a positive U_1 orbit is not
eventually periodic, then

```text
liminf_(L->infinity) D_L = -infinity.                 (B3)
```

More precisely, for infinitely many L,

```text
D_L <= -(8/9)log_2 L+O_(n_0)(1).                    (B4)
```

Proof. A repeated state in a deterministic map forces eventual
periodicity. Thus all n_i are distinct. Every n_i for i>=1 is odd and
prime to three, since `2^k_i n_i=3n_(i-1)+1`. List any L-1 such distinct
positive states in increasing order a_1<...<a_(L-1). The jth positive
integer prime to six is at least 3j-2, so

```text
sum_(i=1)^(L-1) 1/n_i <= sum_(j=1)^(L-1) 1/(3j-2)
                         <= 1+(1/3)H_(L-2),         (B5)
```

for L>=2, with H_0=0. Using log(1+t)<=t, H_m<=1+log(max(1,m)), and
`1+1/(3n_0)<=4/3`, (B1) implies, for all L>=1,

```text
n_L q_L <= C n_0 L^(1/9),
C=(4/3)exp(4/9).                                   (B6)
```

Suppose instead that D_L were eventually bounded below by a constant.
Then q_L would be bounded away from zero and (B6) would put every n_L
below a constant times L^(1/9). Among the first N values there would be
N distinct positive integers below O(N^(1/9)), which is impossible.
This proves (B3), including the unboundedly late quantifier.

For (B4), infinitely many n_L satisfy n_L>=L. Otherwise for all sufficiently
large L, n_L<L, and all sufficiently late states among the first N would
be distinct positive odd integers below N, of which there are at most
ceil(N/2); finitely many early exceptions cannot fill the deficit.
At those indices (B6) gives q_L<=C n_0 L^(-8/9), as required.

The proof is a rearrangement bound, not an estimate of the order in which
residues are visited. Repetitions invalidate it: the fixed orbit n_i=1
has q_L=(4/3)^L and D_L tends to +infinity.

**Elementary dichotomy.** A positive U_1 orbit is eventually periodic iff
D_L is bounded below eventually; equivalently iff D_L tends to +infinity.
For the forward direction, a positive cycle satisfies
`2^K/3^r=product(1+1/(3n_i))>1`, so each traversal adds a positive
discrepancy increment. There are only finitely many phases. The converse
is (B3). This does not prove that the only positive cycle is (1), or that
every orbit is eventually periodic.

## 4. Primary-source strengthening: a full limit on divergent orbits

**CITED input.** Garcia--Tal, *A note on the generalized 3n+1 problem*,
Acta Arithmetica 90 (1999), 245--250,
[original PDF](https://matwbn.icm.edu.pl/ksiazki/aa/aa90/aa9033.pdf),
Proposition 1, Lemma 3, equation (6), and Corollary 1.
The full six-page paper was read; pages 248--249 were also rendered and
visually checked. The shortcut map uses d=2,m=3, satisfying m<d^(d/(d-1)).
The choices of residue -1 and +1 give the positive 3n+1 and 3n-1 maps.
Proposition 1 imports Heppner's estimate; that dependency is **CITED**,
not reproved or formalized here.

The load-bearing consequence is stronger than density zero: there are
constants C>0 and beta<1 such that an infinite orbit's distinct value set
O satisfies `#(O intersect [a,a+X)) <= C X^beta log(2X)` uniformly in
positive integers a,X.
Equation (6) supplies this with beta=max(1-delta_1,delta_2). Its application
to an infinite orbit uses the impossibility of equal-time coalescence
between two different states of that orbit.

The paper has a nonessential overstatement before Lemma 3: positive upper
Banach density does not force every interval to have the stated occupancy.
The actual Lemma 3 assumes occupancy of the particular interval and proves
the bound directly, so that sentence is not used. Its block-shift proof
also permits the nonnegative shift zero when the interval starts near one.
For its stronger claim that the colliding points stay inside the chosen
interval, restrict the shifted set B(k,z) to points from that interval;
the same pigeonhole argument applies. For the orbit bound, even a
collision elsewhere in O would already be impossible. These presentation
issues do not change the cited quantitative consequence.

**PROVED corollary of the cited estimate.** For an infinite orbit,

```text
sum_(x in O) 1/x < infinity.                         (B7)
```

Indeed the shell [2^r,2^(r+1)) contributes at most
`C'(r+1)2^(-r(1-beta))`, a summable bound. A finite orbit also has a finite
set sum, but its time-indexed sum repeats values and diverges; no transfer
from set to sequence is allowed there. On an infinite deterministic orbit
the values are distinct, so the time-indexed sum equals the set sum.
Odd-accelerated states are a subsequence of the shortcut orbit.

For b=1, (B7) makes the positive product in (B1) converge to a finite
P_infinity>1. Distinct positive integer values tend to infinity. Therefore

```text
D_L+log_2 n_L -> log_2(n_0 P_infinity),
D_L -> -infinity.                                  (B8)
```

This is an attributed consequence of the classical quantitative bound,
not a new theorem excluding divergent orbits. It improves (B3) from a
liminf to a full limit. The two possible positive U_1 regimes are now
cleanly separated: eventual periodicity gives D_L -> +infinity; an
infinite orbit would give D_L -> -infinity. Neither case is eliminated.

## 5. The negative-parameter sheet has a different conserved budget

For b=-1, (B2) yields immediately, for **every** positive orbit,

```text
sum_(i>=0)q_i <= 3n_0,   q_L -> 0,   D_L -> -infinity.  (B9)
```

This is inherited from the signed generalization in the earlier
discrepancy note; no distinctness or density input is needed.

There is a useful refinement. If the orbit is eventually periodic, its
values are bounded and q_L -> 0, hence n_L q_L -> 0 and
`sum q_i=3n_0`. If the orbit is infinite, apply (B7) to the 3n-1 shortcut
map. Since `0<1/(3n_i)<=1/3` and
`-log(1-t)<=t/(1-t)<=3t/2`, its product in (B1) converges to a positive
P_infinity<1. Consequently

```text
eventually periodic: sum q_i=3n_0,   lim n_L q_L=0;
infinite orbit:      sum q_i<3n_0,   lim n_L q_L=n_0 P_infinity>0. (B10)
```

Thus on the 3n-1 positive sheet it is **exhaustion of the additive budget**,
not the sign of discrepancy, that distinguishes periodicity from a
hypothetical infinite orbit. The three known cycles all exhaust their
budget. This supplies a precise reason to study that sheet: it has a
positive summable clock and a boundary term measuring what remains.
It does not classify all its cycles or prove universal budget exhaustion.

For comparison, on an infinite positive 3n+1 orbit, (B1)--(B2) give
`sum q_i=3n_0(P_infinity-1)<infinity`; on an eventually periodic one this
sum diverges. Sign reflection transports these statements to the matching
negative half-line, retaining the parameter label.

## 6. Reproduction and formal scope

Run the [exact controls](../../04-computation/experiments/glued_xor_20260921_blueprint.py)
normally and with Python -O. They check literal hostiles, the matrix
relations, the formal certificate inventory, actual product/carry
identities on bounded prefixes, the finite reciprocal rearrangement bound,
and the exact geometric budgets of the three known positive 3n-1 cycles.
No floating-point comparison is used to certify an infinite theorem.
The logarithmic limits above follow from the written proofs and cited
input, not the finite prefix data.

No general discrepancy theorem is claimed to be Lean-verified in this
session. The inherited Lean arithmetic remains valid in its stated scope.
The next formalization target is (B1)--(B2), the finite rearrangement
inequality, and a separate interface for the cited sparsity estimate.
Global Collatz descent and global negative-parameter budget exhaustion
remain **OPEN**.
