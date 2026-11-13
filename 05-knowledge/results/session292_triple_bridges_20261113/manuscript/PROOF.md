# From finite observations to infinite certificates

**Status:** PROVED elementary lemmas and explicit applications; CONDITIONAL
physical program; source-inspired, not a proof of RH or Chan–Pak. The linear
realization argument is standard weighted-automaton linear algebra, re-derived
here. No literature novelty is claimed for it or for Lagrange interpolation.

## Inheritance and precise source use

The supplied [Nielsen–Semita manuscript](https://github.com/startigerjln/Nielsen-and-Semita-Riemann-Proof-in-Lean4/blob/8a79f4a20e21afe8cb27de1283d087b79539537b/unified.tex)
organizes its argument around round-trip semantic fidelity, finite approximation
versus global continuation, and phase/holonomy rigidity. It claims RH and
independence from ZFC; neither is imported here. Its geometric argument explicitly
uses a Phase Rigidity Axiom. We take the questions about reconstruction and lost
information as inspiration, not its asserted conclusions. A Lean `sorry` count
is not an axiom/hypothesis discharge, and no Lean build is asserted in this lane.

The least-used repo antecedent is the knot-context construction
[THM-2961, ordered-continuation poset](../../../../01-canon/theorems/THM-2961-ordered-continuation-poset-and-finite-group-context-tax.md):
full continuation equivalence, not present scalar equality, is the operation
congruence. [S291's endpoint atlas](../../session291_cluster_ports_20261113/cluster/PROOF.md)
retains subset/endpoint data; its attachment guarantee is interface-relative.
The [Collatz/RH coordinate study](../../collatz_procgen_20260923_trunk_rh.md)
already distinguishes a genuine coordinate map from a shared symbol `1/2`.
The hostile below targets *unbounded model dimension*, not the validity of
finite-rank continuation certificates.

## 1. Exact symmetric escape, without any RH input

Fix rational `0<delta<1/2` and let `z=s-1/2`. For `R>0` define

    Q_R(z) = (((z-delta)^2+R^2)((z+delta)^2+R^2))/(R^2+delta^2)^2
           = 1 + [2(R^2-delta^2)z^2+z^4]/(R^2+delta^2)^2.

It has real coefficients, is even, equals one at zero, and has exactly the
four roots `z=±delta±iR`. Thus in the s-plane all four roots lie inside the
critical strip but off its central line. On `|z|<=M`,

    |Q_R(z)-1| <= (2(R^2+delta^2)M^2+M^4)/(R^2+delta^2)^2
               <= 2M^2/R^2 + M^4/R^4.

Consequently these exactly reflection/conjugation-symmetric entire polynomials
converge uniformly to 1 on every compact set while possessing off-central
zeros at heights tending to infinity. This proves only that these stated
symmetries and compact convergence do not imply global zero exclusion.
They are **not** zeta functions, and no assertion is made that they obey
zeta's arithmetic, growth, or functional equation with its gamma factor.

The useful transferred question is: what *uniform complexity bound* prevents
an invisible feature from moving beyond every tested scale?

## 2. A finite continuation certificate when dimension is bounded

Let a finite alphabet A act on row states `x in K^d` by fixed matrices `M_a`,
where K is a field. Let `v_1,...,v_t` be output columns. For a word w let M_w
be its ordered matrix product. Define

    V_j = span{M_w v_i : |w|<=j, 1<=i<=t}.

These spaces increase. If `V_j=V_(j+1)`, each M_a preserves V_j, so all later
spaces equal it. Unless V_0 is zero (when every output is zero), its dimension
is at least one. Strict growth therefore occurs at most d-1 times, and
`V_(d-1)=V_infinity`.

**Finite-observer theorem.** Two states x,y have the same outputs for every
word iff they agree for every word of length at most d-1. More economically,
choose a basis of V_infinity while repeatedly closing under the M_a; equality
on this basis is a complete certificate. Proof: `(x-y)` annihilates the basis
iff it annihilates every M_w v_i. Two systems of dimensions d_1,d_2 can be
compared by block direct sum with bound d_1+d_2-1.

**Guard condition.** If outputs also encode legality, the theorem certifies
legality only when the legal alphabet/domain is actually represented by the
fixed finite-dimensional model. A separate finite guard automaton may be
incorporated by a product-state lift. Arbitrary unbounded physical feasibility
predicates cannot simply be appended as if they were linear outputs. The
point is to establish the carrier first, then use finite tests.

**Incoming physical hostile (INHERITED, not rerun here).**
[The core-descent endpoint result, E2](https://github.com/eliottcassidy2000/math/blob/1ece157001637120f4260478719aef48c2e6d9e9/05-knowledge/results/core_descent_endpoints_20261112/PROOFS.md)
compares `even:0` and `even:2`: their complete ordinary seventeen-query
profiles agree, so every word of the declared twelve-command language has
equal H. Both have H=2663 at the base clock length. But their selected-set
slopes are11292 and17025, giving unit-tail counts22585 and34051. The selected
set O is not a mark in that ordinary observer. This is stronger than agreement
on a long finite prefix: *all* old-language observations agree, yet a newly
admitted operation separates the presentations. It is not a counterexample
to finite-rank completeness for the declared language. The required change
is to retain the selected-set/cover data and prove closure for the enlarged
language, rather than just take more ordinary observations.

## 3. Sharp application to the eight S291 exponential modes

The exact base set in [S291 paired-port arithmetic](../../session291_cluster_ports_20261113/ports/PROOF.md)
is

    B=(1,2,3,4,6,8,9,12).

For any two profiles in their rational span, eight consecutive values at
nonnegative integer lengths determine equality at every length. The evaluation
matrix at m_0,...,m_0+7 is a Vandermonde matrix multiplied by the nonzero diagonal
`diag(b^m_0)`, hence invertible. Equivalently the fixed order-eight recurrence
has characteristic polynomial `prod_(b in B)(X-b)`.

Seven samples do not suffice even for two nonnegative integer exponential
profiles. Set

    W=(-216,1188,-2464,2079,-924,594,-264,7),
    P(m)=1188*2^m+2079*4^m+594*8^m+7*12^m,
    N(m)=216+2464*3^m+924*6^m+264*9^m.

Then P(m)=N(m) for m=0,...,6, while P(7)-N(7)=3991680.
The proof is Lagrange interpolation: if
`w_b=1/prod_(c!=b)(b-c)`, comparison of the leading coefficient in interpolation
of X^m gives `sum w_b b^m=0` for m<7 and 1 for m=7. Here W=3991680*w.

At physical chart lengths m>=2 one may instead use

    W'=(-186624,256608,-236544,112266,-22176,8019,-2816,42),

which gives equality at m=2,...,8 and difference3448811520 at m=9.
Splitting signs again gives positive integer profiles. These are **formal
profiles**, not asserted to be individual tournament counts or lawful sums.
There is no conflict with S291's stronger result that all4096 binary port
selections separate at a single m>=3: that is a narrower coefficient set.

This is a concrete alternative to blindly extending sample size: prove a
fixed response alphabet and dimension, then check a basis. If an operation
introduces new modes, Jordan blocks, or guards, enlarge the carrier before
claiming the certificate persists.

### Sparse-observation upgrade: surviving lengths suffice

The consecutive-sample condition is stronger than necessary. Any eight
**distinct real lengths** determine a profile in this eight-mode span.
This is the elementary Chebyshev-system property of exponentials, not a new
general theorem. For d distinct positive bases, a nonzero sum
`f(x)=sum_j c_j b_j^x` has at most d-1 distinct real zeros. Induct on its
number of nonzero terms: divide by the least base to obtain a constant plus
d-1 exponential terms. The derivative has at most d-2 zeros by induction,
whereas d zeros of f would force at least d-1 derivative zeros by Rolle.
The one-term base case has no zero. Therefore d zero observations imply
all coefficients vanish.

**Exact sieve-to-observer bridge.** The eight lengths
`2,4,5,6,8,10,12,14` all survive the first25-prime filter on99*2^m+23.
Consequently one can insist that observations occur at these arithmetically
admissible lengths without losing profile identifiability. The script checks
both survival and the exact nonzero determinant of their evaluation matrix.
More generally any infinite surviving cylinder contains eight distinct
admissible lengths, hence enough observations for this fixed eight-mode
family. This transfers a *sampling theorem*, not primality or physical
coverage. On an unknown family with unbounded mode dimension, the hypothesis
that made finite sampling complete has again been lost.

## 4. Prime-support compression: exact, but retain absolute size

Write `E=exp(gamma)`, `A(S)=prod_(p in S)p/(p-1)`, and t=log n. For integers n>=2
with nonempty prime support S, `n/phi(n)=A(S)`. Nicolas's corrected observable is

    c(n)=F_S(t)=(A(S)-E log t)*sqrt(t).

Direct differentiation gives

    F'_S(t)=[A(S)-E log t-2E]/(2sqrt(t)).

Thus for all real t>0 the unique maximum occurs at
`t_*=exp(A(S)/E-2)`, with value `2E exp(A(S)/(2E)-1)`.
On the ray `t>=log(rad n)`, clamp t_* to that ray's endpoint. This is an
exact continuous upper envelope for all exponent choices on one prime
support, not a statement that the continuous maximum is attained by an integer.

In particular prime support determines n/phi(n) but not c(n): the size t is
the missing coordinate. Even `{2,3,5}` and `{2,3,7}` produce distinct support
factors15/4 and7/2, while changing exponents preserves either support factor
but changes the corrected observable. This relates to the old
[difference-visibility note](../../difference_visibility_norms_20261004.md):
a prime product density is exact for its specified fiber/sampling law, not
a pointwise realization certificate. The analogy to endpoint profiles is a
loss-ledger statement, not a map from primes to tournaments.

## 5. What can actually carry a Chan–Pak proof

The manuscript's finite-counterexample idea must keep its quantifier type.
With T encoded as a finite tournament, H(T)=c is decidable. Coverage has the
form `for every allowed c, there exists T with H(T)=c`; a finite missing-range
computation does not refute it. A fixed c's nonattainment has a universal
quantifier over tournaments. We make no arithmetical-hierarchy hardness claim.

A computable size bound for an attaining witness would turn each target into
a finite decision. More generally an effective, well-founded target+interface
descent, terminating in physically realized seeds, supplies such a bound by
recursion. The finite-observer theorem can verify that a proposed move preserves
all future outputs **after** a finite-dimensional, physically sound guarded
carrier is proved. It cannot supply the move, its integer positivity, or target
coverage. This is the concrete bridge: finite-rank reconstruction can remove
infinitely many continuation checks from a genuine constructive recursion.
