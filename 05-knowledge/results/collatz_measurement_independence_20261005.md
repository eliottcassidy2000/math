# Measurement independence: finite-bank scope, source-dependent deadlines, and retained carry

2026-10-05, original opus session `opus-2026-10-05-S7`; current scope repair
by Codex after incoming commit `e60612700`.

**Current status. PROVED:** the finite-atom-data obstruction F, the separate
moment lower bound F1, and the source-dependent deadline equivalence in
section 3. **CONDITIONAL:** their application when the stated atom data or
deadline is independently certified. **FINITE-EXACT:** the new minimal scope
audit. The original cone and ledger computations retain their declared finite
universes; large int64 censuses are reported rather than independently audited
for overflow here. **NUMERICAL:** the original floating kernel readouts, mass totals,
and entropy comparisons. **OPEN:** independently positive bounds at every
source and universal Collatz. This is not a canon promotion.

**Correction lineage.** The original Proposition F quantified over every
proved identity while changing a fixed deterministic orbit's fate; its proof
did not construct such alternative models. It is replaced below by an exact
statement about *only* fixed finite atom data and normalization. The original
F1 bound was a lower estimate, not an upper bound on the actual measurement.
The converse from an odd-step deadline reversed an upper weight inequality:
21 has one odd step and weight1/6, below the asserted1/3. Coefficient-cone
exit, actual first descent, and ROOT hitting time are distinct. The original
all-positive-integer equivalence using `E_inf` is withdrawn. The failed
statements remain recoverable at `e60612700`; their mechanisms and strongest
survivors are recorded in `01-canon/MISTAKES.md` and this current text.

[Exact scope audit](collatz_measurement_independence_scope_audit_20261005.out)
([script](../../04-computation/experiments/collatz_measurement_independence_scope_audit_20261005.py)).

## 0. What working backward now gives

An independently established positive source floor forces a finite ROOT
deadline. Conversely, an independently established deadline gives a floor
through the source-dependent
[backward compiler](collatz_backward_measurement_compiler_20261005.md). The source coordinate
cannot be erased: sources `(4^(j+1)-1)/3`, j>=1, all have one odd step, while their
weights `2/((j+1)(j+2))` tend to zero.

A finite bank and normalization alone admit a zero atom at any unlisted
target. This prevents that observation scheme from proving positivity there.
It does not prevent a new arithmetic argument, a global dual inequality, or a
family induction from certifying additional sources. The mixed four-child
family in the backward compiler provides exactly such a separately proved
proper domain. Its scope is not universal coverage.

The exact all-source residual can be stated using actual descent: every odd
`n>1` must have a strictly smaller actual iterate. Coefficient cylinders are
useful sufficient regions only after their carry threshold is retained.
Neither uniqueness of a ROOT word nor uniqueness of a finite optimizer proves
that its certificate exists or has the required sign.

## 1. Inheritance, new tasks as they appeared, and the board

Four Codex notes landed during the previous session and were read as input:

- [source refinement floor](collatz_source_refinement_floor_20261005.md):
  SF1 `lambda(rho(n))^2 >= (160/441) W(n)^3` (sharp nonlinear transport through
  the leaf section); SF4 the evaluation-energy criterion
  `p_m >= 1/(c_0^(-1)+E_m)` with `E_m = sum (1/c_(i+1)-1/c_i)`, equivalent to
  atom positivity; SF5 a proper binary inverse tree of monotonically
  decreasing orbits with the source-only floor `W(n)>=w(3b,9b+1)`,
  `b=bitlength(n)`, recognized by a terminating descent. SF5 is a genuine
  non-recycling floor on an explicit infinite family; the family is defined by
  a guarded monotone-descent construction to5 that certifies it.
- [floor transport deadlines](collatz_floor_transport_deadlines_20261005.md):
  F1 `tau<=N` and `W(n)<=2/((N+1)(N+2))`, sharper than T2's `2/(L+2)`; F2 the
  superlevel set `{W>=epsilon}` is finite with largest element `S^B(1)`; F4
  conditional transport through a verified prefix; the hostile `S^j(3)` with
  `W=4/((j+2)(j+3)(j+4))->0` under inverse siblings (the sibling rays).
- [atom polynomial dual](collatz_atom_polynomial_dual_20261005.md) and
  [localized resolvent floor](collatz_localized_resolvent_floor_20261005.md):
  signed selectors whose exact readout increases to the atom, with the stated
  obligation "an independent positive bound for every source" and the
  admission that the numerical examples reuse completed words.
- [weight threshold receipts](collatz_weight_threshold_receipts_20261005.md):
  the finite superlevel compiler, and the Sun counterexample (THM-4026) as the
  hostile for witnesses escaping in ordinary size.

Concurrent commit `e0777a399` (Codex) integrated and audited the S6 note,
logging a mistakes entry (2026-10-05, refinement-shadow integration audit):
the wordwise rise/slope equivalence fails (165 -> 167), the T3 modulus is a
heuristic blanket comparison because aliases inside the residue cancel part
of the tail, T5 must count visits with multiplicity in general, T1 concerns
one-sided progressions with their lower endpoints, finite census proportions
are not densities, and no `tau = O(log n)` is proved. All accepted; the
statements of this note are made within those scopes.

Closest proved mechanism: Codex's criterion `(8)-(10)` and SF4. Canonical
hostile: the deleted geometric atom (every residue class positive, atom zero).
Corrected near miss: reading a bank-based readout as information about targets
outside the bank. Least-used sidecar: the retained actual source and its guarded relation to a grounded bank. Board: **certified bank / forward
orbit avoidance / kernel readout / deadline / rising cone / Beatty depth /
basin minimum.**

## 2. Proposition F: the observation model must be retained

Let `B` be a finite set of atom indices. Suppose the only supplied constraints
are `p_j=a_j>=0` for `j in B`, and `sum_j p_j=1`. Put
`u=1-sum_B a_j>=0`. The admissible class in this proposition is all probability
sequences satisfying exactly those constraints; it is not asserted to be the
class of canonical Collatz measures.

**F — PROVED, finite-atom-data scope.** For any `m not in B`, the possible
values of `p_m` form `[0,u]`. Indeed distribute the residual mass between m
and one other unlisted index. In particular these constraints alone cannot
force a positive target atom. Additional arithmetic, operator, or global
moment constraints may shrink this class and must be audited separately.

The original orbit-based proof cannot justify a stronger statement. With
`C={5}`, source21 avoids C but reaches1, and `R(C)` is not forward closed:
5 belongs to it but its image1 does not. If C is forward closed and contains1,
then the full predecessor set `R(C)` is already the entire ROOT basin, rather
than a finite epistemic bank. Arbitrary potentials `G lambda'` preserve a
weaker flow equation, not all canonical Collatz counter identities. Changing
an unknown deterministic orbit's fate does not preserve the full operator by
assertion. Transport can validly certify sources outside a supplied finite
bank once their connecting relation is independently checked.

Let `h=h_m(j)=4t/(1+t)^2`, `t=2^(m-j)`, and
`q_d(j)=(9h-8)h^d`. Write `A_B=sum_B a_j q_d(j)` and
`v=sup_(j not in B) h_m(j)`.

**F1 — PROVED, two different lower estimates.** Separate enclosures of the
two moments give the valid lower bound

    L_sep = A_B - 8u v^d <= A_(m,d).

If `m not in B`, then `v=1`, every bank summand is nonpositive, and
`L_sep<=-8u`. This says that *this lower bound* fails; it does not say that the
actual readout is negative. For example B={0}, a_0=1/2, m=1, d=1 and residual
mass1/2 at m give `L_sep=-4` but `A_(1,1)=1/2`.

The sharp lower envelope using the joint information is instead

    inf A_(m,d) = A_B + u inf_(j not in B) q_d(j),

and the analogous supremum uses `sup q_d`. These are sharp as infimum and
supremum, whether or not an endpoint is attained. They can be much tighter
than separate moment bounds. They still cannot force a positive atom outside
B, since F supplies a zero-target completion in this observation class.

If `m in B`, the missing index nearest m determines `v<1`, and
`L_sep -> a_m`. The formula `h_m(max(B)+1)` applies to a contiguous initial
bank, not an arbitrary finite bank. B={0,2}, m=0 is the minimal gap witness:
the missing index1 has h=8/9, larger than h_0(3)=32/81.

The original S7 output computes L_sep with float64 bank weights, not a
certified evaluation of A. Its apparent equality with an atom at finite degree
is rounding. For the actual law, there are infinitely many known positive
leaf atoms; for every finite d at least one has `0<h<8/9`, making the exact
actual readout strictly smaller than the target atom. The numerical output is
retained as an experiment, not as an interval certificate.

## 3. Working backward with unique certificates

For the actual strict ROOT word, put `tau=length`, `L=tau-1`,
`K=sum floor((a_i-1)/2)`, and `N=L+K`. The inherited bound is

    W(n) <= 2/((N+1)(N+2)),       tau<=N.

An independently proved `W(n)>=epsilon(n)>0` therefore supplies a computable
upper deadline. The reverse implication cannot reverse this upper inequality.
For n21 the word is(6), tau1, K2, and W1/6. Along the entire one-step ray
`n=(4^(j+1)-1)/3`, j>=1, W tends to zero with tau fixed.

**Correct converse — PROVED, conditional on the deadline.** If n>1 and an
independent argument proves `tau(n)<=T`, the telescoping identity gives

    A=sum a_i <= floor(log2(n(10/3)^T)).

With `e=#{even a_i}>=1`, the exact identity is
`N=(A+tau-e)/2-1`. The
[backward compiler, sections 2-3](collatz_backward_measurement_compiler_20261005.md)
computes an integer upper bound C(n,T) for N and the explicit floor

    eta(n,T) = 2 / ((C+2) binom(C+1,floor((C+1)/2))) > 0.

The source and deadline are both necessary inputs to this proof. ROOT1 has
weight1 separately. Thus an everywhere computable, independently valid
positive rational floor is equivalent to an everywhere computable, independently
valid ROOT deadline. Existence of such a deadline is equivalent to Collatz:
if all sources reach1, literal search is a total algorithm; this conditional
observation is not a present proof of its totality.

A second valid equivalent is that **every odd n>1 has an actual strict
first descent**. Sufficiency is strong induction on n, and necessity follows
from reaching1. A computable bound on first descent gives a ROOT deadline by
recursing on the actual smaller endpoint. One local descent time is not the
whole ROOT time:41 first descends in two shortcut steps,41->62->31, but takes
40 odd steps to reach1.

### Keep coefficient survival separate

For the shortcut map T(n)=(3n+1)/2 on odd n and n/2 on even n, a prefix of
length j with k odd steps has

    T^j(n)=(3^k n+B_j)/2^j,       B_j>=0.

Actual descent requires `(2^j-3^k)n>B_j`. Define `E_inf` using the coefficient
condition `3^(k_j)>2^j` at every nonempty prefix. Every positive integer in
E_inf has no actual first descent. The converse has not been proved: a
coefficient exit alone omits the carry. Even at ROOT, bits10 give coefficient
3/4 and endpoint1. A hypothetical nontrivial positive cycle would have a
minimum with no strict descent, but its whole-cycle coefficient is below1;
such a cycle would not be excluded merely by ruling out positive E_inf.

[THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
explicitly retains the carry threshold and leaves universal coefficient/actual
stopping-time equality OPEN. Its small exceptional representative is useful
for future certificate synthesis, not permission to erase the exception.
The [sibling dimension ladder, Theorem1(a)](collatz_procgen_20260922_sibling_dimension_ladder.md)
proves the dimension `H(log_3 2)` for this coefficient-survival set. That
valid dimension statement and the coefficient-count heuristic do not establish
an integer-avoidance theorem or a ROOT deadline. A finite-looking cone-exit test
must still be connected to actual descent and then to the inductive suffix.

### What the original census actually checks

The original S6 script compares an exact integer first-descent count in its
finite tested universe with a floating predicted count based on coefficient
prefixes. Its tolerance is `0.02*predicted+50`; the printed ratio1.0000 is not
an exact equality check. It scales that prediction by `2^(head_bits-1)-1` and omits1 from the
actual first-descent census. Exact coefficient-cylinder counts in a complete
dyadic interval use `N(A) 2^(head_bits-A)` when A<=head_bits; this explains
the fractional discrepancy in the original prediction column. The retained figures are finite
observations: maximum287 shortcut steps below2^24, and the printed near-unit
ratios. They prove no universal equality of the two predicates or limiting
dimension statement. The finite regression audit here uses Python integers;
it does not rerun or independently certify the large original censuses.

## 4. Owed items

**Shadow theorem, exhaustive audit (FINITE-EXACT; scope: rising cones only).** For all 458 rising words
of length at most 7: every odd `m` below `4*2^(A+1)` follows `w` with exact
valuations iff `m = x_w mod 2^(A+1)` (2,281,464 checks), and every odd `n`
below `6*3^l` has an integral backward chain along `w` iff `n = x_w mod 3^l`
(2,368,206 checks), with the chain positive, returning to `m<n`. Inside a
rising cone there are no exceptions: integrality is exactly the class
condition and positivity is automatic since `2^a u-1>=1`. The class formula
`c_w 2^(-A) mod 3^l` agrees with `c_w (2^A-3^l)^(-1) mod 3^l`. This audit
does not address descents through non-rising words: for a word with
`2^A>3^l` the chain endpoint is below `n` exactly when `(2^A-3^l) n < c_w`,
a finite initial segment of the class, and such descents exist (165 -> 167).
Whether they ever occur outside every rising cone is open; the census of
section 4 counts all smaller ancestors, the cone series counts rising cones.

**Segment ledger, direct audit (FINITE-EXACT).** For the stopping and
single-rise rules with starts below `2^16`, the incoming sum at every odd
target below `4096`, computed by enumerating the recorded predecessors
directly (52,062 and 43,697 of them), equals the through mass plus the exit
mass, and the defect equals `nu(y)` minus the exit mass.

**Lemma T4b, both directions (PROVED).** For `l>=2` a primitive rising word of
length `l` exists iff `{l log_2 3} < log_2(3/2)`; depth 1 is the separate base
case `(1)` (its fractional part equals the threshold). In integer form the
test is `1 + ceil((l-1) log_2 3) <= floor(l log_2 3)`, which covers `l=1` as
well. Necessity was proved before. For
sufficiency take `a_1=1` and, for `j=1,...,l-1`,
`a_(l-j+1) = ceil(j log_2 3) - ceil((j-1) log_2 3)`, so that the suffix of
length `j` has sum exactly `ceil(j log_2 3)` and is not rising, while the
total `1+ceil((l-1) log_2 3) <= floor(l log_2 3)` makes the word rising. The
letters are 1 or 2, and `ceil(j log_2 3)` is the bit length of `3^j`. S3
builds the word at every admissible depth up to 60; the admissible depths
below 60 are `1,2,4,6,7,9,11,12,14,16,18,19,21,23,24,26,28,30,31,33,35,36,38,
40,42,43,45,47,48,50,52,53,55,57,59,60`, a Beatty set. S4 checks that the new
cone classes vanish exactly at the inadmissible depths through 16.

**Cone series and census (finite proportions, not proved densities).** New cone classes per depth through 16:
`1,1,0,1,0,2,8,0,28,0,124,602,0,2498,0,12319` (2,999,301 rising words
visited); the cone density series reaches `0.46726` at depth 16; the sieve
proportion of `D` below `2^28` is `0.468669` (units `0.703004`), constant to five
digits across the last four dyadic blocks. The remaining `0.0014` is carried by
deeper rising cones and, possibly, by non-rising descents outside every cone.

**Mass split below `2^22` (W head 2.850934).** Root 1.000000 (35.1%), leaves
0.656784 (23.0%), descent set 1.016542 (35.7%), basin minima 0.177608 (6.2%).
Basin minima are 415,230 of 1,398,101 units (29.7%) but carry only 6.2% of the
head mass: they are mass-poor because their atoms are sums over leaves above
them, all larger. This gives a numerical head profile for obligation (iii). No smaller
ancestor in the declared graph is a particular induction obstruction, not a
proof that every other induction or global inequality must fail there.

## 5. What remains

1. Prove actual first descent on more source domains, retaining a proper
   decreasing rank and grounded seeds. For each finite X, eventual vanishing
   of the actual no-descent census over odd3<=n<X is the correct bounded
   obligation. Proving it for every X is equivalent to Collatz. A logarithmic
   deadline is not proved or required by this equivalence.
2. Supply independent moment inequalities or global bounded Poisson duals
   that exclude the zero-target completion at residual sources. Proposition F
   restricts the data-only model, not these stronger arithmetic premises.
3. Preserve the domain of every new family. The
   [mixed four-child construction](collatz_inductive_floor_receipts_20261005.md)
   and [unique finite-degree optimizer](collatz_moment_localizer_feasibility_20261005.md)
   give respectively a proper positivity domain and a canonical optimization
   problem. Neither supplies the missing sign at every source.
4. Improve cone-density tail bounds while keeping their residual ordinary
   representatives and carry. Finite proportions and mass splits remain
   separate from limiting densities and universal coverage.

## 6. Reproduction and scope

The new minimal audit passes89 exact Python-integer/Fraction checks, with
normal, optimized, and saved output agreeing:

```text
python -B 04-computation/experiments/collatz_measurement_independence_scope_audit_20261005.py
python -B -O 04-computation/experiments/collatz_measurement_independence_scope_audit_20261005.py
```

The following original large runs are retained as their producer's record:

[Script](../../04-computation/experiments/collatz_measurement_independence_20261005.py),
[output](collatz_measurement_independence_20261005.out),
[JSON](collatz_measurement_independence_20261005.json):

```text
python3 04-computation/experiments/collatz_measurement_independence_20261005.py --cone-depth 16 --sieve-bits 28 --head-bits 22 --json 05-knowledge/results/collatz_measurement_independence_20261005.json
python3 -O 04-computation/experiments/collatz_measurement_independence_20261005.py --cone-depth 10 --sieve-bits 20 --head-bits 18 --audit-len 5
```

The original saved run reports 4,658,591 explicit checks in 9.8 s (numba for the cone DFS, the sieves, the
counters and the first-descent census; exact integers and Fractions for the
audits). Hostiles: the inadmissible Beatty depths must give zero new classes;
the out-of-bank readouts must be negative at every tested degree; the census
and prefix comparison uses the tolerance quoted in section 3, not exact equality. Limits: the
readouts are float64; large trajectory censuses use int64 and have not been
independently overflow-audited in this repair. No finite statement is read as
an unbounded one, and the original bank experiment certifies no new source.
