# Independent audit: palette orbits and weighted response minimality

**Auditor:** `audit_palette_response_minimality`  
**Date:** 2026-10-25  
**Status:** **PASS for the explicitly scoped theorems.** The prose precision fix
requested in §5 has been applied and independently replayed. The coefficient
and integer refinements in Appendices A–B are independently proved here, not
an unreviewed extension silently attributed to the producer.

Audited producer: `session62-palettes/RESULT.md` and `check.py`, against repository
inheritance at `fde0dc637d81397de2b857e61d1e993ba11231bc`. Initial producer hashes
are recorded in `producer.sha256`; final hashes are in `producer.final.sha256`
and post-fix replay outputs are `producer-final.normal.json` and
`producer-final.optimized.json`. The final RESULT and checker SHA256 hashes
are respectively `cb1d3a5d9b599d623e63aeea7412a3ad184470f8ebf7f20c48a99e1580e87211`
and `98a44167e2191065fd85d3295f01d4dabbea444a13bc98494d1779d3b8bdd30b`. No repository edits or pushes were made by this auditor.

## 1. Orbit counts and exact joint composition

Let `q >= 2`, `0 <= k <= q`, and fix the `k` named palette elements pointwise.
An assignment of `w` ordered variables has an orbit determined by (i) which
variables take which named constants, and (ii) the entire equality partition
of the other variables. At most `q-k` unnamed blocks are possible. This proves

```
N(q,k,w) = sum_j binom(w,j) k^(w-j) sum_{b=0}^{q-k} S(j,b).
```

Here `S(0,0)=1` and `0^0=1`. The boundary cases are correct:
`N(q,0,w)=B_{<=q}(w)`, `N(q,q,w)=q^w`, and `N(q,k,0)=1`.
An alternative independent recurrence counts assignments with exactly `r`
unnamed blocks:

```
t(0,r) = [r=0],
t(w+1,r) = (k+r)t(w,r) + t(w,r-1),     0 <= r <= q-k.
```

It gives the same formula without using the producer's canonicalizer.

The dead/free states contribute two copies of the boundary orbit space.
A forced state contributes one extra colour variable. The path ambient row
budget is therefore exactly `2N(w)+N(w+1)`, and retaining one cycle anchor
gives `2N(w+1)+N(w+2)`. These count the whole ambient state carrier, not a
claim that all rows occur in every graph. The numerical example `q=6,w=6`
is correct: 1,282 and 5,863 versus 373,248 and 2,239,488 labelled rows.

**Required symmetry:** every constraint, action and observation must commute
with the same palette group fixing all named constants. A fixed singleton
letter `{0}` does not descend to the unnamed quotient: forced(0) and
forced(1) are in the same orbit, but the letter kills the former and accepts
the latter. Including 0 as a constant repairs this particular action.
The producer correctly states this restriction; the later fully named
response-rank theorem uses a different context bank.

For composition, a union partition is admissible precisely when its two
restrictions belong to the complete input relations, with fixed constants
kept distinct and no more than `q-k` unnamed blocks. Any representative of
that union partition can be matched by separately renaming the two witnesses,
fixing all marked colours. They then agree on the overlap. This proves
sufficiency; restriction proves necessity. Projection is existential restriction
of the full joint relation. Independent canonical representatives do not
supply every possible relative alignment. The two disjoint universal unary
relations have one orbit each, but their join has both equality and inequality.

For graph-fragment applications the inherited scope remains essential:
all shared vertices must be retained as ports on both pieces, every additional
cross-constraint must be declared, and forgotten ports cannot later be queried.
An abstract relation join does not certify these geometric/decomposition facts.

## 2. Complete correlations, counts, and provenance

The three-port induced unit graph is inherited exactly from
`05-knowledge/results/spindle_tournament_interfaces_20261011g.md`, §3.
The complete boundary relation is all 24 nonconstant triples of three colours,
each with two triangle extensions. It has four equality-partition orbits.
Every pairwise marginal is full, so equality partitions here mean a **complete
set of full joint partitions**, not independent pairwise equality information.
The six-point norm computation and the nonconstant-triple criterion are sound.
The audit independently counts triangle permutations rather than replaying all
`3^6` full assignments.

Credit is adequate and explicit:

* `boundary_palette_compression_20261011.md`, §1 supplies complete boundary
  equality-partition relations and same-boundary gluing/projection.
* `spindle_tournament_interfaces_20261011g.md`, §4 supplies full cross-interface
  amalgamation with distinct port sets and declared extra edges.
* `ports_regular_actions_hn_20261025/hn/README.md`, §§2–5 supplies the Boolean
  automaton, anchored closure and named endpoint-count recurrence.
* Its `RELATIONAL_LIFT.md`, §§1–3 supplies the common-exterior relation and
  named separator row budgets.
* `boundary_context_weights_20261011d.md`, §2 supplies conditional weights and
  the forgotten-new-colour multiplicity. The present Boolean orbit budgets
  do not remove that multiplicity for weighted computations.
* `session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md`, §3 supplies
  the fixed-pair-label `F2^3` frame/parity comparison and `2*b1` loss.

All those paths are relative to `05-knowledge/results/` except where already
written fully. Full permutation holonomy versus sign is sound. A full frame
requires identity holonomy; a single-colour section requires only a common
fixed point. The `q=4` 3-cycle example correctly separates them. The pair-frame
claim is kept within constant coefficients and fixed pair labels.

Actual word products and literal extension counts have nonnegative integer
entries. Rational linear combinations of response matrices need not be literal
word products or nonnegative. Boolean support is recovered by positivity of
integer counts, not by testing a modular or signed sum for nonzero. Already
one full-list edge at `q=3` has six colourings, becoming zero modulo 2 or 3.

## 3. Three different linear ranks over Q

Write `A=J-I`, `r_j` for its jth row, and `M_i=e_i r_i`.
For `i != j`, `M_i M_j=e_i r_j`; for `i=j`, the needed product is `M_i`
itself, not its zero square. Because `det A=(-1)^(q-1)(q-1) != 0` over Q,
the `q^2` matrices `e_i r_j` are a basis of `Mat_q(Q)`. All are genuine
singleton-list products of length one or two. This verifies the full algebra
claim, including `q=2`.

1. **Processed nonempty prefix:** the q-vector of endpoint counts is sufficient.
   Singleton prefixes yield `e_i`, and singleton continuation queries give A.
   Hence the response rank and minimum linear dimension are q.
2. **Whole-word series with `F(epsilon)=1`:** the displayed `(s,v)` model
   gives dimension `q+1`. Its empty/singleton prefix-suffix matrix is
   `H=[1 1^T; 1 A]`. Subtract its first row from every subsequent row to get
   `[1 1^T; 0 -I]`. Thus `det H=(-1)^q`. Every d-dimensional linear model
   factors this matrix through its d-dimensional state space, giving
   `d >= q+1`. This is a complete direct proof; no general classical Hankel
   minimal-realization theorem or literature attribution is needed.
3. **Fully named two-port response:** the q² basis words above span every
   response matrix, and entry queries span the dual. Storing the matrix
   is sufficient and the q²-by-q² basis/query table is nonsingular, giving
   exact dimension q². Named singleton letters and named entry queries are
   indispensable hypotheses for this stated lower bound.

These ranks concern linear realizations over Q, not arbitrary nonlinear
encodings and not the number of deterministic Boolean states. The trace
convention counts proper listed cycles for length at least three, but
`tr(T_epsilon)=q` is not the empty graph count 1.

## 4. Independent and replay evidence

Reproduce the independent test with

```
python3 check_audit.py
python3 -O check_audit.py
```

Both modes currently give identical JSON and **9,866 always-active checks**:

* 520 marked-block recurrence parameter cases (`2<=q<=9`, all k, `0<=w<=9`);
* 60 explicit permutation-group parameter cases (`2<=q<=4`, all k, `0<=w<=4`),
  with complete equality-graph orbit classes and actual state-action budgets;
* 1,568 exhaustive pairs of invariant binary relations for
  `(q,k)=(2,0),(2,1),(2,2),(3,0),(3,1)`, compared through literal witness
  matching over a common variable, then projection;
* 203 named response basis words, directly counted from colour assignments
  for `2<=q<=8`, with rational and prime-characteristic rank controls;
* larger literal Hankel blocks at q=2 (prefix/suffix word length at most two)
  and q=3 (at most one), beyond the decisive singleton minor;
* 4,203 integer row vectors with coordinates in `[-2,2]` and `2<=q<=6`,
  reconstructed by the explicit lattice inverse from Appendix B;
* 16,739 residue vectors used to check the quotient by the identity's diagonal
  row-sum class, plus symmetry, scalar-sidecar, positivity and marginal hostiles.

The test does not import the producer checker. Producer replay is separate:
the initially inspected version gives **156,007**, not the stale handoff
count 156,003. The final precision-fixed version gives **156,010** in both
normal and optimized modes. Both prescribed mutants fail in both modes at
their intended guards before and after the fix; final receipts are in
`mutants.final.json`.
Finite controls certify only the declared universes; the all-q claims have
the direct proofs above and below.

## 5. Precision fix requested and verified

The producer sentence about putting "a weight" on every Boolean state is
ambiguous: an unrestricted `(q+2)`-dimensional weighted realization can encode
the valid `(q+1)`-dimensional model. The intended obstruction must instead
refer to one scalar carried alongside the current deterministic Boolean DFA
state, with homogeneous linear scalar transition/output operations.

Also the inherited vectors `(1,1,0)` and `(2,2,0)` are proportional, so they
show loss of multiplicity in the unweighted DFA but do not by themselves
refute a scalar multiplicity sidecar. A decisive Q-hostile uses q=3 prefixes
`{0,1}` and `{1,2}`, both classified free. Their responses to suffixes `{0}`
and `{2}` are respectively `(1,2)` and `(2,1)`, with determinant -3. Two
histories carried by one state and one homogeneous scalar would give
proportional response rows. The final producer text now states exactly this scoped model, explicitly
permits arbitrary q+2-dimensional weighted automata and nonlinear encodings,
and checks the three nonproportional-hostile facts. The precision fix is
verified; no mathematical blocker remains.

## Appendix A. Exact coefficient-characteristic boundary — PROVED

Let K be any field, q>=2, with counts interpreted in K. Let `W_+` be the
K-linear span of `T_w` for all **nonempty** list words. Then

```
W_+ = {B : every row of B belongs to rowspace_K(A)},   A=J-I.
```

Every `M_L=D_L A` has rows in that rowspace, and every nonempty product
ends on the right with some `M_L`; hence every row of every such product
is a linear combination of rows of A. Conversely the singleton words of
length one or two give every `e_i r_j`, which span precisely the displayed
space. This proves the equality, not merely a finite observation.

If the characteristic does not divide q-1 (including characteristic zero),
A is invertible and `dim W_+=q^2`. If the characteristic p divides q-1,
then `q=1` in K. Solving `Ax=0` gives `x=(sum x)1`, so the kernel is exactly
the one-dimensional all-one line. Therefore `rank A=q-1`,
`rowspace A={v:sum(v)=0}`, and

```
W_+ = {B : B 1=0},             dim W_+=q(q-1).
```

Including the empty word adds I. It is not in `W_+` in the exceptional case
because `I1=1 != 0`, so the exact all-word response-space dimension is

```
q^2                    if char(K) does not divide q-1,
q(q-1)+1               if char(K) divides q-1.
```

This is also the minimum linear interface dimension with named entry queries:
the attainable matrices span the space, the entry queries separate it, and
the space is stable under left multiplication by every `M_L`. The identity
is an allowed initial state. In the exceptional case `q=3,p=2`, the dimension
is 7, not 9.

The full path-counting series still has minimum dimension q+1 over every
field: the same representation is valid and the row-reduced H above has
unit determinant `(-1)^q`. The nonempty-prefix continuation rank is still q
provided the empty continuation is permitted. Its output row `1^T`, together
with `r_i=1^T-e_i^T`, spans all q row vectors in every characteristic. In
exceptional characteristic the singleton suffixes alone span only q-1;
retaining the empty continuation matters for this strengthened claim.

## Appendix B. Exact integer lattice and torsion sidecar — PROVED

Fix q>=2 and m=q-1. Let `L_+` be the **Z-linear span**, allowing signed
integer coefficients, of nonempty-word matrices. Then

```
L_+ = {B in Mat_q(Z) : each row sum is divisible by m}.
```

The proof is explicit. Every row `r_j` of A has sum m. Conversely if an
integer row v has sum `m*t`, choose `a_j=t-v_j`. Then

```
sum_j a_j = q*t-sum_j v_j = t,
sum_j a_j r_j = t*1^T - (a_1,...,a_q) = v.
```

Thus the Z-row lattice of A is exactly the divisible-row-sum lattice.
Every nonempty product has its rows there, while the actual one/two-letter
products `e_i r_j` generate that lattice independently in every row. This
proves the displayed matrix lattice equality.

The surjective homomorphism

```
B -> (row_sum_1(B),...,row_sum_q(B)) mod m
```

has kernel `L_+`. Hence

```
Mat_q(Z)/L_+  is isomorphic to (Z/mZ)^q.
```

Including the empty word gives `L=L_+ + Z I`. The image of I in that quotient
is `(1,...,1)`, generating one diagonal cyclic subgroup. Taking differences
from the last coordinate gives a surjection with exactly that kernel, so

```
Mat_q(Z)/L  is isomorphic to (Z/mZ)^(q-1).
```

Equivalently L consists of the integer matrices whose row sums all have
one common residue modulo m. At q=2, m=1 and both quotients are trivial.
This agrees with the row Smith form `diag(1,...,1,q-1)` but does not import
or need a Smith-form theorem. It exhibits the integer obstruction underlying
the exceptional primes of Appendix A. It is an additive-lattice statement,
not a claim that these signed combinations are literal list-word responses.

## Independent peer review of the new appendices

The producer agent `derive_palette_response_interfaces` separately reviewed
Appendices A–B and returned **PASS**, with no requested edits. It confirmed
the attainable stable rowspace, the named-entry-query lower bound, the role
of empty continuation for prefix rank, the constructive integer coefficients,
and the quotient by the identity's diagonal residue class. It also independently
checked the exceptional field cases `(q,p)=(3,2),(4,3),(5,2),(6,5),(7,2),(7,3)`
and the stated integer row-reconstruction universe. This review is separate
from the auditor's own normal/optimized checks.

## Scope and operational note

No new Hadwiger–Nelson lower bound, uniform exterior width, Euclidean forcing
construction or general graph-colouring algorithm follows. The useful object
is an exact correlated interface and an observer-dependent rank/lattice audit.
The optional characteristic and integer additions have their own complete
proofs and must keep the distinction between signed spans and positive counts.

The read-only startup packet completed. Targeted message inspection via
`python3 agents/processor.py --status --peek --limit 8` was unavailable because
this isolated repository has no `.machine-id`; no identity was registered or
repository file changed. Root owns shared-remote freshness and publication.
