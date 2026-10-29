# Independent audit — Session98 growth packet

Source checkout verified at `5862c94555eafea4e9a0638a17f7c661f8e45d37`.
No source-repository changes or publication were performed by this audit.

## Findings first

**PASS, scoped:** the prime-surplus theorem, fixed-degree projective-response
recurrence, and partial-star classification have valid proofs under their stated
observers. Their checkers and the missing-comparison checker pass in both normal
and optimized Python; paired outputs are byte-identical. The missing-comparison proof note and independent recursive checker also pass.
Its overly broad recurrence-exclusion wording was corrected before approval.

**Additional independent hostile:** exactly three prime tournaments have five
vertices. Direct labelled enumeration, full permutation automorphism testing,
and direct module testing give group orders `1,1,5`. Their surplus-one inflation
counts are `5,5,1`, respectively. Thus two singleton prime families of the same
arity and the same scalar root count have different first growth coefficients.
This is an actual tournament counterexample, not an artificial slot action.
The aggregate value 11 independently corroborates the composition note.

## Proof audit

* **Composition:** factoring x^k uses sum d c_d=k; valuation E(x^d)>=d
  justifies weighted truncation, and canonical prime-root decomposition supplies
  the required orbit model. One singleton is essential. The triangular recursion
  consumes only prior ordinary counts plus the new prime-jet data. Feeding
  Burnside-derived jets back through composition is correctly called a round
  trip, not independent enumeration. No arbitrary-attachment closure follows.
* **Projective recurrence:** exact-r continuation ratios are constant precisely
  when response differences are proportional to j, tested by r unit cycles
  versus one j-cycle and r-j units. In W=w*2^(r h1), insertion cancels the old h1
  term and leaves k(p-1)/2+p*C(k,2)+k*s_p+r'*k. The normalized hierarchy is the
  same gauge applied to distinguished-cycle insertion. Minimality is properly
  restricted to individual-continuation projective responses. The benchmark
  does retain non-target terminal masses: even at p=1 all k<=r are generated.
* **Partial stars:** cycle-coloring Burnside, labelled adjacency double counting,
  and exposure-order counting prove the three mass formulas. The stabilizer
  acts on flagged extensions, not on their forgotten images. For d>=2 the
  center is intrinsic in the missing-edge graph; for d=1 an even automorphism
  group has an involution fixing the complete outside tournament, so it swaps
  only the two missing endpoints and forces twins. All category changes are
  explicitly marked in the note.
* **Missing comparisons (code-level):** an even j-cycle has one reversing
  antipodal pair orbit of length j/2, forced absent. The other internal pair
  orbits have length j; inter-cycle orbits have length lcm with multiplicity
  gcd and do not reverse. This yields exactly closed_profile, with forced
  debt sum_even j/2. The Sylow bound follows from a free action on 2^d completions. The one-hole
  identity uses marked reversal and twin-pair collapse; neither step is valid
  after silently forgetting the pair mark.
* **Limits:** the n^2 exponential lower bound contradicts every factorial-power
  upper bound implied by a fixed polynomial-coefficient recurrence. In the
  normalization proof, all affected pair orbits of an odd permutation have
  length >=3; the support split gives the stated remainder radius. The twisted
  SCC convolution is correctly separated from ordinary analytic multiplication.

## Independent computation and hostile cases

`independent_audit.py` imports none of the packet's algorithms. It performs:

1. Literal permutation/pair-orbit Burnside through n=7, reproducing
   `1,1,1,2,4,12,56,456` without a partition-exponent implementation.
2. Complete ternary oriented-graph orbit census through n=4, checking every
   representative's automorphism order satisfies v2(|Aut|)<=missing edges.
   Missing-edge histograms are n=3: `2,3,1,1` and
   n=4: `4,10,12,10,4,1,1`.
3. Complete five-vertex tournament enumeration with direct nontrivial module
   rejection and automorphism groups, yielding the prime hostile above.
4. C3 marked-versus-unmarked completion counts 4 versus 3.
5. One missing comparison on two vertices, which is fixed by a transposition
   and therefore defeats an odd-permutation-only extension of Burnside.
6. Two tagged singleton child types at C3, giving 4 baseline orbits rather
   than 1; this attacks removal of the unique-singleton hypothesis.

Normal and optimized independent receipts agree byte for byte.

## Reproduction

```
python3 independent_audit.py > normal.json
python3 -O independent_audit.py > optimized.json
cmp normal.json optimized.json
```

The audit also reran the original composition, defect, and recurrence scripts,
and a byte-for-byte local copy `partial_runner.py` of the partial script (which
writes a receipt beside itself). Outputs are `*.normal.json` and
`*.optimized.json` here. The recurrence run uses full default bounds: direct
partitions through30, hierarchy through20, and fixed targets40,60,100.

Finite checks support, but do not replace, the displayed all-order proofs.
No priority, polynomial-time, scalar recurrence, minimal automaton, or universal
structural-growth conclusion is inferred from these results.

## Additional independent marked-pair check

Direct enumeration gives (pair-marked, vertex-marked) tournament counts
through n=5: `(0,1),(1,2),(4,4),(16,12),(88,48)`. With the independent
partial census, `2 p(n,1)=e_n+v_(n-1)` holds through n=4.

An explicit unmarked non-twin loop occurs at n=4: with lexicographic pairs
`01,02,03,12,13,23`, use signs `(-1,-1,+1,-1,-1,-1)` and reverse `13`.
The resulting tournament is isomorphic to the original, although vertices
1 and3 are not twins. This verifies the warning that unmarked metagraph
loops cannot be classified using only the marked-reversal fixed points.

## Proof anchors and additional checks

* [Defect grammar](../defects/MISSING_COMPARISON_GRAMMAR.md), sections2–5:
  signed pair orbits; cycle insertion; free Sylow action; marked reversal.
* [Limits](../limits/LIMITS_AND_NORMALIZATION.md), lines46–57: affected pair
  loss; lines66–78: support3 and remainder split; lines84–90: uncancelled
  fourth-order pole at4 and strict larger radius for the remainder.
* [Composition](../composition/RESULT.md), sections1–3.
* [Projective recurrence](../recurrence/PROJECTIVE_RESPONSE.md), sections3–5.
* [Partial-star classification](../partials/PARTIAL_STAR_CLASSIFICATION.md),
  sections2–4.

The limits checker and defect recursive grammar checker were additionally
rerun in normal and optimized modes with byte-identical receipts. The
recursive checker shares only elementary polynomial arithmetic and the
comparison target with the partition checker; its insertion algorithm is
independent, but it should not be described as a wholly dependency-free test.

## Global meromorphic extension — independent approval

[ALL_POLES_CORE_EXPANSION.md](../limits/ALL_POLES_CORE_EXPANSION.md) is also
**PROVED as stated**. The coefficient-polynomial paragraph now correctly
specifies nonnegative integers after the audit correction.

I independently derived the rational core summand before reading this note:
for a core of support m and c cycles, its normalized coefficient is
`(n)_m/z_mu * 2^(E_mu-C(m,2)) * 2^(-(m-c)(n-m))`.
Summing the falling-factorial series gives exactly equation(1). Since each
core pair orbit has length at least3, the loss bound is quadratic in m.
Summing conjugacy-class sizes first yields the rigorous majorant
`constant*m!*R^m*2^(-m(m-1)/3)`, not an unjustified exchange of infinite
rational sums. This converges for every fixed R. Separating the finite bank
with m-c=2k gives the unique maximal core3^k and its exact nonzero leading
amplitude. Thus the pole at4^k has order3k+1, and infinitely many finite
singularities exclude D-finiteness of this normalized analytic germ.
The explicit H2 bank and finite coefficient extraction also pass.

A code audit caught removable assertions in the original limits checker.
They have now been replaced by explicit RuntimeError checks. Both corrected
limits checkers were copied into this audit directory and rerun normally and
with `-O`; each output pair is byte-identical. Optimized verification now has
active checks. The first pre-correction run had rewritten the limits receipt
with its identical deterministic content; subsequent runs use audit-local
copies to preserve other lanes' artifacts.

All author-claimed novelty limitations remain appropriate. In particular,
the global pole proof is a new analytic deduction within this session from
an explicitly inherited core formula, not evidence of literature priority.

## Bounded extensions: observable quotients and pair-color palettes

### Universal-palette observation quotients

[QUOTIENT_EXTENSION.md](../composition/QUOTIENT_EXTENSION.md) passes independent
proof review. Formal polynomial equality and equality on the finite ordinary
palette grid are distinct observers. Monic one-variable reduction followed by
repeated interpolation proves the exact grid vanishing ideal; the resulting
rank-minimal claim is properly restricted to linear maps of the ambient jet
space. Projection to lower horizon is well defined and surjective because
observation commutes with discarding higher weighted degrees.

`check_observation_independent.py` constructs evaluation matrices with fresh
ordinary polynomial convolution, without importing the observation algorithm.
Modular elimination at prime1000000007 on forty deterministic sampled palettes
at every horizon gives lower bounds matching all reported formal and ordinary
ranks through12. The original exact kernel identities give matching upper
bounds. Independent normal/-O receipts agree, as do reruns of the original
quotient and reachable-witness checkers.

A receipt scaling inconsistency was identified and repaired: the integer
inventory relation has delta24(z3-z1^3), not the unnormalized delta3 previously
printed beside it. The repaired code rechecks the integer relation, its zero
grid observation, and the actual next-horizon difference−18. Its full pinned
census copy is byte-identical to the inherited t12.txt. Thus the family witness
is valid conditional on the explicitly inherited audited census; this session
does not claim to regenerate those prime representatives. Both selected
families have eleven distinct roots, so no hidden multiplicity is needed.

### Pair-color involution pole classification

[PALETTE_POLE_CLASSIFICATION.md](../limits/PALETTE_POLE_CLASSIFICATION.md) passes
independent proof review. In particular, the earlier stronger odd-core bound
cannot be reused: two-cycles fix unordered pairs. The note correctly permits
at most m/2 such pairs, giving total pair orbits<=m²/4 and a Gaussian majorant
`m!R^m q^(-m(m-2)/4)` with ratio tending to zero. Assumptions q>=2 and the
admissible parity restriction on t are explicit.

For t>0, a=k and m<=2k has unique equality core2^k, with leading amplitude
`(2k)!/(2^k k!)*t^k*q^(k²)`. For t=0 only odd cores remain, with polesq^(2k)
and orders3k+1. The no-cancellation and meromorphic continuation arguments
are rigorous. These are finite-color enumeration theorems, not claims about
continuous interpolation or sufficient metagraph states.

The palette checker was copied here and rerun normal/-O with identical active
check receipts. Independent n2 and n3 checks give `(q+t)/2` and
`q³/6+qt/2+q/3`; all admissible q2..4 rows match. The earlier independent ternary
graph census gives42 classes at n4, matching the q3,t1 specialization.
Both bounded extensions are approved; no remaining blockers were identified.

### Synthesis response lemma

[GROWTH_THEORY.md](../GROWTH_THEORY.md) correctly identifies the removable
response space as `span(C-C)^perp`: the ratio of positive response multipliers
is constant exactly when every difference of feasible continuations is
annihilated. A conserved grade provides a removable subspace, but equality
requires the separately stated spanning guard. The four grades are explicitly
not promoted to an unproved universal isomorphism or metagraph carrier.
For absolute clarity, the synthesis now specifies nonzero weight when deriving
observational equivalence of its weighted contributions; alternatively define
response equivalence before multiplication by w. This does not affect any
positive-weight enumeration application here.

## Portable publication subset

Publish only these five audit files:

1. `AUDIT.md`
2. `independent_audit.py`
3. `normal.json` — independent literal graph/action audit receipt
4. `check_observation_independent.py`
5. `observation.normal.json` — independent modular evaluation-rank receipt

Both independent scripts use only Python's standard library, import no other
packet scripts, and need no local input files. The second script embeds the
claimed rank pairs as explicit regression targets, then constructs its own
evaluation matrices. Its finite-field ranks certify lower bounds, not the
symbolic upper bounds; those are supplied by the separately audited quotient
checker. Both scripts run unchanged from any working directory and write only
to stdout. Their normal and optimized outputs were rerun and are identical.

The copied `*_runner.py` files and duplicated integration-check outputs are
local audit execution aids, not publication dependencies. Do not publish them.
The `reviewed.sha256` files are local provenance records, not required assets.
