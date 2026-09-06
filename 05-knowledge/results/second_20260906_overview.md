# Three proved continuations and an actual balanced LRC obstruction

**Status: PROVED in the stated scopes, with independent proof audits and
separately scoped FINITE-EXACT controls. Date: 2026-09-06.** LRC(14), the
general Laurent first-return bound, and general higher-jet Smith partitions
remain **OPEN**. No external priority claim is made.

This session continues the [previous creative pass](creative_20260906_overview.md)
by making its interfaces decide new mathematical cases. The anchor is actual
LRC decoder entry, the niche is carried Laurent cancellation, and the wildcard
is projective mixed-jet interpolation. Three positive results survived
independent audits; a balanced hostile identifies the limit of the new LRC
criterion without confusing that limit with unsafe speeds.

## 1. LRC: remove a unit assumption on explicit larger subclasses

The [entry theorem](second_20260906_entry.md) applies to actual primitive
thirteen-speed rows in the physical box with two decoder components and
full equality `W_(Q,3)=V_dec`, where `Q=91^6`. Write them as `tV union gU`,
with coprime primitive scales, `a=|V|<=6`, and `b=|U|=13-a`.

Earlier work closed five unbalanced split types when U contained one.
The new result permits nonunit U when a pair containing its full maximum
has sufficiently small gcd D:

| Split | Sufficient maximum-endpoint gcd D | Change from inherited work |
|---|---:|---|
| 1+12 | `D<=6,240,321,451` | Removes the larger-unit requirement on this subclass |
| 2+11 | `D<=76,388,115` | Recovered inherited endpoint-gcd result |
| 3+10 | `D<=698,294` | New subclass; crude cofactor substitution gave 175,558 |
| 4+9 | `D<=4,854` | New subclass; crude substitution gave 725 |
| 5+8 | `D<=26` | New subclass; crude substitution gave 2 |

The conclusion is an actual full-row phase with clearance at least `1/14`,
relative to the named inherited proper-component LRC and gcd suppliers.
Strict clearance is proved in the constructive scale branch. There are six
actual equality-entry controls with U lacking both a unit and a coprime pair
at its maximum, including a balanced control for the native criterion.

The new ingredient is an elementary **sharp real cofactor inequality**.
For any n-vertex tree with positive edge coefficients of sum at most C,

```text
min rooted cofactor <= C^(n-1) (n-2)^(n-2)/(n-1)^(n-1),   n>=3.
```

Averaging only over leaf roots retains edge-cut leaf proportions and gives
an even sharper topology-dependent bound. The star attains the uniform
real-cofactor constant. Division by integer content transfers the inequality
to primitive kernels; sharpness is not claimed for distinct primitive
integer labels. The resulting bound controls the smaller minimum, while
the exact signed coefficient box preserves the full larger maximum in the
physical phase interval.

The [independent audit](second_20260906_entry_audit.md) recomputes maximal
minors by Bareiss elimination, enumerates trees by connected edge subsets,
rebuilds the six actual entries, and selects their phases independently.

## 2. Laurent: close the six-first-channel endpoint-33 family

For every integer `g>=17` with `gcd(g,33)=1` and every three nonzero complex
coefficients, let

```text
f(z)=alpha*z^(-33)+beta*z^(2g-33)+gamma*z^(3g-33).
```

The [new theorem](second_20260906_laurent.md) proves that the first nonzero
constant-term moment is **exactly g or 2g**. Both possibilities occur for
every allowed g. Thus the first detection is strictly below the support
width `3g` even when the opposite endpoint grows without bound.

There are six complete first-return channels and twelve doubled channels.
The first phase polynomial has five simple negative roots by inherited
THM-4436, with its precise parameters checked. The actual doubled response,
including its inverse-phase lower carry, becomes a polynomial in the
five-dimensional first-root quotient. Every characteristic coefficient of
its multiplication operator has strictly positive coefficients after writing
`y=x-1=g-17>=0`: **all 155 coefficients** are frozen exactly. Consequently
every real first-root response is strictly negative and hence nonzero.

The [independent audit](second_20260906_laurent_audit.md) starts from complete
literal count fibres, computes characteristic coefficients by principal
minors, and reconstructs the polynomial identities by exact interpolation
with independently proved degree bounds. It matches all155 coefficients.
The producer uses symbolic polynomial arithmetic and trace identities;
the audit uses only Python's standard library.

This proof also explains the failed entry through the previous derivative
cone: an actual normalized endpoint15 row has an exact separating functional
against its entire cone. A positive square-phase multiplier guess separately
fails at x=10000. Neither failure refutes noncancellation. Retaining the
actual lower-carry coefficient makes the successful quotient certificate
possible. Endpoint39 and a uniform all-channel theorem remain open.

## 3. Smith factors: classify the full mixed (2,2,1) projective observer

The [complete classification](second_20260906_smith.md) gives all five
Smith exponents at every prime for the lattice of homogeneous binary forms
of degree four, observed at three pairwise nonparallel primitive integer
directions with two value-and-derivative observations and one value observation.

The weighted bracket metric determines the partition except in one
dyadic adjacent-depth shape: the doubled-pair bracket has valuation `e+1`,
and both simple/doubled brackets have valuation e. At outer depth `e>=2`, a single
intrinsic unit bit epsilon determines the complete list:

```text
(0,0,e+2,3e+2-epsilon,4e+epsilon).
```

The bit comes from a projective bracket ratio modulo four, and is invariant
under reference changes, primitive unit rescalings, lawful projective maps,
and exchange of the doubled directions at exactly the active depth. At
shallower depth the naive bit can change, but it has no effect on the
Smith partition. Infinity and the three-class dyadic chart exception are
included.

This turns an inherited counterexample into a classification and an actual
full-kernel consequence: the two equal-metric branches differ in kernel
size by a factor of two precisely at integer precisions
`3e+2<=N<=4e`, an interval of e-1 levels. This concerns the complete kernel,
not an isolated pair of factors. The [root audit](second_20260906_root_audit.md)
independently verifies full ordinary integer Smith forms and the precision
windows. General multiplicities are not classified here.

## 4. Balanced actual entry does not force the current phase inequalities

The [balanced obstruction](second_20260906_decoder.md) gives two actual
6+7 equality entries satisfying the physical box and every retained
large-subset gcd profile, including complement words. Their smaller primitive
minimum is `185,370,716,505`, above the sufficient balanced threshold
`60,843,134,147`. All six maximum-endpoint pairs fail the new native phase
inequality despite paying its relation-coefficient gate. Both simple
component-grid comparisons also fail.

The mixed-parity entry has an explicit safe phase `11/23` with clearance
`3/23`. All denominators2 through22 fail weak safety, and every odd
sixteenth phase has clearance1/16. Its actual decoder, full mixed-support
exclusion, box, and seven-subset gcds are independently reconstructed in
the [root audit](second_20260906_root_audit.md). Thus retaining more of the
same gcd profile cannot by itself prove the failed sufficient inequalities.
A new phase argument could still certify these rows and the residual class.

## 5. Connection board and remaining decisive tests

| Live concept | Source -> target and preserved predicate | Lost information / retained sidecar | Effect on the other lanes and next test |
|---|---|---|---|
| Rooted cofactors | Weighted tree -> primitive component minimum; preserves kernel equations | Integer content can lower the primitive minimum; retain leaf-cut proportions | LRC cutoffs improve; the balanced witness prevents claiming automatic closure. Search distinct primitive shapes, not raw cofactor equality alone. |
| Signed crossing radius | Short-relation exclusion -> scale bound -> protected phase grid | A chosen pair forgets the full maximum; retain that maximum and coprime physical scales | Connects the previous walk interface to actual safety; all-six-pair failure asks for a different balanced phase object. |
| Actual coefficient normalization | Complete carried fibre -> first-root quotient; preserves every cancellation locus | Quotient forgets off-root values; retain the phase monomial and lower carry | The cone failure motivates the actual operator, as endpoint gcds must retain actual physical scales. Test endpoint39 on its complete fibre. |
| Quotient characteristic coefficients | Actual multiplication operator -> signed values at every real first root | Coefficients alone do not ensure real eigenvalues; retain the proved real-root supplier | Turns an infinite Laurent family into a bounded identity certificate. Seek an all-channel positivity mechanism, not a longer parameter scan. |
| Residual determinantal ideals | Full mixed observer -> two unit factors plus three residual factors | Determinant alone loses intermediate ideals; retain the first ideal and inverse denominators | A three-factor residual makes full classification possible; larger jets require additional ideals before claiming full kernels. |
| Active unit residues | Projective brackets -> intrinsic dyadic bit -> complete kernel precision | A reference can alter shallow raw residues; retain the active depth threshold | Like the Laurent large-x hostile, a coordinate matters only at its actual scale. Determine which residue and ideal survive for the next mixed multiplicity. |

The proofs and exact failures are the current truth sources. The
[mistakes ledger](../../01-canon/MISTAKES.md) records the failed implications
and their repairs; the prior session remains provenance for inherited
interfaces. No tournament representation was forced where no intrinsic
binary observable was needed. No new meta-pattern is promoted from this
single continuation.

## 6. Frozen evidence and scope

| Bundle | Producer exact checks | Independent path |
|---|---|---|
| LRC cofactor/entry | 134,477 weighted trees; six actual entries; 884,724 gates | 26,800 trees from literal minors; full six entries; 207,160 gates |
| Endpoint33 | All155 symbolic coefficients; literal x1,...,24 controls | Literal fibres x1,...,51; principal minors and degree-bounded interpolation; 1,438 gates |
| Complete mixed Smith | 8,732 full rows; 648 projective/reference controls; 215,537 gates | 1,059 ordinary integer matrices; 4,236 prime partitions; full precision windows e2,...,15 |
| Balanced decoder | Two complete actual rows and all retained profile filters | Literal graph, rank, both mixed-support orientations, 3,432 seven-subsets, and actual phases |

Each source uses checks that remain active with Python optimization, and
normal/optimized frozen outputs agree. Full reproduction commands are in
the individual notes. The [manifest](second_20260906_manifest.json) hashes
every source, proof, audit and frozen output. The Laurent producer and root
audit use SymPy1.9; all other session verifiers use the standard library.
Finite banks control implementations and hostiles; the universal results
rest on the displayed proofs and, for Laurent, a complete polynomial identity
certificate with a proved degree bound.
