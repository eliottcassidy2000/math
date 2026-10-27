# Independent audit: bounded-arity spectral parse and Paley exception budgets

Session75, 2026-10-27. Auditor owns this directory only; no shared repository
edits. Repository source pin: `96e35642110104cc547516a63f37eebcf930f707`.

**Verdict: PASS for the all-shape trace theorem, flat-target dual bound,
rank budget, fixed-size native-word arc tariff, finite-field Paley Gram
identity, C3-tower formula, optional positive-weight bounds, and finite exact
join/C3 arc distances P3=0, P7=3, P11=11.** Two minor
quantifier/coordinate clarifications were requested and incorporated; neither
changes a finite inequality. The independent finite-distance supplementary
review is recorded in Section 9.
No historical novelty claim is made or audited.

## Scope and inheritance

Reviewed producer source: `/tmp/session75-spectral/RESULT.md`, with the
reviewed byte hash to be recorded in `REVIEW_BINDING.json` at close-out.
The primary checks in `audit_checks.py` neither import nor execute producer
code. They depend only on Python's standard library.

The closest recovered mechanism is THM-1960,
`01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md`:
regular-child spectral splitting, including its explicit irregular-child
failure. Its tower specialization is inherited, not a novelty claim.
The native operation and score potential were recovered from THM-4810,
`01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md`.
Session74's six/15 native-complement distances are not interchangeable with
arc-edit distance to a join/C3 grammar. Targeted source copies are retained
beside this report for review only, with source provenance in the binding.

Anchor: all-shape spectral cap. Niche: independently certified finite-field
and regular-tower controls. Wildcard/hostile: a single native complement has
an actual rank-four correction already at order four.
The active concept board is disjoint modules; entropy-balanced cuts;
flat-target duality; physical-pair XOR costs; weighted/regular spectral data.

## 1. Parse typing, unbalanced trees, and genuine modules

The theorem uses actual labelled leaves and a disjoint substitution tree.
After suppressing unary nodes, every internal node has between two and b
nonempty children, where the **effective** bound satisfies b>=2. A transitive
macro-leaf or unbounded transitive quotient may be expanded into binary joins.
It must not remain an arbitrary undecomposed large leaf in the proof. A claimed
bound b=1 is not compatible with such nontrivial expansion: replace it by
max(b,2). The n=1 case is separate, with trace norm zero.

For n>=2 set T=bn/(b+1). The root has size n>T. Descend to a child larger
than T while possible. Sizes decrease strictly and T>=4/3, so the terminal
large node is internal. Every one of its children has size at most T, and
its largest child has size m>T/b=n/(b+1). This gives

    n/(b+1) < m <= bn/(b+1) < n.

Every subtree leaf set is a genuine module of the *whole* tournament: its
orientation to an outside vertex is set at their first distinct ancestor
children. Thus the selected child is a module, not merely a numerical split.
This handles binary combs, arbitrary unary insertions after suppression,
n<=b, and highly unequal subtree sizes. No original-tree depth bound is used.

For every induced subset, intersect all child leaf sets, delete empty ones,
restrict each quotient, and suppress unary nodes. This reconstructs exactly
the induced tournament with arities still at most b. A restriction need not
lie in the original bank; the proof uses the larger hereditary class of all
quotients of order at most b. The complement of the selected module need not
itself be a module. Its induced-parse membership is the load-bearing fact.

For the module cut, the cross matrix is

    X = [[0, 1_m s^T],[-s 1_m^T,0]],  s in {+1,-1}^{n-m}.

Its two nonzero singular values are both sqrt(m(n-m)). The outside vertices
need not all see the module in the same direction: only each outside vertex's
uniformity on the module is needed. The triangle inequality therefore gives
the claimed recurrence with cross cost 2sqrt(m(n-m)).

## 2. Entropy constant: independent analytic proof

Write H(x)=-x ln x-(1-x)ln(1-x). For 0<x<1/2, the sign of the negative
derivative of sqrt(x(1-x))/H(x) is the sign of

    D(x)=2x(1-x)H'(x)-(1-2x)H(x)
        =(1-x)ln(1-x)-x ln x.

Extend D continuously to x=0. Then D(0)=D(1/2)=0 and

    D''(x)=(2x-1)/(x(1-x)) < 0.

Strict concavity gives D(x)>0 inside the interval. Therefore
2sqrt(x(1-x))/h2(x) decreases strictly on (0,1/2), is symmetric under
x -> 1-x, and is at most its value at alpha=1/(b+1) throughout
[alpha,1-alpha]. The entropy identity

    n log2 n - m log2 m -(n-m)log2(n-m) = n h2(m/n)

closes the induction, including singleton pieces. This proves the stated
C_b n log2 n cap at arbitrary parse depth. The numerical grid in the checker
is only a non-load-bearing sanity control; this derivative proof is the
all-order justification.

The producer's rational C3<13/12 comparison and n>=512 sufficient threshold
also check out. The derivative of (n-1)/(sqrt(n)ln n) is positive when
(n+1)ln n>2(n-1), in particular for n>=512. The finite starting inequality
is an exact integer comparison after squaring positive quantities.

**Quantifier clarification:** the finite bounds allow arbitrary b>=2.
The asymptotic half-of-all-arcs and n/4-O(sqrt(n)log n) conclusions hold
for **fixed b**, or under a separately stated growth control. Taking b=n
admits the target itself as one quotient and defeats an unqualified uniform
asymptotic. This clarification was sent to the producer.

## 3. Flat-target edge and rank costs

For skew tournament matrices A and B, each agreeing unordered pair contributes
2 to the Frobenius inner product and each disagreeing pair contributes -2.
Thus, with m physical unordered disagreements,

    <A,B> = n(n-1)-4m.

If AA^T=nI-J, A has operator norm sqrt(n) and n-1 singular values sqrt(n).
Nuclear/operator duality gives

    m >= [n(n-1)-sqrt(n)||B||_*]/4.

This bound is uniform under all relabellings and global converse of B. It
remains valid if the inner product is negative. Integer ceilings and a maximum
with zero are appropriate. The leading constant 1/4 in n^2 is best possible
for a baseline class containing all transitive tournaments, since a random
order has expected n(n-1)/4 backward pairs.

For E=A-B of rank r, let W=ker(E) intersect 1-perp, of dimension
d>=n-1-r. Choose orthonormal v_1,...,v_d in W. Then u_i=A v_i/sqrt(n)
are orthonormal, and

    sum_i u_i^T B v_i = d sqrt(n).

The matrix sum_i u_i v_i^T has operator norm one, so duality gives
||B||_* >= d sqrt(n). This independently proves the rank bound without
assuming an interlacing claim for a nonsymmetric perturbation.

If E is actually a sum of k rank-at-most-two matrices, rank(E)<=2k.
That is a genuine **additive matrix-rank currency**. It is not an arc cost,
module-occurrence count, or automatic cost of a syntactic operation.
Uniformly oriented disjoint complete bipartite reversals have rank two but
can reverse quadratically many physical pairs. Arbitrary nonuniform
rectangles or arbitrary internal module replacement do not have this cap.

## 4. Native complements and free symmetries

For a legal directed Hamilton path P with signed path matrix L,

    K(tau_P T)=-K(T)+2L.

Against a jointly transported baseline, one move converses that baseline
and XORs exactly the n-1 physical path edges into the disagreement set.
After t fixed-size moves the disagreement support therefore has size at
most t(n-1), even with arbitrary legal paths and repetitions. Relabel both
the baseline and the exception support whenever a free relabelling occurs;
free converse changes the sign of both matrices and no support size.

**Coordinate clarification:** with free relabelling the final baseline is a
relabelled B or B^op, not literally B or B^op in the original fixed labels.
The invariant spectral cap makes this harmless once stated. This clarification
was sent to the producer. No modularity of intermediate states is needed.
Initial exceptions add their initial support size. Size-changing operations
need final propagated physical-pair tariffs and are outside the simple bound.

**Exact hostile:** at n=4, start with the transitive tournament and the path
(0,1,2,3). The correction K(tau_P T)+K(T)=2L has rank four. Hence one native
complement must not be invoiced as one rank-two correction. The rank result
and the t(n-1) physical-pair tariff remain independently valid.

The sharper path nuclear-norm formula was algebraically reviewed as well:
the eigenvalues of the Hermitian gauge of iL are 2cos(j*pi/(n+1)). Summing
absolute values yields 2cot(pi/(2(n+1)))-2 for odd n, and the same expression
with csc for even n. Its asymptotic native nuclear tariff is 8n/pi.

## 5. Finite-field Paley and the regular C3 tower

Here is an independent proof of the character Gram identity. For distinct
a,b in a finite field F_q, z -> (z-a)/(z-b) bijects F_q minus {b} onto
F_q minus {1}. Quadratic multiplicativity gives

    chi(z-a)chi(z-b)=chi((z-a)/(z-b))  for z != b.

Summing gives sum_z chi(z-a)chi(z-b)=-1. The diagonal sum is q-1.
When q=3 mod 4, chi(-1)=-1, so A_xy=chi(y-x) is skew and AA^T=qI-J.
This proof applies over the field of any allowed prime-power order. It does
not apply to the ring Z/qZ for composite q.

The independent exact field control uses F3[x]/(x^3-x+1), a different
irreducible cubic from the producer. Rootlessness over F3 proves
irreducibility. The checker additionally verifies every nonzero multiplication
map is a permutation, the character's multiplication law, skewness,
regularity and all 729 entries of the Gram matrix. The naive construction
using residues modulo 27 is explicitly rejected: it is not even skew.

For the C3 tower on 3^k vertices, contrast vectors supported on each prefix
block of length p have squared singular value 3^(2(k-p)-1), with
multiplicity 2*3^p, p=0,...,k-1. Together with the all-ones null vector they
span the whole space. Consequently

    ||K(R_k)||_* = sum_p 2*3^p * 3^(k-p-1)*sqrt(3)
                = 2k*3^(k-1)*sqrt(3) = (2k/sqrt(3))*3^k.

The checker builds this tournament by the first unequal ternary digit,
not by the producer's recursive matrix construction. At k=1,2,3 it checks
the complete contrast eigenbasis using exact matrix-vector arithmetic and
full rank modulo a prime. Full rank modulo a prime implies real independence;
it is not being used to infer a floating-point spectral match.

For odd k, the flat-target formula gives the producer's exact integer bounds:

| k | n | arc disagreements at least | native moves at least |
|---|---:|---:|---:|
| 1 | 3 | 0 | 0 |
| 3 | 27 | 54 | 3 |
| 5 | 243 | 9,234 | 39 |
| 7 | 2,187 | 988,524 | 453 |

The values are lower bounds, not attainability claims. Both endpoint classes
are regular, including all their relabellings and converses. Their sorted
score vectors, score potential and degree variance agree, so these particular
degree-based endpoint differences do not detect the spectral separation.

## 6. Optional positive vertex weights

For K_a=diag(sqrt(a)) K diag(sqrt(a)) and W=sum a, a module cut has cross
nuclear norm 2sqrt(W_A W_B)<=W_A+W_B. Count-balanced recursive cuts have
depth at most ceil(log_{(b+1)/b} n), and each original vertex weight is
charged at most once at a given depth. This proves the safe W-times-depth
bound. It neither assumes mass-balanced cuts nor substitutes weights for
physical leaf multiplicities in the unweighted theorem.

The weighted pairing identity and operator estimate follow directly:

    <DAD,DBD>=W^2-sum a_v^2-4 sum_disagree a_u a_v,
    ||DAD||op <= ||D||op^2 ||A||op = (max a_v)sqrt(n).

These claims do not address signed source laws, ROOT positivity, or arbitrary
clone merging. The distinction is appropriately explicit in the producer.

## 7. Independent finite verification and reproduction

Run from any directory with Python 3.8+:

    python3 -B audit_checks.py > normal.json
    python3 -O -B audit_checks.py > optimized.json
    cmp normal.json optimized.json

No external package, producer import, repository file, network lookup,
or current-working-directory assumption is used. Every guard uses
`require`/`RuntimeError`, not an optimization-removable assertion.
The script was also copied alone to a fresh temporary directory; normal/-O
outputs matched each other and the retained original output byte-for-byte.

Finite universes and exclusions:

* All 2,154 ordered plane shapes with 2<=n<=7 and b in {2,3,4}; one
  deterministic-seeded quotient orientation per shape, not all possible
  quotient decorations. Every selected cut and exact cross Gram is checked.
* Every nonempty induced restriction for those tested n<=4 parses; selected
  module/complement restrictions at larger n. Quotients are restricted by
  their actual positions, not silently replaced by bank representatives.
* Maximally unbalanced binary combs through n=32, unary root wrappers,
  and n=4 under b=4,5,9,100.
* A 10,000-point floating entropy-ratio grid, explicitly not a proof.
* Exact GF27 and tower spectra through n=27, as described above.
* All 206 legal native moves across all labelled tournaments of orders 2,3,4;
  involution and n-1 support tariff; rank-four native hostile.
* 100 deterministic-seeded legal six-move words at n=6 with interleaved
  free relabelling/converse, checking the jointly transported support tariff.
* 128 deterministic-seeded tower27 relabelling/converse controls against
  Paley27. Observed disagreements lie in [147,200], all above 54.
  Sampling checks the implementation, not uniform optimality or attainment.

The all-order statements come from the proofs, not extrapolation from these
finite controls. All targeted controls passed normally and under `-O`.

## 8. Boundaries and remaining frontier

This audit establishes a spectral **obstruction** and honest exception
currencies. It does not establish reachability by native words, an exact
minimum word length for the regular tower examples, LRC(14), arbitrary
overlapping-port decomposition, a finite-bank classification, or historical
priority. The method's new useful coordinate is the large flat-target inner
product, which survives regular-to-regular comparisons where scores agree.

The next focused interface, if desired, is an explicitly typed grammar of
non-modular repairs with computable final arc/rank tariffs, not a reuse of
untyped module-occurrence counts.

## 9. Supplement: independently recomputed finite join/C3 distances

`independent_distances.py` recomputes the exact optimal **maximum pair
correlation** with each prime Paley target, without using the producer's code
or supplied lower-bound vector during the optimization. For a vertex set S,
let R(S) be the maximum sum of target_ij * baseline_ij over its unordered
pairs. For a binary root, the cross reward is |cross(A,B)|. For a cyclic
ternary root it is |cross(A,B)+cross(B,C)+cross(C,A)|. These absolute values
choose between the two quotient orientations. The recurrence maximizes this
reward plus the already-computed child rewards over every partition.

Completeness follows from the root of a lawful parse. A quotient with at most
three vertices is either a binary join, a cyclic triple, or a transitive
triple that compiles into two binary joins. Conversely each transition creates
a lawful parse on disjoint child sets. Nonempty proper children make the
recurrence terminate. We enumerate unordered parts canonically by requiring
the first part to contain the smallest label, and the second part to contain
the smallest remaining label. Thus all vertex allocations, unbalanced shapes,
and parse depths are covered. Relabellings are already included because the
partitions range over all named subsets; converse is included by the two signs.

Finally, since a disagreeing pair changes its reward from +1 to -1,

    minimum_edits(S) = (binom(|S|,2)-R(S))/2.

This independently gives:

| target | exact distance | all subset optima | binary roots | ternary roots |
|---|---:|---:|---:|---:|
| P3 | 0 | 8 | 6 | 1 |
| P7 | 3 | 128 | 966 | 1,701 |
| P11 | 11 | 2,048 | 86,526 | 611,501 |

Only after all independent optima have been computed does the script read the
optional producer certificate. Every one of the 2,184 subset values matches.
The script separately decodes each attaining witness, checks exact distinct
leaf labels, lawful arity/type/sign, every unordered orientation, and the final
disagreement count. The P7 witness disagrees exactly on {0,2}, {3,6}, {4,5};
the P11 witness's eleven pairs are recorded in `distances.normal.json`.

Reproduce with a copy of the reviewed `certificate.json` from the producer:

    python3 -B independent_distances.py --certificate certificate.json > distances.normal.json
    python3 -O -B independent_distances.py --certificate certificate.json > distances.optimized.json
    cmp distances.normal.json distances.optimized.json
    python3 -B test_certificate_corruption.py certificate.json
    python3 -O -B test_certificate_corruption.py certificate.json

The optimizer also runs without `--certificate`, recomputing the minima
independently. Normal/-O full outputs agree byte-for-byte. Four independently
authored corruption cases are rejected in both modes: changed final minimum,
changed interior subset optimum, illegal quotient orientation, and duplicate
physical leaf. These are rejection controls for this certificate interface,
not a claim that every malicious input or schema deviation was fuzzed.

The producer's `check.py` was separately read and replayed normally and with
`-O`; its full receipts agree. It checks all Bellman lower inequalities using
a restricted-growth partition enumerator and direct flip costs, then checks
attaining parses. This is a complementary verification direction to the
auditor's recomputed optimum. All entries of the producer's `SHA256SUMS`
verified at the audited byte version. The code review also checked its GF27,
clone/module, and radius-two controls. The producer's corruption suite is
distinct from the independently authored suite above.

The producer's final prose additionally clarifies that a compressed DAG of
disjoint template *occurrences* may be expanded, but identification of their
physical vertices is not licensed; and that the alternative rank proof uses
rectangular domain restriction by an isometry, not an ambiguous same-subspace
compression. Both clarifications are correct and have been reviewed.
