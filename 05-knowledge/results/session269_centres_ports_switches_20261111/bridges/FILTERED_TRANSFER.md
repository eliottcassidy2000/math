# Local certificates for an infinite row budget

**PROVED. Session269, 11 November 2026.** This is a forward-state supplement to
Session267's reduced rational-cone theorem, not a construction of the missing
Jacobian source carrier. It also identifies exactly when a nonnegative transfer
system can be tested by its support. No literature-priority claim is made for
the underlying weighted-graph potential argument.

## 1. The carrier and its certificate

Let K be a field, let `A(t,x)` be an N by N polynomial matrix, and let `b(t,x)`
and `c(t,x)` be polynomial column vectors. Every monomial in A has positive
t-degree. The unique formal solution and scalar output are

    y = b + A y,                 R = c^T y.

For a monomial `t^r x^j` put `w = 2j-r`; the weight of a polynomial is its
maximum monomial weight, with `w(0)=-infinity`. States may have different
budgets. If real numbers `h_1,...,h_N` satisfy

    w(b_v) <= h_v,
    w(A_vs) + h_s <= h_v       for every nonzero A_vs,
    w(c_v) + h_v <= B          for every nonzero c_v,

then **every monomial of R has weight at most B**.

**Proof.** Expand `y=sum_(k>=0) A^k b`. Positive t-degree makes each fixed
coefficient a finite sum. An individual summand is a path: the input monomial
starts with weight at most its state's budget and each transition carries a
state within the next budget. The final observer adds at most `B-h_v`.
Consequently every summand lies in the desired cone. Addition can cancel
monomials but cannot create one outside the cone. This proves the claim over
any field, even when coefficients have signs. There is no inference from a
finite prefix. QED.

For an output `R=sum_m F_m(x)t^m`, taking B=5 gives exactly

    deg F_m <= floor((m+5)/2).

This is Session254/267's homogeneous rigid-debt degree bound. An explicit
source-derived state system plus the displayed local inequalities would pay
all those rigid debts at once. Scalar payment, compatibility, entry and
termination remain separate obligations.

## 2. Positivity makes the criterion exact

Suppose now that coefficients are nonnegative real numbers. Replace each
nonzero monomial in each matrix entry by a directed edge `s -> v` carrying
its weight. Add an entrance edge for each monomial of b and an exit edge for
each monomial of c. Remove states not lying on any entrance-to-exit path.

For a nonempty trimmed carrier, the following are equivalent:

1. Every output monomial has weight at most B.
2. Every entrance-to-exit path has weight at most B.
3. Every directed cycle in the trimmed state graph has nonpositive weight,
   and every simple entrance-to-exit path has weight at most B.
4. The state-budget inequalities of Section1 have a solution.

In (3), simple means no repeated internal state; entrance and exit edges are
retained, including their weights. In the empty case R=0 and the assertion is
vacuous.

**Proof.** Positive coefficients imply that each path's monomial survives in
R, irrespective of other paths with the same exponents. Thus (1) and (2) are
equivalent. A positive cycle lies on an accepting path and can be repeated
indefinitely; it contradicts (2). Conversely, removing a nonpositive cycle
from a path does not decrease its weight, reducing to a simple path. This
proves equivalence with (3). Define h_v to be the largest weight of an
entrance-to-v path. It is finite when all cycles are nonpositive and attained
on a simple path. It majorizes every input weight and satisfies
`h_v>=h_s+w(e)` along every transition. Appending an exit shows
`h_v+w(exit)<=B`. This is (4); Section1 supplies the last implication. QED.

The certificate can be found by a longest-path/difference-constraint
algorithm. It is a finite all-row decision for a **given positive carrier**,
not a bounded test for an unknown carrier. Keeping different state budgets is
essential: a legal cycle may contain a positive edge paid by a negative edge.

## 3. A nontrivial two-state example

Consider

    y_0 = 1 + t y_1,
    y_1 = t x y_0.

State budgets `(h_0,h_1)=(0,1)` are valid, and the output y_0 has B=0.
The positive edge has weight1 and the returning edge has weight-1. Directly,

    y_0 = 1/(1-t^2 x),      y_1 = tx/(1-t^2 x).

A method demanding that every transition have nonpositive weight would reject
this legal carrier. The retained state offset repairs it. Replacing `tx` by
`tx^2` creates a positive cycle: the output becomes `1/(1-t^2 x^2)` and
`2deg F_(2k)-2k=2k` is unbounded. A large initial delay can hide this failure
past any prescribed finite row.

For signed carriers necessity fails. Set

    A = diag(tx,tx),   b=(1,1)^T,   c=(1,-1)^T.

Both reachable cycles have positive weight, but R=0. A support graph alone
cannot certify the true reduced rational denominator. One must either prove
an observable positive realization, retain exact cancellation, or use only
the sufficient direction of Section1. The same caution applies to complex
Stein classes and Hamiltonian source corrections.

## 4. What transfers between this session's three problems

**Planar Jacobian.** Source rows -> a forward polynomial state recurrence ->
weighted transitions -> state budgets. The map preserves actual coefficient
support only after the recurrence has been derived from the source equations.
The theorem supplies a finite certificate to aim for; it does not assert that
every source has such a finite recurrence. Signs destroy necessity but not
sufficiency. Existing Session267 proves the output's reduced-denominator
criterion; this note gives a compositional local certificate before elimination.

**Hadwiger–Nelson.** A fixed finite geometric gadget strip has a nonnegative
transfer matrix indexed by boundary color assignments. Summing or forgetting
interior colors is literal path composition. A state disappears only when
there is no accepted coloring path, not because its positive weight becomes
small. Section2 applies to extra additive statistics on those paths, but the
Euclidean coordinates and injectivity must be retained separately. A transfer
matrix for an arbitrary graph does not supply a plane unit-distance embedding.

**Chan–Pak.** Endpoint-labelled path-cover states and the interacting cell
matrix of Session268 have nonnegative counts. They admit support/cost
potentials when the geometry of their chosen construction is a finite-state
path system. However, a prescribed **integer H value** is an arithmetic landing
condition, not support nonemptiness or a degree bound. A support certificate
therefore cannot prove the odd-spectrum conjecture. The new arithmetic lane
retains the low binary digits and a physical constructor instead.

The transferable operation is a state-dependent budget under composition.
The predicates remain distinct: source support, coloring support, and exact
integer landing. The critical sidecar is correspondingly source provenance,
Euclidean incidence, or integer height beyond a residue.

## Sources and inheritance

- Session254, `jacobian/ALL_ROW_DEPTH.md`: the all-row homogeneous source flag.
- Session267, `planar/EXTENSION.md`: reduced rational support and delayed-pole
  controls; `planar/PROOF.md`: Gaussian rotation need not preserve the flag.
- Session268, `automata/PROOF.md`: actual positive transfer matrix on rooted
  interacting tournament cells, with scalar dynamics distinguished from state.
- `00-navigation/META-PATTERNS.md`, “Normalize repeated-response towers before
  scalarizing” and “Retain the product and shared labels before taking scalar
  totals”: used here without adding a new universal meta-pattern claim.
