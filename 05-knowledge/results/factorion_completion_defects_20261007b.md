# Starting with a completed dynamics and retaining its exact defects

**PROVED:** the all-base finite-core reduction, cycle-preserving histogram
carrier, global rank, and completion-defect identities. **FINITE-EXACT:**
the complete base-6 and base-10 cycle portraits and independent controls.
**OPEN:** a Collatz all-source rank or authenticated absorbing core.

The digit-factorial examples give a useful fully solved comparison model.
One can prove universal arrival at a finite terminal portrait, then modify
one edge of every unwanted cycle to make every integer reach 1. Keeping
the changed edges visible shows exactly what must be repaired to transfer
that completed proof back to the original dynamics.

## 1. Inheritance and the target of the analogy

The closest proved interface is [finite grounded receipt networks,
sections 2–3](paper_local_global_grounding_20261007.md), which retains a
source-to-ROOT boundary rather than just cycles of compatible receipts.
The [child-closure synthesis](collatz_child_closure_synthesis_20261007.md)
adds guarded paid exits and complete generated families, while leaving
arbitrary-source membership open. The prior [Poisson patch
theorem](collatz_poisson_patch_20261005.md) already requires the transported
forcing when eliminating interior states; we do not claim that principle
as new here.

Anchor: explicit proof obligations left by an assumed completed model.
Niche: classify every digit-factorial orbit using a small future-sufficient
carrier. Wildcard: compare ordered run clocks with unordered digit data.
The concept board is **source / future carrier / finite core / rank /
terminal cycle / marked repair**. The canonical hostile is the decimal
fixed point 145: a modified proof that it reaches 1 cannot be exported
through its unchanged self-loop. The least-used sidecar is the set of
artificial edges, with their exact residual charge.

The listed decimal fixed points are the classical factorions
[OEIS A014080](https://oeis.org/A014080); periodic decimal points are
[A188284](https://oeis.org/A188284). The base-6 list is recorded in
[A260651](https://oeis.org/A260651). The enumeration below is an independent
certificate for the toy-model transfer, not a new priority claim about
those known lists.

## 2. A complete finite reduction for every base

For a base `b>=2`, let `F_b(n)` be the sum of the factorials of the
canonical base-b digits of a **positive** integer `n`; `0!=1` is retained.
Set `M=(b-1)!`, and let `D>=2` be the least integer such that

\[
b^{D-1}>DM.                                      \tag{1}
\]

For every `d>=D`, `b^(d-1)>dM`, since the ratio of the left/right sides
increases by `bd/(d+1)>1` at each increment. A d-digit source therefore
satisfies

\[
F_b(n)\le dM<b^{d-1}\le n.
\]

In fact its digit count strictly falls. After at most
`max(0,digits_b(n)-(D-1))` such steps the source has at most `D-1` digits.
Its next image belongs to

\[
\mathcal C_b=\{1,\ldots,(D-1)M\},                \tag{2}
\]

which is forward invariant: `(D-1)M<b^(D-1)`, so every integer in that
interval again has at most `D-1` digits. Every positive orbit thus reaches
a cycle in an explicit finite graph. This is the all-source bridge that
a finite cycle census alone would lack.

## 3. Exact compression: the future does not need the digit order

Let `H(n)=(c_0,...,c_(b-1))` be the digit-count vector and

\[
E(c)=\sum_{j=0}^{b-1}c_j j!,\qquad T(c)=H(E(c)).
\]

Then

\[
F_b=E\circ H,\qquad T=H\circ E,\qquad
H\circ F_b=T\circ H.                             \tag{3}
\]

The finite histogram carrier consists of all vectors with total count
between 1 and `D-1`, excluding vectors supported only on digit 0. Every
such vector is realized by a positive numeral of that length by putting
a nonzero digit first. Its size is

\[
|\mathcal H_b|=\binom{b+D-1}{b}-D.               \tag{4}
\]

Evaluation sends every vector into (2), so `T` preserves this carrier.
Its image `I_b=E(H_b)` is a still smaller invariant integer set.

**Cycle preservation is exact.** On a histogram cycle, `E` cannot identify
distinct vertices: equal evaluations would have equal next histograms,
contradicting distinct vertices of a simple cycle. Equations (3) then
give a bijection, preserving periods, between the histogram cycles and
the integer cycles. The discarded digit order cannot affect the next
integer, but it does affect the original source: 169 and 196 have the
same histogram and next value, while remaining different integers.

This compression has a direct connection contract: source integers map to
histograms; immediate future and terminal cycle are preserved; source
identity and the first-step clock are lost; retain the original numeral
and its evaluated edge when exporting a certificate.

The corresponding Collatz quotient needs a different predicate. An
unordered multiset of valuations does **not** determine the next state:
words `(1,4)` and `(4,1)` have the same length and cost but carriers
`(9n+5)/32` and `(9n+19)/32`, with different native guards. The digit-order
quotient is licensed by (3); the analogous carry-forgetting quotient is
not. This explains which part of the attractive finite-carrier mechanism
can transfer and which part needs a new proof.

## 4. The complete portraits, including every cycle

Two independent computations agree: enumeration of the entire histogram
carrier, and a direct integer functional-graph computation on every
integer in (2). The latter uses the digit recurrence directly, without
histograms.

| Base | D | Integer core bound | Histograms | Integer image size | Terminal cycles |
|---|---:|---:|---:|---:|---:|
| 6 | 5 | 480 | 205 | 97 | 4 |
| 10 | 8 | 2,540,160 | 19,440 | 8,503 | 7 |

In base 6 the only cycles are the four fixed points `1,2,25,26`, where
the labels are decimal. In base-6 notation the last two are `41_6` and
`42_6`, explaining `4!+1!=25` and `4!+2!=26`.

The decimal cycles are exactly

\[
\begin{gathered}
(1),\ (2),\ (145),\ (40585),\\
(169,363601,1454),\ (871,45361),\ (872,45362).
\end{gathered}
\]

The finite carrier therefore retains four fixed components and three
nontrivial cyclic components. Identifying only the four fixed points
would leave genuine terminal obligations unresolved.

For a finite functional graph, its cycle polynomial is
`prod_cycles(1-z^period)`: trees contribute only zero eigenvalues to the
transition operator. The two portraits give `(1-z)^4` and
`(1-z)^4(1-z^3)(1-z^2)^2`. This spectral summary preserves cycle periods
but loses basin sizes and source membership; the histogram/evaluation
maps and actual source are still needed for a certificate.

### A genuine global rank, rather than a finite experiment extrapolated

Let `tau(h)` be the distance of a histogram to its terminal cycle in the
finite carrier. For a positive integer `n`, put
`e=max(0,digits_b(n)-(D-1))`. The following nonnegative integer rank is
explicit:

\[
R(n)=
\begin{cases}
(|\mathcal H_b|+1)e,&e>0,\\
0,&n\text{ is on an enumerated integer cycle},\\
1+\tau(H(n)),&\text{otherwise}.
\end{cases}                                      \tag{5}
\]

It strictly decreases at every step outside the cycle set. In the first
case the digit count falls, and the multiplier dominates the finite
interior rank. Inside the carrier, `tau` falls by one until a histogram
cycle; a noncyclic integer with a cyclic histogram lands on the associated
integer cycle in one more step. This proves termination of the certificate
emitter for **every** supplied positive integer.

## 5. Complete the model first, then compute its residual honestly

Choose one representative `c` from each cycle not containing 1. Define
the modified map `G` to equal `F_b` except that `G(c)=1` at those marked
representatives. The proven finite-core entry and complete portrait show
that **every positive integer reaches 1 under G**.

This construction uses the user's proposed direction: begin with a
completed modified dynamics, then compare it with the original one.
Let `V(n)` be its exact number of steps to 1. On the finite integer core
the original-map residual is

\[
\delta(n)=\mathbf1_{n\ne1}+V(F_b(n))-V(n).        \tag{6}
\]

Away from the marked representatives the residual is zero. At a marked
cycle of period `p`, `V(c)=1` and `V(F_b(c))=p`, so **delta(c)=p**.
The exact remaining obligations are:

| Base | Marked source : residual |
|---|---|
| 6 | `2:1, 25:1, 26:1` |
| 10 | `2:1, 145:1, 169:3, 871:2, 872:2, 40585:1` |

The implementation computes `V` at every core integer and verifies that
these are the entire defect support. Its maximum in the core is 8 for
base 6 and 59 for base 10; no uniform all-integer deadline is inferred.

For any function `V`, not merely the chosen one, summing (6) around a
nonroot original cycle gives

\[
\sum_{n\in\text{cycle}}\delta(n)=p.              \tag{7}
\]

The potential terms telescope. Consequently these defects cannot all be
removed by changing weights while keeping the original dynamics. The
modified completion is valid; the proposed transfer to original arrival
at 1 is false, with 145 as a one-state hostile. Exactly one edge per
nonroot cycle is necessary and sufficient to make this finite portrait
all-to-1: without altering an edge on a cycle, that cycle persists.

## 6. What this changes in the Collatz strategy

Starting with completion is useful when it produces **marked, checkable
repair obligations**. For Collatz, the earlier common-future graph allows
proof-equivalent routing edges, and [the new deletion-budget
compiler](collatz_completion_paper24_23_20261007b.md) distinguishes the
depth needed for payment from the depth actually available in a supplied
parent. The [new two-anchor phase
transfer](collatz_completion_anchor_20261007b.md) adds actual guarded
receipts in a previous entry complement. These are legitimate repairs
because their actual integer words authenticate the replacement.

Three separate interfaces remain:

1. **Entry:** prove every supplied source reaches the proposed core or a
   smaller authenticated obligation. The factorial digit bound supplies
   this here; no such global Collatz bound has been supplied.
2. **Complete terminals:** retain all terminal components, including
   nontrivial cycles or missing branches, rather than only recognizable
   fixed points. A finite Collatz abstraction may also have spurious
   cycles; deciding their actual realizability is an additional task.
3. **Repair:** replace each artificial completion edge by an authenticated
   route or a well-founded child dependency. Keep the original source and
   the exact residual during elimination. A small average defect is not
   a proof that an individual ungrounded component disappeared.

Logically, the condition `P(n) := n=1 or P(F(n))` admits the everywhere-true
solution for any deterministic map. Actual reachability is its **least**
solution generated from 1 by finitely many backward steps. A rank or a
finite authenticated proof tree is what converts the assumed completion
into this least solution. The factorial cycles make the difference
visible without any unproved conjecture.

The productive target is therefore to compile the artificial boundaries
into explicit arithmetic receipts and prove a global entry/rank lemma.
The current Collatz family results address parts of that task; the
factorial portrait supplies a fully solved test model for every transfer.

## 7. Reproduction and controls

[Script](../../04-computation/experiments/factorion_completion_defects_20261007b.py)
and [saved output](factorion_completion_defects_20261007b.out).

```text
python -B 04-computation/experiments/factorion_completion_defects_20261007b.py
python -B -O 04-computation/experiments/factorion_completion_defects_20261007b.py
```

Normal and optimized runs agree on **3,049,061 explicit exact checks**;
saved normalized-LF SHA256:
`170e1fcf628b251ee74d46e495d54c52b84db8d65d019ad1920eae2f33ed5ee3`.
The declared universes are the full two integer cores, the full histogram
carriers and their image sets, all sources 1 through 299, and sources
`b^k-1` for `k=1,2,7,8,31,1000`. The public certificate interface
recomputes the portrait; forged terminal lists and numeric type aliases
are rejected. No certificate-discovery cutoff is used in production:
termination follows from (5).
