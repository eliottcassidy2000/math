# Induced count spectra obstruct arbitrary pure-birth returns

**PROVED, scoped, using a complete finite head:** none of the first five
successor returns on `405*2^k+71`, starting at881, can be accomplished by
adding vertices alone, even with arbitrary comparisons, arbitrary final
root selection, arbitrary initial same-count replacement, and unlimited
transitive padding. **FINITE-EXACT:** the fourteen finite target spectra.
The exclusion does not cover interleaved changes of existing arcs, vertex
deletions, or count-changing preparation. Chan–Pak remains **OPEN**.

Inheritance: [READY](../../two_entry_funnels_20261103/READY.md) supplies the
unique-start-parent classification. [FAMILIES](../../multicore_metagraph_20261113/FAMILIES.md)
and its `next_target.json.gz` supply the complete target heads. The older
homogeneous-inflation exclusion is weaker: the present calculation permits
every induced embedding, not just a uniformly substituted quotient. The
least-used sidecar is the *whole induced Hamiltonian-count spectrum*.

## 1. The hereditary sidecar

For a tournament T define

    C(T) = { H(T[S]) : S subset V(T) },     H(empty)=1.

If U is an induced subtournament of T, then C(U) is a subset of C(T).
In particular H(U) belongs to C(T). This is independent of vertex labels
and of the marked root. Every sequence of vertex additions retains its
initial tournament as an induced subtournament, however the new arcs are
chosen and however often the root is changed.

The scalar H alone also increases weakly under induced inclusion. For one
new vertex x, insert x immediately before the first vertex of an old
Hamiltonian path that x beats; if there is none, append x. The preceding
vertex, when present, beats x. Deleting x recovers the old path, so the
map is injective. The spectrum is stronger than this scalar inequality:
881<1691 passes the inequality and fails the spectrum test below.

This is an exact structural sidecar, not a sufficient embedding test:
matching one count or even inclusion of count spectra does not assert a
graph embedding. It discards adjacency, multiplicity, and root information.
The receipt also retains subset-size ranges and count multiplicities.

## 2. All-order unit-ready targets reduce to finite strong targets

Write F_r(T) for the number of ordered spanning two-path covers whose first
path is nonempty and ends at r; the second may be empty. Unit readiness
means F_r(T)-H(T)=1, not merely a deletion/start difference at arbitrary
indegree. The inherited classification says a strong unit-ready graph is
the canonical newborn over a strong unique-start parent (with the small
singleton boundary handled separately).

Here is the padding reduction directly. Decompose T into its totally
ordered strong components. Let A be the union strictly before the root
component C and B the union strictly after it. Then

    H(T)=H(A)H(C)H(B),
    F_r(T)=P_A(2) F_r(C) H(B),

where P_A(2) counts two labelled path bins, allowing either to be empty.
Set H(empty)=P_empty(2)=1. If A is nonempty, P_A(2)>=2H(A), while
F_r(C)>=H(C) by cutting each Hamiltonian path just after r. Consequently
F_r(T)-H(T)>=H(T), impossible for unit readiness when H(T)>1.
Thus A is empty. The remaining equality is

    1=H(B) (F_r(C)-H(C)).

Both factors are positive integers. Hence H(B)=1 and B is transitive;
C is strong and unit-ready. If C were a singleton its defect would be0,
so C is nontrivial. Therefore every target with H>1 is exactly

    C -> TT_a,    a>=0,

with its unit-ready root in C. Every induced subtournament of this graph
has the same count as its intersection with C, because the dominated
transitive factor has count1. Thus C(C -> TT_a)=C(C) for all a.

The completeness of the finite heads uses the inherited unique-start
grammar together with the strong-order minimum in
[Busch, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).
For target N the parent count is at most (N+1)/2. The inherited head records
the resulting maximum parent orders13,14,16,17,18 for the five targets
below. These bounds and the exhaustive ranker are dependencies, not newly
proved enumeration here.

## 3. Complete target spectra

Both forward extension and backward deletion dynamic programs count every
subset of every strong target. They agree entry by entry. Literal
permutation counts on two smaller induced subtournaments of each target
give a third control. The full target count and actual readiness are
recomputed; the inherited parent supplies a positive spectrum-membership
control. A C3-to-H5 birth supplies a positive growth control.

| Initial count | Ready target | Strong target classes | Largest order | Contains the initial count? |
|---:|---:|---:|---:|---|
|881|1691|3|12|No|
|1691|3311|3|13|No|
|3311|6551|4|15|No|
|6551|13031|2|15|No|
|13031|25991|2|17|No|

The receipt contains literal rows, the complete count histograms, subset
table hashes, and the inherited input hashes. No subset-size filter is used.
For example, the three1691 targets have predecessor-count gaps
(859,1099), (699,1055), and (851,959), respectively; all miss881.

**Conclusion.** For each row, no tournament of the initial count is an
induced subtournament of any unit-ready target of the specified count,
at any order. Thus even an unrestricted initial same-count preparation,
followed by any finite number of arbitrary births, cannot perform that
return. This includes an old root retained throughout and an arbitrary
old or new root selected at the end.

The five-row exclusion is not an assertion about every later405-ray step.
It also does not assert that the target counts are absent: all fourteen
strong witnesses are present in the receipt.

## 4. How a repair can escape

If an initially retained source is to disappear as an induced object,
one must delete one of its vertices or change at least one arc between
two of its vertices. Editing only arcs incident to new vertices cannot
help. A count-changing preparation before growth also escapes, because
the resulting induced source need no longer have the forbidden count.

The separate terminal-newborn budget is consistent with this: a unit-ready
newborn has a unique-start parent, so

    H_new=H_parent+F_parent-1 >= 2H_parent-1.

Reaching1691 requires H_parent<=846. Extra induced births cannot lower an
881 starting count. This inequality is inherited, not a new theorem; the
spectrum test additionally closes the arbitrary retained-root loophole.

The new selective-star repair in [repairs](../repairs/PROOF.md) deliberately
changes old arcs and lowers the parent count before its final birth. It
therefore addresses the mechanism of the obstruction rather than merely
enlarging a birth search.

## Reproduction

From the repository root:

```text
python -B 05-knowledge/results/growing_t_returns_20261113/containment/count_spectrum.py
python -B -O 05-knowledge/results/growing_t_returns_20261113/containment/count_spectrum.py --output receipt.optimized.json.gz
```

All decisions use explicit checks rather than removable Python assertions.
The int64 tables are safe here: induced counts are at most the full count
25991 by the injection above, and each endpoint count is nonnegative and
at most its subset total. This receipt does not regenerate the inherited
unique-start isomorphism census.
