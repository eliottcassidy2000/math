# Defect deficiency and bounded collision motifs

**Status: PROVED statements 1–3 below; FINITE-EXACT numerical checks.
Compressing these strata into a scalar recurrence without enumerating their
weighted-core orbits remains OPEN in this approach.**

This note uses the exact weighted straightening and cap theorems in
[the preceding sparse-defect proof](../defect_boundary_20261020/sparse_theory.md).
The new coordinate explains how shared endpoints, shortened gaps and unused
budget enter the same hierarchy. It is a presentation theorem with an explicit
canonical-core condition, not an invariant of an arbitrary unnormalized template.

## 1. An exact loss equation

Fix budget \(k\), cap \(K=k+1\), and a canonical pair-free weighted core
\(Q[w]\) with \(f(Q[w])\le k\). Choose a block order of weighted backward
cost at most \(k\), which exists by weighted straightening. Define:

- \(r\): number of distinct backward core arcs;
- \(t\): number of their distinct endpoints;
- \(h\): number of coordinates with capped weight \(K\);
- \(u\): number of the \(t+1\) endpoint gaps that are empty or contain a
  light coordinate.

Every heavy coordinate is outside the defect endpoints: a backward arc
incident with weight \(K\) would already cost more than \(k\). Every gap
contains at most one core vertex. Indeed, consecutive vertices in a gap
have no backward incidences, so they are a homogeneous pair, contradicting
the canonical core's pair-free property. Therefore \(h=t+1-u\).

**Proposition 1 (PROVED).** The deficiency of this stratum is exactly
\[
 D:=2k+1-h=2(k-r)+(2r-t)+u.                 \tag{1}
\]
All three terms are nonnegative. Writing
\[
 q=k-r,\qquad \delta=2r-t,
\]
the weighted excess obeys
\[
 \sum_{\{a,b\}\in E(D_{\rm arc})}(w_aw_b-1)\le q.             \tag{2}
\]

Here \(D_{\rm arc}\) is the ordinary graph of backward arcs on its incident
endpoints, retaining their linear order and arc identities. Equation (1)
is just \(h=t+1-u\), with the weighted budget supplying the inequalities.
Different good orders can have different \(q,\delta,u\); their sum is fixed.

The distinction from an arbitrary raw template matters. Clone normalization
can merge gaps, absorb endpoints or create a newly heavy coordinate. One must
recompute \(h\) after normalization. The checker records these changes rather
than treating the source deficiency as a conserved quantity.

## 2. Fixed deficiency means finitely many exceptional motifs

The defect graph has no isolated vertices. A connected component with
\(v\ge3\) vertices and \(e\) edges satisfies
\[
 2e-v\ge v-2\ge v/3.
\]
Summing over such components shows that all nonmatching components together
have at most \(3\delta\) vertices and at most \(2\delta\) edges.

An isolated matching edge whose endpoint product exceeds one consumes at
least one unit of the excess in (2). There are at most \(q\) such decorated
matching edges, contributing at most \(2q\) further vertices. Every endpoint
weight is at most \(q+1\), since its incident edge alone contributes at least
\(w_a-1\) to (2).

**Proposition 2 (PROVED).** At fixed deficiency \(D\), every good canonical
presentation consists of:

1. an arbitrary ordered matching of unit-weight endpoints;
2. exceptional defect components and decorated matching edges on at most
   \(3\delta+2q\le3D\) endpoints;
3. at most \(D\) exceptional gaps.

Every defect-endpoint weight is at most \(\lfloor D/2\rfloor+1\).
The nonheavy gap weights still range up to \(k\); their values and positions
are indispensable data. This is a bounded exceptional *shape* theorem.
It does not bound the number of matching patterns, identify isomorphic
presentations, or make the unrestricted recurrence scalar.

This gives a precise use for the user's four-vertex graph comparison.
Both \(P_4\) and \(K_{1,3}\) have three edges, four endpoints, and
\(\delta=2\); both first occur in deficiency two. Two disjoint \(P_3\)'s
also have \(\delta=2\). A triangle has three edges but only three incident
endpoints, so \(\delta=3\). The isolated vertex in “triangle plus isolated
vertex” is a gap coordinate, not a defect endpoint. Complementarity of
these ordinary graphs therefore does not preserve this arithmetic grading.

## 3. The complete deficiency-two menu

For \(D=2\), (1) gives exactly the following canonical presentation types:

| Unused arc budget \(q\) | Collision excess \(\delta\) | Exceptional gaps \(u\) | Defect graph |
|---:|---:|---:|---|
| 0 | 0 | 2 | \(k\) matching edges |
| 0 | 1 | 1 | one \(P_3\), plus matching edges |
| 0 | 2 | 0 | one \(P_4\), one \(K_{1,3}\), or two \(P_3\)'s, plus matching edges |
| 1 | 0 | 0 | \(k-1\) matching edges |

In the first three rows every endpoint has weight one. In the last row
all endpoints have weight one, or exactly one has weight two. Empty defect
graphs include the transitive boundary case \(k=1\). Each exceptional gap
can be absent or have any weight \(1,\ldots,k\).

**Proposition 3 (PROVED).** Generate these presentations, normalize their
clone chains, cap at \(K\), canonicalize the weighted core, and keep only
states with \(h=2k-1\). The resulting set is exactly that stratum.

Soundness follows from the displayed weighted costs. For completeness,
take any state in the stratum and choose its good canonical block order;
equation (1) places it in one of the four rows. The \(\delta=2\) list is
exhaustive because a connected nonmatching component then has at most four
vertices: on four vertices it must be a tree, hence a path or claw; the
only split possibility is two three-vertex paths.

The generator is a set recurrence, not a sum of disjoint families.
Normalization and alternative good orders can identify presentations from
different rows. The computation explicitly measures these collisions.

## Reproduction and next question

Run:

```text
python 05-knowledge/results/defect_transitions_20261020/strata_check.py
python -O 05-knowledge/results/defect_transitions_20261020/strata_check.py
```

The script generates deficiencies zero through two at \(k=1,\ldots,4\),
compares the entire retained set against the independent complete endpoint
compiler at \(k=1,2\), and recovers the previously proved leading two
coefficients. It uses explicit exceptions, active under optimized Python.
Its outputs are in [strata_check.out](strata_check.out).

The first three eventual coefficients, multiplied by \(2^k k!\), are:

| Budget \(k\) | Coefficient of \(n^{2k}\) | Of \(n^{2k-1}\) | Of \(n^{2k-2}\) | States with \(h=2k-1\) |
|---:|---:|---:|---:|---:|
| 1 | 1 | -3 | 4 | 9 |
| 2 | 1 | -6 | -1 | 264 |
| 3 | 1 | -9 | -15 | 5,365 |
| 4 | 1 | -12 | -38 | 102,750 |

The first two columns are inherited all-budget theorems and provide controls.
The independently audited [local collision proof](strata_audit.md) now gives
\[
 [n^{2k-2}]a_{n,\le k}=
 \frac{k(17-9k)}{2^{k+1}k!}.                              \tag{3}
\]
**PROVED for every \(k\ge1\)**: the bracket refers to the eventual counting
polynomial. The four numerical cases are controls; the proof classifies
all canonical identifications inside the exceptional light gaps.

At deficiency two the owner rule is explicit: a unit directed triangle
with one external defect port belongs to the matching/two-empty-gap
description; its two occupied-\(P_3\) descriptions are discarded. Surviving
mergers belong to the all-heavy \(k-1\)-matching family. All other templates
are recovered by the heavy-gap order and their local component structure.

The next target is to extend this ownership rule to deficiency three and
then to arbitrary deficiencies. The full successor-overlap quotient from
[budget.md](budget.md) remains a separate, harder compression problem.
