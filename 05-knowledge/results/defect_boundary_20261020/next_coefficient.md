# The second coefficient of fixed-feedback tournament counts

**PROVED below**, with **FINITE-EXACT** independent template checks through
budget four. No external novelty claim. This extends the fixed-budget
classification in [sparse_theory.md](sparse_theory.md); it is not a recurrence
for unrestricted A000568 at bounded cost.

Let \(a_{n,\le k}\) be the number of isomorphism classes of tournaments with
minimum feedback arc count at most \(k\). For every fixed \(k\ge1\), the
eventual counting polynomial begins

\[
 \boxed{\displaystyle
 a_{n,\le k}
 =\frac{n^{2k}-3k\,n^{2k-1}}{2^k k!}
       +O_k(n^{2k-2}).}
 \tag{1}
\]

For \(k=1\), this is \((n^2-3n)/2+O(1)\); for \(k=2\),
\(n^4/8-3n^3/4+O(n^2)\), agreeing with the separately compiled
all-order formulas. The proof identifies every capped core orbit that
can contribute to these two coefficients.

## Inherited mechanism and the precise remaining issue

Put \(K=k+1\). A capped core has a pair-free tournament \(Q\), positive
weights \(c_v\le K\), and weighted minimum feedback at most \(k\).
A vertex of weight \(K\) is called heavy.

The previous proof gives:

- a canonical homogeneous-chain quotient with positive weights;
- a cap stabilizer fixing every heavy vertex individually;
- an orbit contribution \(x^v/(1-x)^h\), where \(h\) counts heavy vertices
  and \(v\) is the sum of capped weights;
- \(h\le2k+1\);
- exactly \(M=(2k-1)!!\) orbits with \(h=2k+1\), each having
  \(v=2k^2+5k+1\).

Therefore the next coefficient requires only the number of orbits with
\(h=2k\). The delicate point is injectivity when a heavy separator is
removed or demoted: the few neighboring light vertices might acquire
internal symmetries. The directed-triangle case below is a real symmetry,
but it does not identify distinct templates.

## 1. Complete classification of the \(h=2k\) stratum

Choose a core order with weighted backward cost at most \(k\).
No backward arc meets a heavy vertex. Mark the endpoints of the backward
arcs, and let their number be \(t\). Each unmarked interval is homogeneous;
pair freeness permits at most one core vertex in such an interval.
The \(2k\) heavy vertices occupy distinct unmarked intervals. Hence

\[
 2k-1\le t\le2k.
\]

If \(t=2k\), there are exactly \(k\) disjoint backward arcs, forming a
perfect matching on the ordered endpoints. Every endpoint weight is 1:
the \(k\) positive edge products already sum to at least \(k\).
There are \(2k+1\) intervening and exterior gaps; exactly one is nonheavy.
It is either occupied by a light unmarked vertex or empty.

If \(t=2k-1\), there are again exactly \(k\) backward arcs and every
endpoint has weight 1. The defect graph has no isolated vertices, total
degree \(2k\), and \(2k-1\) vertices. Thus one vertex has degree two and all
others degree one: one three-vertex path plus \(k-2\) disjoint edges.
All \(2k\) gaps are heavy.

This gives precisely the following cases.

| Family | Defect graph | Gap change | Core size |
|---|---|---|---:|
| A | Matching on \(2k\) ordered endpoints | Demote one gap from \(K\) to any weight \(1,\ldots,k\) | \(4k+1\) |
| B | Matching on \(2k\) ordered endpoints | Omit one gap, subject to the exclusion below | \(4k\) |
| C | One \(P_3\) and \(k-2\) disjoint edges on \(2k-1\) ordered endpoints | All gaps heavy | \(4k-1\) |

For family B, an internal omitted gap cannot lie between endpoints paired
to each other. Such endpoints become homogeneous, and together with the
neighboring heavy vertices normalize to one chain; the heavy count drops.
All other B templates remain pair free.

Family C is absent at \(k=1\). All three families have a displayed order
of weighted cost \(k\), so they are admissible.

## 2. Pair freeness and injectivity after losing a separator

Heavy vertices induce a transitive tournament. Its order is intrinsic and
preserved by every weight-preserving isomorphism. Each light vertex has
a unique gap position relative to that order because every arc involving
a heavy vertex agrees with the chosen good order.

In C, every internal gap contains one light vertex. Heavy order therefore
identifies every vertex, and the path-plus-matching graph is recovered
exactly. Pair freeness follows from the same intervening-heavy and
defect-partner witnesses used for the top matching templates.

In B, omitting an exterior gap leaves a single exterior light vertex.
Omitting an internal gap creates one gap with two light vertices.
If they were matched, the stated exclusion applies. Otherwise they form
a transitive pair with distinct outside matching partners. Their internal
orientation fixes their order, and each outside partner distinguishes
the pair. Thus the missing gap and the complete matching are recoverable.
All other vertices remain separated by heavies, and neighboring
light/heavy pairs are distinguished by the light vertex's partner.

In A, the underlying unweighted core is unchanged from its top matching
template, so it remains pair free. An exterior demotion produces a gap
with two light vertices; their transitive orientation recovers the order.
An internal demotion produces a gap with three light vertices:

- If the two endpoint neighbors are not matched to each other, the
  three vertices induce a transitive tournament. Its unique order recovers
  the middle unmarked vertex and both endpoint positions.
- If they are matched, the three vertices induce a directed triangle.
  Both endpoint vertices have no outside defect partner. A rotation can
  change which vertex is regarded as unmarked, but all such choices
  describe the same ordered template: one internal matching edge and
  the same outside matching. If the demoted weight exceeds 1, the unique
  nonunit weight additionally identifies the unmarked vertex.

The three-vertex symmetry therefore does not produce distinct parameter
tuples with the same capped orbit. Different demoted gap positions are
distinguished by the intrinsic heavy order and the location of the
enlarged light gap. Different demoted weights are distinguished by the
weight multiset of that gap.

Consequently the A, B, and C parameterizations are injective. Their core
sizes are different, so the families are disjoint. The completeness
argument in section 1 shows that no \(h=2k\) orbit is missing.

## 3. Counting the three families

There are \(M=(2k-1)!!\) matchings on \(2k\) ordered endpoints. Thus

\[
 N_A=k(2k+1)M.
\]

Before exclusions B has \((2k+1)M\) templates. Summing the number of
adjacent matched endpoint pairs over all matchings gives
\((2k-1)(2k-3)!!=M\). Each is exactly one excluded internal gap, so

\[
 N_B=2kM.
\]

For C, choose three vertices, choose the path center, and match the rest:

\[
 N_C=\binom{2k-1}{3}\,3(2k-5)!!=(k-1)M\qquad(k\ge2),
\]

with \(N_C=0\) at \(k=1\). Hence the total number with \(h=2k\) is

\[
 \boxed{\displaystyle
 N_{2k}=M(2k^2+4k-1).}
 \tag{2}
\]

The shifts are also explicit. Put \(B=2k^2+4k\). An A state with
demoted weight \(a\) has shift \(B+a\); B states have shift \(B\); C states
have shift \(B-1\). Therefore the entire \(h=2k\) stratum contributes

\[
 \frac{M}{(1-x)^{2k}}
 \left((2k+1)\sum_{a=1}^k x^{B+a}
       +2k x^B+(k-1)x^{B-1}\right).
 \tag{3}
\]

## 4. Extracting the second coefficient

The top \(M\) states each contribute
\(x^s/(1-x)^{2k+1}\), with \(s=2k^2+5k+1\). Their coefficient is
\(\binom{n-s+2k}{2k}\), whose first two polynomial coefficients are

\[
 \frac1{(2k)!},\qquad
 \frac{(2k+1)/2-s}{(2k-1)!}.
\]

Each state with \(h=2k\), regardless of its shift, contributes
\(n^{2k-1}/(2k-1)!\) to the next coefficient. States with fewer heavy
vertices have degree at most \(2k-2\) and cannot affect it. Therefore

\[
\begin{aligned}
 [n^{2k-1}]a_{n,\le k}
 &=\frac{N_{2k}+M((2k+1)/2-s)}{(2k-1)!}\\
 &=-\frac{3M}{2(2k-1)!}
 =-\frac{3k}{2^k k!},
\end{aligned}
\]

proving (1).

## Finite-exact controls

[next_coefficient_check.py](next_coefficient_check.py) constructs A, B, C,
and top templates independently, then applies the general compiler's
canonical clone normalization and weighted-isomorphism routines.

| Budget \(k\) | Top | A | B | C | Total with \(h=2k\) |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 3 | 2 | 0 | 5 |
| 2 | 3 | 30 | 12 | 3 | 45 |
| 3 | 15 | 315 | 90 | 30 | 435 |
| 4 | 105 | 3,780 | 840 | 315 | 4,935 |

The checker verifies injectivity after exact normalization, the distinct
core sizes, heavy counts, excluded-gap collapses, and both polynomial
coefficients. At \(k=1,2\), its top two strata agree as sets with the
independent complete finite compiler. A separate positive hostile control
finds the three actual rotations of the internal A triangle while retaining
one template.

```text
python 05-knowledge/results/defect_boundary_20261020/next_coefficient_check.py
python -O 05-knowledge/results/defect_boundary_20261020/next_coefficient_check.py
```

All checks use explicit exceptions. Normal and optimized output agree
byte-for-byte. The universal result follows from the complete structural
classification, not extrapolation from these four budgets.
