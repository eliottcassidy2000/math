# An exact transition from feedback budget \(k\) to \(k+1\)

Date: 2026-10-20. **PROVED below:** a finite recurrence on capped weighted
core orbits, using one actual arc reversal, local chain splitting, and
canonical merging. **FINITE-EXACT:** its first two transitions agree as
sets with the independently compiled budget-one and budget-two classes.
No external novelty claim or bounded scalar recurrence for A000568 is made.

The recurrence starts from one capped core and requires no independent
input census of the next budget's cores. It does retain an increasing
finite collection of weighted core objects. This is a structural budget
recurrence, distinct from the already implemented cycle-response recurrence
that computes unrestricted A000568 through order 200.

## Inheritance, corrections, and the live board

The closest mechanism is [the fixed-feedback compiler](../defect_boundary_20261020/sparse_theory.md):
canonical transitive clone chains, exact weighted feedback straightening,
the cap \(k+1\), and the sharp core bound \(4k+1\).
[The second-coefficient classification](../defect_boundary_20261020/next_coefficient.md)
shows precisely where omitted gaps and newly shared defect endpoints first
matter.

The incoming [future-response packet](../tournament_future_boundary_20261020/README.md)
distinguishes three carriers: Boolean extension support, numerical
continuation weights, and marked object transport. Here the carrier
preserves actual arc-reversal descendants up to the target cap.
Its two-run law is a nearby hostile: the same weighted triangle can forget
the original run split and endpoint ports.

The corrected historical clique-separator scalar division is not inherited.
It failed because the boundary occupancy state was discarded. The
analogous tempting shortcut here would keep only the number of heavy
coordinates, or reverse an entire quotient arc. Both fail by explicit
four-vertex controls below.

The live board is: one actual arc; two marked endpoints; chain offsets;
cap refinement; canonical merging; intrinsic minimum feedback.
The anchor is the exact budget transition. The niche is the difference
between an individual lift and its capped family. The wildcard is the
sharp four-vertex jump in core size. No Collatz consequence is inferred.

## 1. States and their semantics

Let \(f(T)\) be the minimum number of backward arcs in an ordering of the
vertices of a tournament \(T\). Write

\[
 T=Q[w_1,\ldots,w_b]
\]

for its canonical quotient by maximal transitive clone chains. Thus
\(Q\) has no homogeneous pair, and each vertex \(i\) is replaced by an
ordered transitive chain of positive length \(w_i\).

For \(K=k+1\), define the capped state

\[
 \kappa_k(T)=[Q;c],\qquad c_i=\min(w_i,K),
\]

where brackets mean isomorphism of the weighted core. Let \(S_k\) be the
finite set of such nonempty states with \(f(T)\le k\).
The inherited cap theorem says that this is equivalent to
\(f(Q[c])\le k\), and every expansion of a coordinate with \(c_i=K\)
remains admissible.

A coordinate \(c_i<K\) specifies an exact chain length. A coordinate
\(c_i=K\) specifies the entire range \(w_i\ge K\); it does not specify
one chain of length exactly \(K\).

The empty tournament is handled separately. For nonempty tournaments,

\[
 S_0=\{[K_1;1]\}.
 \tag{1}
\]

This single capped state represents every nonempty transitive tournament.

## 2. One arc reversal needs at most four additional core vertices

Let \(x,y\) be the two vertices whose arc will be reversed.

- If they belong to different chains, split each selected chain at its
  marked vertex into an initial gap, the singleton mark, and a final gap.
- If they belong to the same chain, place them in their intrinsic chain
  order and split into an initial gap, singleton \(x\), middle gap,
  singleton \(y\), and final gap.
- Leave every other chain intact.

Empty gaps are omitted. Every nonempty gap remains a transitive module:
the only modified edge joins the two singleton marks. The tournament
after reversal is therefore an inflation of this refined tournament,
with at most \(b+4\) vertices. Normalize its clone chains to obtain its
canonical weighted core.

This is exact even when the new maximal chains join pieces formerly
belonging to different parent chains. Those mergers must be computed;
they are not guessed from the source quotient.

In particular, if \(b'\) is the new core order,

\[
 b'\le b+4,\qquad |b'-b|\le4.
 \tag{2}
\]

The second inequality applies the first one also to the inverse reversal.
It is sharp: start with a transitive tournament on five ordered vertices
and reverse the arc from the second to the fourth. The child is
\(K_1\to C_3\to K_1\), with core order five rather than one.

The same local construction works without retaining a global optimal
order. Only the intrinsic orders inside the selected parent chains are
used.

## 3. Exact finite refinement of an unbounded parent weight

Put \(L=k+2\), the target cap. A selected parent chain contains
\(s=1\) or \(s=2\) marked singleton endpoints. Let
\(g_0,\ldots,g_s\) be the lengths of the intervening gaps.
For an untouched chain use \(s=0\), so there is just one gap.

Enumerate gap caps

\[
 d_j=\min(g_j,L)\in\{0,\ldots,L\}.
\]

Keep precisely the tuples satisfying

\[
\begin{cases}
 s+\sum_j d_j=c_i,&c_i<K,\\[2mm]
 s+\sum_j d_j\ge K,&c_i=K.
\end{cases}
\tag{3}
\]

These conditions are exact:

- For an old light chain, the total length is smaller than \(K<L\), so
  no gap can reach the target cap; all entries are exact.
- For an old heavy chain, if a gap reaches \(L\), its lower bound alone
  exceeds \(K\). If none reaches \(L\), every gap is exact. Therefore
  capping the gaps neither creates nor removes a solution to the old
  lower-bound condition.
- Conversely every retained tuple has a realizing parent: take each
  \(g_j=d_j\). Entries equal to \(L\) may additionally be increased.

An untouched heavy coordinate therefore has exactly two target options,
\(K\) and \(L\). This refinement is necessary even if no arc is reversed.

Build the bounded refined tournament using weights \(d_j\) on nonempty
gaps and weight one on marks. Reverse the marked singleton arc. Then
normalize clone chains and cap their summed weights at \(L\).

Capping commutes with those sums:

\[
 \min\!\left(L,\sum_j g_j\right)
 =\min\!\left(L,\sum_j\min(L,g_j)\right).
\tag{4}
\]

The normalization depends only on the refined unweighted tournament,
which is unchanged by increasing a nonzero gap length. Thus this finite
construction returns exactly the capped descendants of all lifts of the
parent state.

**Quantifier boundary.** A capped parent represents a family. The
successor relation is existential over its actual integer lifts and
their endpoint choices. It does not assert that every individual lift
has the same descendants.

## 4. The recurrence and its completeness proof

For \(s\in S_k\), let \(R_k(s)\) contain:

1. every no-reversal target-cap refinement of \(s\);
2. every target from choosing two distinct actual endpoints, distributing
   them among one or two parent chains, applying (3), reversing their arc,
   and normalizing as above.

Use one weighted-isomorphism representative for each output. Then

\[
 \boxed{\displaystyle S_{k+1}=\bigcup_{s\in S_k}R_k(s).}
 \tag{5}
\]

**Soundness.** Every branch is realized by some parent \(T\) with
\(f(T)\le k\), by (3) and the old cap theorem. Reversing one actual arc
changes the cost of each total order by exactly one in absolute value,
so

\[
 |f(T')-f(T)|\le1.
\]

Hence \(f(T')\le k+1\). Normalization preserves the tournament; the target
cap theorem preserves membership in the budget class.
Carry branches are sound for the same reason without the reversal.

**Completeness.** Take any state in \(S_{k+1}\) and its finite representative
\(U=Q[c]\), capped at \(L\). If \(f(U)\le k\), its old cap state belongs
to \(S_k\) and a carry branch recovers the target cap. Otherwise
\(f(U)=k+1\). Choose an optimal order and reverse one of its backward
arcs. The resulting parent \(T\) has an order with cost \(k\), hence its
old cap state lies in \(S_k\).
Mark the endpoints of the reversed arc in the canonical chains of this
parent. Their actual gaps satisfy (3); applying the local construction
and reversing the arc back returns the chosen target state. QED.

The choice of a canonical source representative loses no branches:
weighted-core isomorphisms carry chains and their intrinsic orders to
chains and intrinsic orders. Every chain pair and every allowed offset
is enumerated. Child isomorphism is nevertheless a separate quotient,
and distinct branches or distinct parents can have the same child.

Every stage is finite. If the parent core has \(b\le4k+1\) vertices,
a deliberately crude bound on raw branches is

\[
 2^b\left(1+\frac{b(b+1)}2(k+3)^4\right).
\]

There are at most four selected gap coordinates, each with values
\(0,\ldots,L\); untouched old heavy coordinates each have two refinements.
The bound is not a complexity estimate for the complete iteration:
the number of parent states also grows.

### Optional fixed-order guard

For a requested tournament order \(N\), let \(m\) be the sum of the
refined gap caps and the singleton marks, **before** clone normalization.
That branch has a realizing actual tournament of order \(N\) exactly when

\[
 m\le N\quad\text{and}\quad
 \bigl(m=N\ \text{or some gap cap equals }L\bigr).
 \tag{6}
\]

If no gap reaches \(L\), every length is exact and the only mass is \(m\).
If some gap reaches \(L\), increasing it realizes every integer mass at
least \(m\). Such an increase can only occur in an old heavy chain, whose
lower bound remains satisfied. The flip preserves actual length.

A target canonical chain can become heavy by merging several fixed
fragments even when no source gap reaches \(L\). This does not make that
particular raw branch realizable at larger orders; (6) must therefore be
applied before normalization. For a fixed \(N\), the initial state still
represents the transitive tournament on \(N\) vertices, and this guard
retains every length-\(N\) branch needed by the completeness proof.

## 5. Ordered defect grammars versus intrinsic minimum budget

For a *marked ordered* tournament, its defect graph records precisely
which arcs point backward. Its active vertices are the defect endpoints,
and the gaps count all other vertices. Adding one previously absent
defect has three cases: join two old endpoints; extract one new endpoint
from a gap and join it to an old one; or extract two new endpoints from
one or two gaps. This is a reversible exact grammar before quotienting.

The displayed defect count is not necessarily \(f(T)\). Different orders
of one tournament may have different defect graphs and different counts.
Canonical clone normalization also discards the displayed global order.
Equation (5) bypasses the need to preserve that order: soundness uses the
one-arc Lipschitz bound, while completeness uses an optimal order only
as an existence witness.

It would be invalid to count raw ordered descriptions and divide by a
constant. The recurrence is a union of weighted-core isomorphism classes,
not an arithmetic sum of its branch counts.

## 6. Exact checks and hostile controls

[budget_check.py](budget_check.py) computes (5) starting from (1), without
using the next budget's endpoint compiler as input. The latter is loaded
only for an independent equality test.

| Transition | Parent states | Raw carry | Raw same-chain flips | Raw cross-chain flips | Distinct target states |
|---|---:|---:|---:|---:|---:|
| \(0\to1\) | 1 | 2 | 27 | 0 | 20 |
| \(1\to2\) | 20 | 51 | 2,624 | 4,144 | 2,691 |

Both resulting sets agree exactly with the independent finite compiler.
Every target also passes an intrinsic weighted-feedback computation.
The literal refinement construction is independently compared with
reversing every actual arc in every tournament class through order six:
991 controls, including both sharp core changes \(+4\) and \(-4\).
The fixed-order refinement guard has 1,100 additional exact boundary
controls at budgets zero through two.

Four hostile controls identify what cannot be dropped.

1. **A cap is not an exact length.** The old budget-zero state
   \([K_1;1]\) includes every transitive chain. Length three permits the
   same-chain reversal producing \(C_3\); the literal one-vertex
   representative has no arc to reverse. Length two additionally requires
   the untouched-heavy refinement from cap one to cap two.
2. **The endpoint offset matters.** In \(C_3[2,1,1]\), write the doubled
   chain as \(x_0\to x_1\), and let \(y\) be the next vertex on the directed
   triangle. Reversing \(x_0y\) gives core order three; reversing \(x_1y\)
   gives core order four.
3. **A quotient arc is many actual arcs.** Reversing the entire quotient
   arc in that example gives a transitive tournament, with core order one.
   It is neither of the two singleton-arc transitions above.
4. **Heavy locations matter.** Fix the same unweighted core
   \(Q=C_3\to\text{sink}\). Give one vertex cap weight two and the other
   three cap weight one. Both choices have core order four, one heavy
   vertex, and total cap mass five. A heavy sink has 66 distinct target
   states; a heavy triangle vertex has 95.

At \(1\to2\), the number of distinct parents of a child is not constant:

| Distinct parents | Number of target states |
|---:|---:|
| 1 | 570 |
| 2 | 2,026 |
| 3 | 32 |
| 4 | 43 |
| 5 | 15 |
| 6 | 4 |
| 19 | 1 |

The large transitive cap can absorb reversals from many different parents.
This supplies a concrete obstruction to a uniform parent-divisor rule;
it does not preclude a richer weighted incidence formula.

## 7. Counting consequence and remaining compression target

The inherited fixed-budget counting formula now has a recursive input bank:

\[
 A_k(x)=1+\sum_{[Q;c]\in S_k}
              \frac{x^{\sum c_i}}{(1-x)^{\#\{i:c_i=k+1\}}}.
\]

Thus (1), (3), and (5) give a complete structural recurrence generating
every \(A_k(x)\). Every tournament of order \(n\) has feedback at most
\(\lfloor n(n-1)/4\rfloor\), so unrestricted A000568 is recovered from
that budget layer. This is exact but does not assert computational
superiority to the already available cycle-response method.

The unresolved compression problem is specific: replace the union of
individual target weighted-core orbits by a smaller state that preserves
the full successor overlap pattern. Counts of parents, heavy coordinates,
or total weight do not suffice. Endpoint locations and canonical merger
data are the sidecar restored by the present construction.

Reproduce:

```text
python 05-knowledge/results/defect_transitions_20261020/budget_check.py
python -O 05-knowledge/results/defect_transitions_20261020/budget_check.py
```

All checks use explicit exceptions rather than assertions. Normal and
optimized runs produce identical output.
