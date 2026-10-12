# Eleventh-root thirds contexts: keep the ports and test the whole relation

2026-10-12. **PROVED, scoped and self-audited:** marked degree elimination,
the full lattice colouring, and the displayed arithmetic identities.
**FINITE-EXACT:** complete geometric graphs and their complete six-colour
boundary relations. **OPEN:** six-colour forcing in Q(zeta55), the full
field chromatic number, and the plane six/seven question. No independent
peer audit, literature priority, or formalization is claimed.

## 1. Inheritance and the changed test

Put L=Q(zeta55), alpha=zeta5, gamma=zeta11, beta=sqrt(-11),
K=Q(alpha,beta). The original marked ports remain exactly

\[
 P_j=\sum_{i<j}\alpha^i\quad(0\le j\le 4),\qquad
 R=(0,-\eta,3/4,-\bar\eta),\quad\eta=(-3+\sqrt{-55})/8.
\]

Their union B has eight distinct points. Its complete graph is a C5 and a
C4 meeting at 0. The [one-layer theorem](eleven_port_layers_20261011.md)
already rules out boundary strengthening by arbitrary outside-K satellites
on at most five K-ports. The [component theorem](cyclotomic_unit_component_20261011.md)
already proves universal unit-path coverage of its exact localization,
including 1/3, and the [field barrier](rhombus_overlap_elliptic_20261011.md)
already gives chi(K)<=5. We do not present a new reachability theorem or
field escape as six-colour forcing.

The nearest positive mechanism is complete boundary equality-partition
gluing. The hostile is an actual separating colouring surviving a large
field-escaping graph. The corrected near miss is treating failure of one
additive colour map as a chromatic lower bound. The least-used sidecar is
now a **colour-extension elimination order with all ports frozen**.

The incoming [return-thickness packet](marked_return_thickness_20261012a.md)
explains why six ports can retain 203 equality patterns whereas two ports
only retain equality/separation. We keep the eight existing ports instead
of replacing them by a convenient pair. The new test also retains all
cross-edges when two separately harmless contexts are united.

The shared concept board and Collatz inheritance are recorded in the
[source 127 companion](collatz_paid_exit_127_20261012.md). The analogy is a
checked elimination with a witness lift. A graph-colour extension is not a
Collatz convergence theorem.

## 2. An exact port-preserving elimination rule

Let G be a finite simple graph with marked vertices B and palette k. If
v is unmarked and its current degree is less than k, then

\[
                     R_k(G;B)=R_k(G-v;B).                \tag{1}
\]

One inclusion is restriction. For the other, at most k-1 colours occur on
the already coloured neighbours, so a spare colour extends the colouring
to v. Keep the deleted vertex, its actual neighbour list, and its position
in the deletion order. Reversing the order supplies a colouring lift for
every retained boundary assignment.

Iterate until every remaining unmarked vertex has degree at least k.
This is the usual degree-peeling argument, used here relative to immutable
ports. It preserves the **entire boundary relation**, not merely chromatic
number or the number of colourings. It also need not preserve extension
counts: those are a different invariant.

For k=6, an empty interior after peeling proves every proper boundary
colouring extends. A nonempty interior is only a necessary place to search
for further forcing. It is not a lower-chromatic certificate. The sharp
abstract hostile is a hub joined to six differently coloured marked leaves:
the degree-six hub must remain and really blocks that assignment. With
five leaves it can always be deleted. These are abstract controls of the
lemma, not claimed six-colour forcing graphs in L.

The cost model is explicit: a deletion consumes one finite vertex and
exports a guaranteed extension, whereas a Collatz dependency consumes a
proved rank decrease and needs a grounded child. In both cases ignoring
the protected source or boundary would change the question.

## 3. The previous mixed-generator graph has the entire bare relation

The [previous graph](eleven_port_layers_20261011.md) has 8,436 vertices and
24,238 complete unit edges, with generator

\[
 u=(1-\alpha-\gamma)/(1-\alpha^{-1}-\gamma^{-1}).
\]

Its coordinates, denominator and complete edge hashes are unchanged here.
Freeze B and peel at k=6. Exactly874 vertices remain, including the eight
ports. This removes 7,562 vertices without losing any boundary information.

Enumerate every canonical equality partition of B using at most six
colours and respecting all nine actual port edges. The complete counts,
classified by the number of occupied colours, are

| Used colours |3|4|5|6|
|---|---:|---:|---:|---:|
| Proper equality partitions |30|180|247|113|

There are 570 in total. Every one has a proper six-colouring of the retained
874-vertex core. The output stores all 570 explicit core colour vectors.
Every vector is independently restored through the deletion order and
checked on **all 24,238** original edges. The SAT solver supplies positive
witnesses only; no unrecorded UNSAT conclusion is used.

Thus this complete, three-layer mixed-generator graph has exactly

\[
                       R_6(G;B)=R_6(G[B];B).             \tag{2}
\]

It preserves all 327,600 labelled proper six-colour assignments. Indeed
570 counts equality patterns, while labelled assignments multiply each
column by 6!/(6-j)!; alternatively the two cycle counts give 3120*630/6.
A three-colouring alone would not have proved (2). Keeping a proper
six-colouring of an arbitrary precoloured boundary is the stronger test.

## 4. A denominator-three unit and its nonlinear colour survivor

The old [field-escape construction](unit_defect_frontier_20261011f.md),
section 3, already uses b=(5+sqrt(-11))/6. We recover its sign/conjugate

\[
                  v=(-5+\beta)/6=-\bar b,
\qquad v\bar v=1,\qquad 2+v+\bar v=1/3.                \tag{3}
\]

So this is an inherited unit with a new retained connector use, not a new
rotation discovered from scratch. The five fixed points

\[
                  0,1,2,2+v,1/3
\]

form an induced four-edge unit path. All six nonconsecutive pairs are
checked not to have distance 1. It gives an explicit short path for a
point whose reachability was already known from the stronger component
theorem. No minimal path-length claim is made.

There is an eleventh-root presentation of the same object. The factors
of Phi11(X) modulo 3 are

\[
 X^5-X^3+X^2-X-1,\qquad X^5+X^4-X^3+X^2-1.
\]

For the first factor f, N(f(gamma))=3^5, and exact reduction gives

\[
 \frac{f(\gamma)}{f(\bar\gamma)}
 =\frac{\gamma^7-2\gamma^5+\gamma^4+\gamma^2+\gamma+1}{3}
 =\gamma^5\bar v.                                    \tag{4}
\]

The norm and last equality are checked exactly. Equation (4) identifies
the construction as a root rotation of a K-unit. It can leave K through
that root; it does not supply a chromatic increase. This simplification
retains useful information about the construction. It exits the earlier coefficientwise
3-integral additive-colouring domain.

That exit does **not** force four colours. Let

\[
                  \Lambda=\{A/3+Bv:A,B\in\mathbb Z\}.
\]

For a displacement of this form,

\[
              |A/3+Bv|^2=((2A-5B)^2+11B^2)/36.          \tag{5}
\]

The only unit coordinate directions are

\[
                 \pm(3,0),\quad\pm(0,1),\quad\pm(5,1).
\]

To prove completeness, (5) forces |B|<=1, after which the integer equation
has exactly those six solutions. Every additive map Lambda->F3 vanishes
on the unit displacement 1=3*(1/3), so no additive three-colouring exists.
But the **nonlinear** periodic rule

\[
 c(A,B)=c_5(A\bmod 5)+B\pmod 3,\quad c_5=(0,1,1,2,0)      \tag{6}
\]

colours every unit edge in the full infinite lattice. Steps(0,1) and(5,1)
change the colour by 1; step(3,0) follows a properly coloured C5 in the
five-entry list. The odd closed unit walk 5+3v+3bar(v)=0 has length 11, so
the lattice is not bipartite. Consequently chi(Lambda)=3 exactly.

Equations(3)–(6) give a concrete all-depth hostile to “defeat the linear
colouring, therefore force more colours.” The missing coordinate in that
argument is the possibility of a nonlinear palette rule. No colouring of
all L is inferred from this lattice.

## 5. Deeper eleventh-root contexts and their actual union

Let

\[
 X=\{a/3+bv:-4\le a\le 4,\ -1\le b\le 1\},\qquad
 H=B\cup\{p+\gamma^j x:p\in B,\ 0\le j<11,\ x\in X\}.
\]

Identify equal coordinates and authenticate every unit pair. This preserves
all eight original ports and uses points beyond their immediate unit
circles. The graph has 2,278 vertices and 6,405 edges, with an explicit
three-colouring. More strongly, its entire unmarked interior peels at
palette 6: its relative core is exactly B. Hence every proper six-colouring
of B extends, with no search over colour assignments needed for the proof.
The checker nevertheless verifies all 570 canonical assignments and their
full lifts as a cross-control.

Next form the **actual geometric union** of H with the old graph G.
They share 116 points, not just the eight ports. Its complete graph has

\[
 10,598\text{ vertices},\quad 30,506\text{ unit edges},
 \quad 88\text{ edges absent from the two input edge sets}.          \tag{7}
\]

Those88 new constraints are recorded explicitly. They cannot be discarded
by applying abstract gluing only at B. A complete norm computation retains
them, along with every input edge and every coordinate identification.
An independent positive three-colouring is checked on the whole union.

Even after those seams are added, degree peeling leaves **exactly the same
874 vertices and the same induced core edges as G**. Thus (1) and the 570
stored witnesses prove

\[
                      R_6(G\cup H;B)=R_6(G[B];B).        \tag{8}
\]

This is a changed-generator composition with real extra edges, whose
failure is explained by a complete extension certificate. It is stronger
than merely observing that its chromatic number stays three. Both isolated
and composed candidates preserve every boundary pattern needed to defeat
the proposed forcing conclusion.

The next useful supplier must alter this retained six-core relation.
A new denominator, a new field generator, a long path, or extra cross-edges
can still be entirely absorbed by the degree-extension mechanism. First
check whether an alleged stronger context survives that cheap test; then
ask it to exclude one specific surviving boundary partition while retaining
a protected lift. Merely increasing the point count is not that test.

## 6. What2047 and 127 actually encode here

The identity 2^11-1=23*89 has a precise eleventh-root role:

\[
 N_{\mathbb Q(\gamma)/\mathbb Q}(2-\gamma)=\Phi_{11}(2)=2047.
\]

Both primes have exact order 11 for 2; both are 1 mod 11 and split completely
in Q(gamma). Directions

\[
 u_j=(2-\gamma^j)/(2-\gamma^{-j})\quad(1\le j\le 10)
\]

therefore have denominators naturally supported on this norm. The exact
sum and a constructive denominator identity are

\[
 \sum_{j=1}^{10}u_j=14849/2047,\qquad
 1799-248\sum_{j=1}^{10}u_j=1/2047.                     \tag{9}
\]

The latter is a4,279-edge unit walk, without any optimality claim. It adds
an explicit receipt to the already known localization theorem. The trace
formula follows by summing the finite geometric inverse of 2-gamma^(-j);
the checker also computes each of the ten directions independently.

In the separate Collatz graph, the actual 2047 trajectory reaches 127 after
41 odd steps, total valuation 69, and odd charge 27. Its marked prefix is
saved. Multiplying23 and 89 is not a rule composing their Collatz proofs;
(9) supplies no such multiplicativity or conjugacy.

The other prompt has a different, equally exact interpretation:

\[
127=10^2+3^3=N_{\mathbb Q(\sqrt{-3})/\mathbb Q}(10+3\sqrt{-3}).
\]

It gives the norm-one rotation

\[
                         (73+60\sqrt{-3})/127.          \tag{10}
\]

But sqrt(-3) is **not** in L: adjoining zeta3 doubles the cyclotomic degree
from 40 to 80. Thus this attractive 127-based rotation leaves the prescribed
field and its triangle-free geometry. It cannot be used as an L forcing
step merely because 127 also labels the Collatz source. This is a concrete
field-side guard, not a dismissal of the numerical identity.

## 7. Scope, reproducibility, and the connection ledger

| Source object | Map and preserved predicate | Lost data / retained sidecar |
|---|---|---|
| Geometric point set | Complete unit graph | Exact coordinates and all cross-edges |
| Marked graph | Degree-below-six elimination | Reverse colour lift, original ports |
| Old plus thirds contexts | Actual coordinate union |116 coincidences and 88 new edges |
| Unit factor modulo 3 | Ratio in(4) | Root rotation versus true field enlargement |
| Additive colouring obstruction | Relation3*(1/3)=1 | Nonlinear colourings such as(6) survive |
| Mersenne norm 2047 | Unit denominators and trace | No Collatz orbit homomorphism |
|127 as square plus cube | Eisenstein norm(10) | Coordinate field changes |

[Program](../../04-computation/experiments/paid_exit_boundaries_20261012.py),
[complete output](paid_exit_boundaries_20261012.json). All coordinate arithmetic
uses the exact 40-dimensional power basis of Q(zeta55). The two modular
filters at 331 and 661 are authenticated ring maps and cannot discard a true
edge. Every survivor receives an exact full norm check. A separate small
unfiltered control covers all marked ports, the connector and the entire
small lattice patch. The large old graph's hashes agree with its inherited
artifact before any new boundary assertion is accepted.

The output retains all 570 old core colour vectors, every deletion order,
complete graph hashes, the 88 new seam edges and positive colour lifts.
Ordinary and `-O` runs agree byte for byte, with 17,578,663 always-active
checks in the new driver and 36,933 in the imported exact geometry engine.
The Collatz companion specifies its compiler and ledger universes.
Reproduction commands, from the repository root:

```text
python3 -B 04-computation/experiments/paid_exit_boundaries_20261012.py
python3 -B -O 04-computation/experiments/paid_exit_boundaries_20261012.py
```

Source SHA256:
`04a8d928662442e92dd3543b81cc783696e8d9d72222f4b6e948ec26bdded9c6`.
Output SHA256:
`631a879f09344a82b3869e10a0f537f0f42ba57f50425efedea3f77f3270b30e`.
These counts include repeated checking of every edge in every retained
colouring lift; they do not strengthen the quantified scope of the proofs.

The next geometric obligation is a realizable context that removes a
specific surviving six-colour port pattern after this exact reduction.
The independent charge task similarly needs an authenticated reduction
at every remaining original source. The shared research move is to test
whether the new operation changes the full retained obligation; the two
mathematical predicates and their unresolved universal steps remain distinct.
