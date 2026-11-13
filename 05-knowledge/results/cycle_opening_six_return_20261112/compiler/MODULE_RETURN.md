# A uniform successor return preserving the physical root and incoming port

Status: PROVED in notes; FINITE-EXACT literal and independent module controls. Chan–Pak remains OPEN. This is a full-count graph replacement, not a claim that its individual arc flips preserve the count or that it factors into directed cycle-opening primitives.

## 1. The source family and the excursion

Use the seven-vertex parent

\[
P=(78,12,24,16,35,79,30),
\]

where each row is an outgoing-neighbor bit mask on vertices 0,…,6. Replace vertex 1 by a nonempty tournament M, retaining all its external incidences. Write

\[
u=H(M),\qquad v=P_M(2),
\]

where \(P_M(2)\) counts ordered pairs of possibly empty directed paths partitioning M. Add z beaten only by old vertex 4. The [full-count proof](../counts/PROOF.md) establishes

\[
H(S(M))=45v+7u,\quad \delta_z=u,\quad H(S(M)-z)=21v+3u.
\]

Here \(\delta_r=H(T-r)-s_r(T)\), with \(s_r\) the number of full Hamiltonian paths starting at r. This is the indegree-one interface defect; it is not the root's starting-path count itself. The same parent also has a two-arc neutral port switch with defect \(6v+u\), furnishing the growing-defect entry to this source.

Apply two ordinary indegree-one births, first x beaten only by z, then y beaten only by x. The same physical y is the marked root below. The inherited two-birth identity gives

\[
H(R(M))=135v+23u=3H(S(M))+2u,
\]

and its complete three-coordinate state \((H-d,\delta,d)\), where \(d=H(T-y)\), is

\[
(66v+10u,\;24v+5u,\;69v+13u).                 \tag{1}
\]

Thus the excursion has increased the defect. When M is transitive of order k, u=1 and v=2^k, so its target is exactly \(3H+2\).

## 2. The target construction and full-count proof

Start with the directed triangle \(0\to1\to2\to0\), marked at 0. Form a one-start funnel with tail

\[
M\to TT_3\to\{w\}.
\]

Every vertex of M and w points to the new funnel root x; x points to the three middle tail vertices. The triangle points to the whole tail; x points to triangle vertex 0, and the other two triangle vertices point to x. Finally add y beaten only by x. Call the resulting rooted tournament Q(M).

For M=TT_k the legal funnel word is \(1^k0001\). The inherited exact tail coefficients are

\[
L_2=9\cdot2^k,\quad B_1=9,\quad D_1=0.
\]

For example B1=9 follows by partitioning the tail before/after the root: the all-before case, one adjacent 10 cut, and seven choices whose first after-vertex is a zero and whose last before-vertex is the final one. The triangle state is \((3,1,3,12,3,7)\). Hence its new parent has

\[
H=27\cdot2^k+27,\quad E(2)=108\cdot2^k-3,
\]

and the appended ready child has

\[
H(Q(TT_k))=H+E(2)-1=135\cdot2^k+23,
\quad \delta_y=1.                              \tag{2}
\]

This all-k identity upgrades to arbitrary M. In a fixed uniform module context every full count, external-root starting count, and external-root deletion count is a finite linear combination of \(K_j(M)\), the counts of ordered j nonempty paths partitioning M. On TT_k, \(K_j=j!S(k,j)\); these functions are linearly independent (their values for k=1,…,J form a triangular matrix with nonzero diagonal). Therefore equality for all transitive sizes identifies every coefficient, not merely a fitted scalar value. Since \(P_M(2)=2K_1(M)+K_2(M)\), equation (2) yields the universal identities

\[
H(Q(M))=23u+135v,\qquad \delta_y=u,
\qquad H(Q(M)-y)=27u+27v.                       \tag{3}
\]

The complete returned state is consequently

\[
(108v-4u,\;u,\;27v+27u).                       \tag{4}
\]

All these tournaments are strong: their singleton-module base tournaments are strong, and substituting any nonempty module into a strong tournament preserves strong connectivity. The receipt checks the base graphs directly. Equations (1)–(4) control the whole Hamiltonian count, including every sector omitted by a conditioned cycle-opening argument.

## 3. A labeled return with unchanged root star

At k=1 give the raw tournament's old parent vertices names 0,…,6, followed by z=7,x=8,y=9. The module is vertex 1. Give Q its construction order: triangle 0,1,2; module 3; zero-tail vertices 4,5,6; final one-tail vertex 7; x=8; y=9.

Use this exact map from raw names to target names:

| Raw name | 0 | M=1 | 2 | 3 | 4 | 5 | 6 | z=7 | x=8 | y=9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Target name | 0 | M=3 | 5 | 6 | 7 | 2 | 4 | 1 | 8 | 9 |

For general M match every internal module vertex with the same internal label. Pull Q back along this map. The resulting replacement changes precisely the six fixed arcs

\[
\{0,4\},\{0,7\},\{4,5\},\{4,7\},\{4,8\},\{5,8\},
\]

and all three uniform boundaries between M and old vertices 4,6,8. Every other arc, including every internal module arc and **every arc incident to y**, is retained. In particular y and its sole predecessor x remain the same physical vertices. There are exactly \(6+3|M|\) flipped arcs, independent of the internal structure of M.

This is an actual graph replacement with a fully specified retained-name map. It does not infer a return from matching scalar counts alone. A second realization found by the sibling search changes four fixed parent arcs and retargets the incoming root port, then migrates the module. It offers a different interface tradeoff; see [the return search](../returns/MODULE_MIGRATION.md).

## 4. Iteration boundary and next obstruction

Combining the growing-defect neutral switch and this successor replacement gives the uniform path

\[
\begin{array}{c}
H=45\cdot2^k+7,\quad\delta=6\cdot2^k+1
\\ \longleftrightarrow \\
H=45\cdot2^k+7,\quad\delta=1
\\ \xrightarrow{\text{two births}} \\
H=135\cdot2^k+23,\quad\delta=24\cdot2^k+5
\\ \xrightarrow{\text{named replacement}} \\
H=135\cdot2^k+23,\quad\delta=1.
\end{array}
\]

The second return is uniform in k and restores the precise +2-ready defect. It is **not yet an indefinitely repeatable grammar**: Q(M) has a different parent template from S(M). At the next two births the target is \(405v+71u\), and the defect becomes \(108v-3u\). No universal next replacement is proved here.

For u=1, the target \(135v+23\) cannot re-enter the original numerical family \(45v'+7\), even allowing an arbitrary integer v', by reduction modulo45. Nor does the next \(405v+71\) re-enter \(135v'+23\), by reduction modulo135. A successful iterative construction needs genuinely different templates or a changed module interface.

Allowing arbitrary M gives defect u rather than 1. Thus the universal module theorem does not manufacture a +2 offset for a module with H(M)>1. This is the exact information lost by replacing the module with its order alone.

## Verification

Run `python -B 05-knowledge/results/cycle_opening_six_return_20261112/compiler/module_return.py`. It independently counts complete paths on TT1,…,TT7 and three nontransitive modules, checks the physical root/port and every changed arc, and verifies both complete states. The largest literal graph has 16 vertices. The sibling `counts/cross_audit.py` independently checks both return realizations on transitive, cyclic, and five-core modules. The universal proof is coefficientwise and does not extrapolate these finite controls.
