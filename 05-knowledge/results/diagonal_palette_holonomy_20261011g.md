# Diagonal three-phases, reflection and the missing seam data

2026-10-11. **PROVED:** the nine-colour lower bound under the stated global
symmetry assumptions, its sharp square-tiling control, and the marked seam
criterion. **FINITE-EXACT:** permutation and signed-seam controls.
**OPEN:** an unrestricted six-colouring of the unit-distance plane. The
symmetry obstruction does not apply to a construction that breaks its rules
on the diagonal or uses position-dependent palette actions.

## 1. Inheritance and scope

The closest proved interface is the binary gluing theorem in
[pants boundary holonomy](pants_boundary_holonomy_20261011b.md), together with
[the three bipartite parts formulation](plane_cycles_pants_and_six_seven_20261011b.md).
The older [glued chessboards](glued_chessboard_rings_scaffolds_20261006.md)
and [chessboard weave](chessboard_weave_20261006.md) retain seam identifications
and show why equal local counts need not preserve global colourings.
The current task's diagonal is a geometric reflection axis, not merely an
unoriented edge in a quotient picture.

Closest mechanism: a marked parity transport. Hostile: individually bipartite
half-plane pieces can glue to a unit triangle. Corrected near miss: reversing
a colour cycle is not the same action as reflecting across its direction of
translation. Least-used sidecar: the stabilizer of each actual point, including
points on the diagonal. The board is **three-phase orbit / colour-pair bit /
reflection / fixed axis / complete crossing edges / component frame**.

This applies the META-PATTERNS cards 'Preserve the selected side, not only the
walls' and 'Retain local incidence through counting'. The anchor is the
six-colour plane construction, the niche its symmetry algebra, and the wildcard
the sharp nine-colour model. No new general meta-pattern is claimed.

## 2. Translation and reflection have to satisfy the same relations on colours

Let R(x,y)=(y,x) and T(x,y)=(x+h,y+h), h nonzero. Geometrically RT=TR.
Suppose a colouring c has fixed palette permutations r,s satisfying

\[
 c(Tp)=r(c(p)),\qquad c(Rp)=s(c(p))\quad\text{for every point }p.
\]

On the used colours, necessarily sr=rs and s squared is the identity.
Moreover every diagonal point has a colour fixed by s. A point at distance
1/2 from the diagonal and its reflected point are unit separated, so a proper
unit colouring also uses at least one colour moved by s.

Write six colours as (j,b), j in Z/3 and b in Z/2. Advancing all three pairs
in the same direction means r(j,b)=(j+1,b). The tempting reflection rule
s(j,b)=(-j,1-b) reverses that advance and does not commute with it. This rule
belongs to a dihedral relation, whereas these two geometric maps commute.
A half-turn about a point on the diagonal reverses the translation and is a
different possible geometric symmetry; it must not be confused with R.

Retaining the side b and using opposite advances,
r(j,b)=(j+(-1)^b,b), does restore commutation with that s. But s has no fixed
colour, so no point of the diagonal can obey the same rule. A special seam
rule is therefore a mathematical change of the construction, not an ignorable
measure-zero convention.

## 3. Sharp obstruction: this global three-phase symmetry needs nine colours

**Theorem.** Suppose a proper colouring of the whole Euclidean unit-distance
plane has the above fixed translation and reflection actions, and every used
colour lies in an r-orbit of length exactly three. Then it uses at least nine
colours. This is sharp within the stated symmetry class.

**Proof.** A colour appearing on the fixed diagonal is fixed by s. Commutation
then makes its entire three-element r-orbit fixed by s. There is also a used
colour moved by s, from a reflected unit-separated pair. Its r-orbit cannot be
preserved setwise by s: a permutation commuting with a transitive three-cycle
on that orbit is a power of that cycle, and the only such involution is the
identity. Thus s exchanges this orbit with a distinct three-element orbit.
Those six colours are disjoint from the fixed orbit. Total: at least nine.
This proof uses every point of the plane, with no measurability assumption.

**Sharp positive control.** Let a=3/5 and colour the half-open square
[ai,a(i+1)) times [aj,a(j+1)) by (i mod3,j mod3). There are nine colours.
Two points in the same square have distance less than a sqrt(2)<1. Different
squares of the same colour have some index difference of absolute value at
least three, hence a coordinate gap at least2a>1. The colouring is proper.

Translation by (a,a) acts as (i,j)->(i+1,j+1), with three orbits of length
three. Reflection exchanges i,j, fixes the three diagonal colours and swaps
the other two phase orbits. Half-open boundary ownership is unchanged by
coordinate exchange, so both symmetry formulas hold on boundaries as well.

At six colours, direct enumeration independently gives18 permutations
commuting with two disjoint three-cycles. Only four are involutions: identity
(six fixed colours) and three with no fixed colours. Identity fails on a
reflected unit pair; the other three fail on the fixed axis. At nine colours
there are162 commuting permutations and ten involutions; the nine nonidentity
ones fix three colours, matching the construction.

This is **not** a lower bound of nine for the unconstrained plane. In particular,
changing translation's cycle type to3+1+1+1 admits a commuting reflection that
fixes the three-cycle and one other colour while exchanging the remaining two.
That removes this symmetry contradiction only; it supplies no proper metric
colouring. Position-dependent rules or a separately coloured seam are other
ways to leave the theorem's hypotheses.

## 4. What must replace a simple matching across the diagonal

A six-colouring can be written as a three-valued pair label j and a bit b.
For each fixed j, its induced unit graph must be bipartite; this is the
inherited exact three-parts reformulation. Cut the plane into its two open
half-planes and the axis, and suppose actual binary colourings beta_C are
supplied on every connected component C of each pair-labelled piece.

Each component admits an independent frame flip x_C. For every actual unit
edge pq crossing pieces with equal pair label, impose

\[
 x_{C(p)}+x_{C(q)}=1+\beta_{C(p)}(p)+\beta_{C(q)}(q)\pmod2.       (1)
\]

Keep parallel constraints, component identities and axis vertices. A global
bit assignment exists iff every finite cycle of this labelled component
multigraph has zero sum. Necessity follows by summing (1). For sufficiency,
choose a root in each component and transport a bit along a finite path; the
cycle condition gives path independence. For arbitrary infinite graphs this
is an ordinary ZFC existence argument, not an effective or measurable selector.
For finite input it is an executable check and exports the actual flips.

Here is a geometric minimal hostile in orthonormal coordinates (u,v), where
v=0 is the diagonal. Put

\[
 A=(-1/2,-\sqrt3/6),\quad B=(1/2,-\sqrt3/6),
 \quad C=(0,\sqrt3/3).
\]

All three pairwise distances are one. The lower piece is the edge AB, and
the upper piece is isolated C, each bipartite. Give A,B local bits0,1 and C
bit0, all with the same pair label. The two crossing edges demand both
x_lower+x_upper=1 and x_lower+x_upper=0. Pairing mirror points or checking one
edge per component pair would miss this obstruction. It is detected by a
two-edge cycle of opposite seam labels.

The promising generalization of the user's diagonal is therefore a labelled
constraint boundary, not just a bijection between point sets. The preserved
object is the complete same-pair edge parity. Its finite obstruction is an odd
cycle, while a successful positive construction must supply the pair labels
and satisfy all such constraints, including those on the diagonal itself.

## 5. Reproduction and limits

Run [diagonal_palette_holonomy_20261011g.py](../../04-computation/experiments/diagonal_palette_holonomy_20261011g.py)
with ordinary Python and python -O. The declared finite universe consists of
all720 six-colour permutations; the162 parametrized commuting nine-colour
permutations; rational test points including negative square boundaries;
and all signed seam systems on at most three vertices, with loops and both
possible labels independently present on each vertex pair. Literal binary
assignments independently check the seam solver. Exact-type hostile inputs
are rejected. The all-plane assertions are the proofs above, not an inference
from these finite test sets.

Normal and optimized output agree at **9,210 always-active checks**. The saved
[output](diagonal_palette_holonomy_20261011g.out) has LF-normalized SHA256
`37c641724f65029f2ca1153ae24fa12b83e2fd025c977acc699c68a45f75cb9a`.
