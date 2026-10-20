# Four-black normalization of every strong five-vertex end atom

Status: **PROVED**, using a complete small finite base and the uniform
known-arc argument below. This is an end-context theorem, not arbitrary
interior-SCC normalization.

Let A be any nonempty tournament with a source. Let B be any strong
tournament on five vertices. Then there is a path of four black
Hamiltonian-path reversals from A→B to a tournament isomorphic, by a
permutation of B alone, to

    A → 1 → 1 → C3.

A is restored as a literal adjacency matrix. There is no upper bound on
its size or any restriction on its other SCCs. Converse supplies the
corresponding left-end theorem with a sink context. When A itself is in
the three-atom language, [inventory.md](inventory.md) therefore joins the
source to the already connected golden component at its order.

## Inheritance and a necessary path choice

The preceding [R4 dissolution](../scc_basins_20261020d/arithmetic.md) supplies
the closest mechanism: even moves can restore arbitrary context arcs.
The new hostile is A=1→R4 in one five-vertex reduction: the standard
context path makes the second move blue. The correct theorem chooses a
different path in this context. Ignoring that choice would give a false
universal formula even though every path remains directed.

The least-used coordinate is a known source arc paired under positional
reversal with an untouched backward cross arc. Its asymmetry certifies
blackness in arbitrarily large contexts, even when endpoint scores alone
do not settle it.

## Six exact templates

Enumerate fixed-path five-vertex tournaments by forcing b0→b1→...→b4
and choosing the other six arcs. The 64 presentations include 50 strong
ones, representing exactly six ordinary isomorphism classes. The
independent count in [spectrum.md](spectrum.md) agrees. The following row
masks specify one representative of each class; bit j in row i means
bi→bj.

| index | five adjacency rows |
|---:|---|
|0|14,28,24,16,1|
|1|10,28,25,16,1|
|2|6,28,24,17,1|
|3|14,20,24,18,1|
|4|10,20,25,18,1|
|5|6,12,24,17,3|

The last is the regular C5. In the following table, a numeral i means bi;
P is a Hamiltonian word of A starting at its source, and P* is its
reverse. Each row lists the four successive reversal paths.

| index | first | second | third | fourth |
|---:|---|---|---|---|
|0|P,1,2,4,0,3|P*,4,2,3,0,1|3,2,1,0,4,P|2,3,4,0,1,P*|
|1|P,1,2,0,3,4|P*,4,3,0,2,1|0,1,2,3,4,P|4,0,3,2,1,P*|
|2|P,0,1,2,3,4|2,1,0,P*,4,3|0,1,2,4,P,3|4,0,2,1,3,P*|
|3|P,0,2,3,1,4|1,3,2,0,P*,4|3,1,2,4,P,0|2,1,3,4,0,P*|
|4|P,0,1,2,3,4|1,0,P*,4,3,2|0,1,4,P,2,3|4,1,0,3,2,P*|
|5|P,0,1,2,3,4|0,P*,4,3,2,1|3,0,2,4,P,1|2,0,3,4,1,P*|

Every row is an identity of boundary-port path words. Every internal
path edge of A is toggled four times; no other internal A edge is
touched. Every unwanted cross-boundary edge is toggled twice. The
remaining finite local edge ledger produces two singleton SCCs followed
by a triangle, with A preceding all of them. This ledger uses only the
path arcs of A and the original forward cross arcs. Thus directedness
and the literal target hold for every A and P, regardless of source or
blackness; the m=1 port identification is checked separately.

For example, row 5 gives exactly

    A → b1 → b4 → C3(b0,b2,b3).

Its cross-boundary cancellations are the pairs {last(P),b0},
{first(P),b4}, and {last(P),b1}. An independent complete ledger and
checker are in [five_atom_audit.md](five_atom_audit.md).
The other exact labelled local targets are retained in
`five_atom_probe.json` and checked by `five_atoms_verify.py`.

## Uniform blackness of the two middle moves

Write P=(a0,...,a_(m-1)) and u=a0. The first and final tournaments both
have a singleton as their first SCC, and respectively a five-vertex
strong SCC and a triangle as their last SCC. Neither is self-converse,
so the first and fourth paths are black. Colour is unchanged when a
path-reversal edge is traversed backwards.

For m≥7 the following table proves blackness of the middle paths.
The ordered pair in its left column is always an arc at that step;
the reversed positional image in its right column is always a nonarc.
A blue path would equate them. The inequalities in the indices ensure
that these are distinct context vertices and none of the cross arcs
listed is an altered boundary edge.

| template / step | arc | reversed positional image, a nonarc |
|---|---|---|
|0 or 1, step 2 or 3|u→a_(m-5)|b4→a_(m-6)|
|2, step 2|u→b3|b2→b0|
|2, step 3|u→a_(m-3)|b4→a_(m-4)|
|3, step 2|u→b1|b4→b3|
|3, step 3|u→a_(m-3)|b4→a_(m-4)|
|4, step 2 or 3|u→a_(m-1)|b4→a_(m-2)|
|5, step 2 or 3|u→a_(m-3)|b4→a_(m-4)|

For the context pairs, the source edge is not a consecutive path edge
in the first step, or has been restored by the second step. The cross
edge still points from the context to B. For the two entirely local
image pairs, inspection of the specified five rows and first path gives
the displayed values. Computing positional reversal in the four-word
table gives exactly the listed images. This establishes all m≥7 without
enumerating context tournaments.

## The small base and its sole repair

Relabel a chosen path of A as (0,...,m-1). With vertex 0 a source, every
nonconsecutive pair among vertices 1,...,m-1 is free. The complete
universe at m=1,...,6 therefore has sizes

    1,1,1,2,8,64.

All six templates are directed and have the required literal target
throughout this universe. All are black except template 0 on

    A=(30,12,24,16,2),   P=(0,1,2,3,4),

where only step 2 is blue. This is the explicit source followed by R4.
Replace P by

    (0,2,3,4,1).

The same template then has four black moves and the same literal A in
its target. This completes the small base and the theorem.

The proof uses an explicit finite classification of six atoms and
77 path presentations for the small context base. It does not extrapolate
from a finite absence of blue paths: the unbounded part is the fixed-arc
contradiction above. `five_atom_lift.py` additionally tests all 1,101
source/path presentations through m=7 and preserves the failing
unrepaired witness. `five_atoms_verify.py` checks the repair and the
parametric certificate table.

## Exact boundary of the consequence

This normalizes every five-vertex strong atom at an end with the stated
context. It does not supply a move for A→B→W with an arbitrary nonempty
right context W, nor does it permit every larger SCC to be replaced by
small atoms. These are the remaining operations needed for full terminal
normalization. The new theorem and the independent outer-blue detours
address different obligations: normalizing a vertex and replacing a
specified blue edge are not equivalent tasks.
