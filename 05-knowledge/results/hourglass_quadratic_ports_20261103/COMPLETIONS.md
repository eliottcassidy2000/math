# Four-arc atoms and all-size strong-completion recurrences

**PROVED in the stated labelled completion model; FINITE-EXACT controls.**
The user's clarification concerns an oriented four-vertex graph obtained
by deleting two arcs of the strong tournament. It does not uniquely specify
which arcs. We classify all fifteen deletions, rather than identify a
particular drawing by its name. The seven-cell magic-square hourglass is
a different recovered object and is not a dependency here.

## H1. The complete small boundary

Fix the strong tournament with arcs

    0->1, 1->2, 2->3, 2->0, 3->0, 3->1.

Its fifteen two-arc deletions give nine oriented-graph isomorphism types.
The full labelled catalogue, multiplicities and completion counts are in
verification.json. There are two underlying graph shapes: a four-cycle
when the deleted edges are disjoint, and a triangle with one pendant edge
when they share a vertex. The three four-cycle cases are:

| Deleted arcs | Retained orientation | Strong completions of the two missing pairs |
|---|---|---:|
|0->1, 2->3|a source-to-sink path of length3 plus its direct shortcut|1|
|1->2, 3->0|two source-to-sink paths of length2|2|
|2->0, 3->1|a directed four-cycle|4|

For the middle, diamond-shaped orientation, the missing source/sink arc
must point back toward the source; the other missing arc is free. For the
directed cycle every completion remains strong. The first case has exactly
one successful choice, as inspection of the four completions shows.
Thus 1,2,4 count different guarded completion sets, not three independent
bits of a common orbit.

The middle case has a literal arithmetic source. From an ordered odd pair
(a,b), independently apply A:(a,b)->(a+2,b) and B:(a,b)->(a,b+2).
The four states x,(Ax),(Bx),ABx form exactly this diamond, because AB=BA.
The two intermediate states are distinct and have no intrinsic dominance
orientation. A strong tournament completion adds the return arc ABx->x
and one of the two middle-state orientations. That return is an added
operation, not the original arithmetic's forward evolution. On the
unordered-pair diagonal, folding a and b identifies the two middle states;
ports and the multiplicity of the two paths must then be retained.

## H2. The quadratic supplies a literal four-arc atom

The complete rational preperiodic set of x^2-29/16 is inherited from
THM-4139, `rational-three-cycle-order-six-lift-and-horizontal-carrier`,
with its chart limitations in THM-4146. On odd integer numerators put

    F(a)=(a^2-29)/4, B={-7,-5,-3,-1,1,3,5,7}.

The sign quotient a->|a| is a genuine dynamical factor, with map

    1->7->5->1, 3->5.

This is an oriented graph on four vertices: a directed triangle plus a
pendant incoming arc. Its missing pairs are {1,3} and {3,7}. Exactly three
of their four orientations make a strong tournament. The sole failure is
3->1 and3->7, making3 a source. All three successful completions have the
same strong four-vertex isomorphism class, but different named ports.
This type occurs in exactly three deletions of H1's fixed tournament.

The sign quotient has fibre sizes2,2,2,2. In contrast, the eight fixed-path
tilings have class fibres5,1,1,1, and converse merging gives5,2,1. Therefore
no bijection of the eight source states identifies these quotient maps.
Equality of their source cardinalities does not erase this obstruction.
The partial atom is an actual connection; a canonical tiling/class
identification is not supplied. NATIVE.md gives a separate chosen encoding
that does intertwine the allowed operations.

If the observer reads all future signed values F(a),F^2(a),..., the four
magnitudes are the exact minimal deterministic carrier: equal magnitudes
have equal next value, and the four next values are distinct. Observing
the current signed value instead needs all eight states. This is a scoped
four-state result, not minimality for all tournament observers.

## H3. Every rooted functional attachment has strong completions

Start with a directed cycle of length at least3, orient its remaining
pairs arbitrarily, and add vertices in parent-before-child order, prescribing
one outgoing arc from each child to its parent. At a stage with j existing
vertices, require at least one of the other j-1 old vertices to beat the
new child. Then the new tournament is strong: the old source is strong,
can enter the child, and is reachable from the child through its parent.
Conversely, a child with only outgoing edges cannot make this stage strong.

There are exactly 2^(j-1)-1 choices at this stage. Starting with the fixed
C3 gives exactly

    product_(j=3)^(n-1) (2^(j-1)-1)

labelled completions for which every insertion prefix is strong, for every
specified increasing rooted attachment forest. This counts a construction
subfamily. Later vertices can repair an earlier nonstrong prefix, so it is
not the total number of final strong completions. The first values3 and21
are construction counts, not forbidden Hamilton-path values.

## H4. Exact recurrence for every finite partial orientation

Let D be any loop-free partial orientation with at most one direction per
unordered pair, on a finite named set V. Define

    C_D(U)=2^(binom(|U|,2)-|D[U]|), C_D(empty)=1,
    S_D(empty)=0.

The exact number of strongly connected tournament completions on U is

    S_D(U)=C_D(U)
      - sum_(empty!=A proper subset U, no D arc U\A -> A)
              S_D(A) C_D(U\A).                         (1)

Proof: every nonstrong tournament has a unique source strongly connected
component A. All cross arcs must point from A to its complement, possible
exactly under the displayed guard. The inside of A has S_D(A) choices;
the outside has C_D(U\A); cross arcs have one choice. These cases are
disjoint and exhaust the nonstrong completions. Subtract from all completions.

This is an application of the inherited ordered-SCC decomposition, not
a new discovery of that decomposition. It gives an exact O(n^2 3^n)
arithmetic-operation algorithm with O(2^n) integer states, avoiding the
2^binom(n,2) full orientation census. Bit costs are separate. No assertion
of global algorithmic novelty or optimality is made.

For a functional graph consisting of one directed cycle and attached
in-trees, every proper permitted A avoids the cycle and is descendant-closed:
if a parent belongs to A, all its children must belong to A. Thus the
operation keeps an ancestry/ideal guard, not merely vertex and edge counts.
There may still be exponentially many such ideals.

## H5. Two uniform families and their exact generating quotients

Fix a named C3. Add m named vertices, each with its one prescribed outgoing
parent arc. Let B_m=2^(m(m+3)/2) be the number of all completions: there are
m+3 prescribed arcs among m+3 vertices.

Let u_k=2^binom(k,2) count all labelled tournaments, s_k their strong count,
and t_0=1, t_k=2^binom(k-1,2) for k>=1 count tournaments with a prescribed
Hamilton path. Let j_k be the strong subset of the latter. The standard
source-component decompositions give

    s_k=u_k-sum_(i=1)^(k-1) binom(k,i) s_i u_(k-i),
    j_k=t_k-sum_(i=1)^(k-1) j_i t_(k-i).                 (2)

For the fixed path, its SCCs are contiguous intervals in path order,
which explains the absence of a binomial coefficient in the second law.
Set s_0=j_0=0. Now:

* Leaves attached directly to any cycle vertices have total strong count

      L_m=B_m-sum_(k=1)^m binom(m,k) s_k B_(m-k).

  Every subset of leaves is an allowed source-component set, and no arcs
  among those leaves were prescribed. Parent choices on the cycle do not
  change this formula.
* A single chain of m vertices into one cycle vertex has total strong count

      Q_m=B_m-sum_(k=1)^m j_k B_(m-k).

  The allowed sets are the terminal suffixes of the chain; each carries a
  prescribed Hamilton path. There is one suffix at each size.

Equivalently, as formal series (no analytic convergence is claimed),

    sum L_m z^m/m! = (sum B_m z^m/m!)/(sum u_k z^k/k!),
    sum Q_m z^m    = (sum B_m z^m)/(sum t_k z^k).

The counts start

    m       0  1   2    3      4
    leaves  1  3  24  414  14280
    chain   1  3  28  479  15863.

This is an all-size distinction between unordered leaf selection and
ordered chain suffixes. At m=2 it supplies the smallest example in these
two families where final strong counts differ despite identical size,
arc count and stage-strong construction count21.

The bridge to the tournament catalogue is exact:

    s_k=sum_(strong classes T) k!/|Aut(T)|,
    j_k=sum_(strong classes T) H(T)/|Aut(T)|.

These use different weights on the same strong classes. They do not recover
individual H values or extend A000568 by themselves. The open compression
question is which descendant-ideal response data suffice under repeated
grafting; depth or a bare LCA tree shape without named ports is not assumed
to suffice.
