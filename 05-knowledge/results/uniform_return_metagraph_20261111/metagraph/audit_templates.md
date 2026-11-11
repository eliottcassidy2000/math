# Independent audit of the fixed-defect inflation theorem

**PROOF AUDIT: PASS.** This note independently reviews
[BOUNDED_TEMPLATES](../defect/BOUNDED_TEMPLATES.md), including its exceptional
vertex bound, positive inflation injection, uniform box threshold and
isomorphism quotient. It adds the automorphism refinements below. It does
not claim a bounded-defect source for every Hamilton count.

## 1. The exceptional-vertex argument counts distinct positive objects

Use the conventions of the reviewed proof: the root is r, its inports
are I, d=|I|>=2, and the ordered residual SCC chain is R. Contact letters
record the inports beaten by a residual vertex. The all-one letter is
ordinary. The positive defect is the sum of I-starting Hamilton paths
and the pair capacities c_ac. Each tagged odd cycle occurring alone has
weight two in the appropriate capacity.

For d>=3 a mixed-contact singleton x gives a four-vertex path with
inport endpoints and x as its only residual vertex. The choice of two
inports on one side and one on the other is possible precisely because
d>=3. The directed edge between the two inports determines which of the
two valid path orders to use. Closing with the two-edge tagged ear gives
a five-cycle. Different x have different supports.

A zero-contact singleton x cannot leave toward an earlier residual SCC,
an inport, or r. Strongness therefore gives a later residual vertex y
beating some inport c. For any other inport a, a,x,y,c is a directed path.
The resulting cycle contains exactly one zero-contact residual vertex,
namely x: y beats c. Thus cycles selected for different zero-contact
vertices are distinct even when their chosen y is shared. They are also
distinct from the preceding one-residual-vertex cycles.

For d2, two ordered residual singletons with the same mixed contact give
the corresponding four-vertex path between the two inports. Distinct
pairs have distinct supports; the two contact types are disjoint; and
these cycles contain no zero-contact residual vertex. Consequently the
reviewed bounds add distinct atoms in one positive sum. They do not
require disjointness of cycles within one packing. This distinction is
essential and is correctly retained in the proof.

## 2. Nonempty transitive inflation is monotone for the actual defect

Fix the finite skeleton and the set of nonempty ordinary slots. In each
slot, the ordinary vertices form a transitive module. Insert one new
smallest clone immediately before the old smallest clone in its path.
Any predecessor of the old smallest clone is outside the module: a
later clone cannot point to it. External incidence is uniform, so the
new path is valid. If the old clone begins its path, the insertion simply
changes that unprescribed start. No prescribed start or end can be that
clone, because all prescribed vertices lie in I. The argument applies
equally to an I-starting Hamilton path and to each ordered two-path
cover defining W_ac. Deleting the new clone recovers the original object
without changing which ordered bin contained it.

This proves monotonicity separately for every positive term, hence for
the exact defect b. There is no subtraction or appeal to a difference of
two monotone quantities. The distinction between empty and nonempty slots
is retained. Strongness is invariant under these nonempty module
inflations because the unchanged quotient contains r and its inports.

## 3. The threshold is uniform over the other coordinates

Let a deterministic finite automaton recognizing the fixed skeleton's
exact-defect strong sources have m states. For any prefix state and any
fixed suffix, membership as one ordinary block length varies is periodic
from length m onward. On positive lengths it is also an interval:
strongness is constant, while the actual integer defect is nondecreasing.
An interval which is periodic on that tail has constant membership there.

The same m works for every prefix and every suffix. Thus changing any
coordinate from a value >=m to m preserves acceptance while all other
coordinates are held arbitrary. Applying this successively proves the
finite-grid/Cartesian-box assertion; there is no hidden dependence of
one coordinate's threshold on the sizes of other blocks.

An independent explanation is available directly in the positive run
automaton. Under repeated identical all-one letters, each run coordinate
can only move to the one fixed last-vertex state for that letter. Its raw
transition graph has no directed cycles except self-loops. Sufficiently
long walks can be injected into one-step-longer walks by inserting a
canonical self-loop at the first visited state with a loop. Raw register
counts are eventually nondecreasing, and their bounded capped values
eventually stabilize. The monotonicity-plus-DFA argument above is shorter
and does not require unpacking those registers.

## 4. Automorphisms are constant on each nonempty-slot skeleton

This is a consequence of the actual rooted structure, not an assumption
that the retained ports are individually named. A rooted automorphism
preserves I setwise and preserves R. The strong components of R have an
intrinsic linear order, so every component position is fixed. Every
singleton in an ordinary run is fixed individually; a finite transitive
tournament has trivial automorphism group. No two retained residual
components can be exchanged, even if they are isomorphic.

The remaining possible automorphisms are exactly the finite simultaneous
permutations of I and the retained nontrivial residual cores which
preserve the finite contact data. Their conditions do not depend on any
positive ordinary-slot length. Conversely each such finite automorphism
extends uniquely by fixing every ordinary vertex. Therefore both the
abstract rooted automorphism group and its action on the retained
skeleton are constant throughout each box.

Consequences, all for fixed branching defect b:

* Rooted classes of any specified rooted automorphism-group isomorphism
  type have a generating function with denominator a power of 1-x, and
  their order counts are eventually polynomial.
* The same holds after weighting each rooted class by 1/|Aut(T,r)|.
  Thus the labeled count with a distinguished root, divided by n!, is
  eventually polynomial. With a prescribed root label, divide instead
  by (n-1)!. The factorial itself is not a polynomial growth claim.
* Fixed numbers of additional distinct named vertices give the same kind
  of eventual polynomial enumeration. Burnside counts each named tuple
  by the falling factorial of the number of fixed vertices; inflation
  merely adds individually fixed vertices. Summing a fixed-degree
  polynomial over a Cartesian box preserves the stated denominator form.
* The full cycle-index contribution of a skeleton automorphism acquires
  only additional p_1 factors as the ordinary blocks grow. Nontrivial
  cycle inventory remains finite and fixed.

Hamilton-path-frame weights are different: H can grow exponentially even
inside a single fixed-defect box. No eventual-polynomial frame count is
claimed. Forgetting the root can also identify different eligible roots;
unrooted enumeration requires a separate root-forgetting audit.

If a skeleton has e exceptional singletons and c nontrivial residual
components, it has at most e+c+1 ordinary slots. Its ordinary rooted
class count therefore has eventual polynomial degree at most e+c.
The exceptional bounds and the residual H bound make this an effective
bound depending only on b. It is generally far from sharp.

## 5. Source compatibility and stopping boundary

The classification retains an actual rooted graph and all contact data;
it does not collapse sources having equal H or equal scalar defect.
Each free coordinate is a specified physical transitive-module inflation
at a specified residual position. Thus the boxes are valid infinite
families of source interfaces. They do not imply that a desired outgoing
move preserves a box, preserves b, or returns to any ready interface.

The existing complete b2, b4, b5 and b6 classifications are concrete
controls: eventual constants, a finite stratum, and the b6 linear law
5n+2 all have the predicted form. The proof is the positive-object and
finite-state argument above, not extrapolation from these controls.
