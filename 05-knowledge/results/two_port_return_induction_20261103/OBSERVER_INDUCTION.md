# Two-port inverse faces and an exact marked-source compiler

**PROVED:** an explicit four-face necessary condition for adaptive outputs;
the positive C3 observer ray has only its original point in the strong-source
image; and a terminal-port liveness criterion compiles any compatible named
ancestry into an actual two-port program. **EXACT FINITE CERTIFICATES:**
command ranks and image equations, and the stated small hostile controls.
All-target prime/7p induction remains **OPEN**.

The main positive addition is an exact source-entry test when an actual
chronological ancestry is supplied. The negative addition is that even a
strong child with an actual integral graph predecessor can fail entry into
the C3-seed language. Inverse matrix equations need both strong-parent and
terminal-port guards.

## 1. Inheritance and the live board

The closest proved mechanisms are
[LIVE_PORTS](../adaptive_state_induction_20261103/LIVE_PORTS.md), its 24-coordinate
positive carrier and sharp 17-query observer, and
[COMMON_OBSERVER](../adaptive_state_induction_20261103/COMMON_OBSERVER.md),
which extends the common response quotient to all actual two-marked sources.
The inherited 999 coverage is not rerun or extended here.

The canonical hostile is the order-six equal-rooted-profile pair from
[ADAPTIVE](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md), whose
older-root responses differ. The corrected near misses are identifying
positive integral observer vectors with actual sources, and applying the
unmarked future-live-width criterion while forgetting the prescribed final
ports. The least-used sidecar is the chronological lifetime of a **name
required at the end**, even when that name is never again used as a root.

| Live concept | Hypothesis / cheap probe | Outcome |
|---|---|---|
| Inverse images | Six singular command maps might have useful finite image tests | Four rank-13 spaces cover all last-command images; their semantic faces are explicit. |
| Actual strong cone | Positive integral, odd query data might certify entry | The actual weak H9 scalar-ray state has no strong representative. |
| Physical predecessor | An actual predecessor plus strong child might suffice | A marked strong H5 child has a weak predecessor but no C3-seed program. |
| Terminal liveness | Include final ports in the old future-root width | A two-color lifetime test gives an exact constructive criterion. |
| Prime/7p induction | Invert an observer and compile the source | Compilation is solved given ancestry; choosing realizable predecessors is still open. |

The explored representation changes were `query vectors → rooted profiles`,
`linear image → actual strong image`, and `root word → terminal-marked
lifetimes`. No large scalar-support census or unconstrained formal-state
search is used.

## 2. Command images: four large branches, not unique inverses

Use the inherited six commands a=3r+d, with root slot r∈{0,1} and discarded
slot d∈{0,1,2}; d=2 forgets the newborn. Let A_a be their exact 17-by-17
rational query matrices, acting on row vectors. Exact row reduction gives

    rank(A_0,...,A_5)=(13,4,13,4,13,13),
    Row(A_1)⊂Row(A_4),   Row(A_3)⊂Row(A_0).           (2.1)

Thus every nonempty program's output is in

    Row(A_0) ∪ Row(A_2) ∪ Row(A_4) ∪ Row(A_5).        (2.2)

`observer_induction_probe.json` supplies a primitive integer basis of the
right kernel of each matrix. The corresponding equations Yz=0 exactly
test membership in its rational row image. They are coefficient identities,
not fits to sampled source vectors. For a rank-13 command, a rational
preimage, if it exists, is an affine four-dimensional set; a reset command
has a thirteen-dimensional preimage fibre. All six maps are singular.

The containments in (2.1) are **linear image** containments. They do not
authorize replacing command 1 by command 4, or command 3 by command 0, on
an actual source. An inverse must still choose an integral physical state
and prove its realization by a compatible program.

The adaptive query basis remains valid on arbitrary actual two-marked
sources: its 17 columns are independent on the old 24-coordinate carrier,
while the common universal 25-coordinate continuation space has dimension
17. They are therefore a basis of that universal space too. The six query
matrices agree there by uniqueness of coordinates. This is the inherited
common-observer theorem, not an assumption that every source arose from C3.

The adaptive and common static query bases also have the same **integral
content** on actual sources. The exact full-carrier certificate
`observer_induction_basis.json` gives

    X_static = X_adaptive B,   det(B)=−1/4,

with denominator least common multiple two for both B and B⁻¹. The selected
17-row carrier minors have determinants −16384 and 4096 respectively; both
changes of coordinates are checked on all 25 carrier rows, not just the
selected minor. Every actual query is an odd Hamiltonian-path count, so
each vector has odd gcd. From 2B and 2B⁻¹ integral, each gcd divides every
coordinate of the other vector. Thus the two query gcds are equal for every
actual two-marked tournament. This justifies using the same continuation
content in the adaptive and common-static conventions; it does not equate
their formal nonnegative cones or certify source realization.

## 3. A semantic four-face invariant

For a retained vertex v in an actual source T write

    h=H(T), s_v=# Hamiltonian paths starting at v,
    q_v=H(T−v), F_v=# two-slot covers (P,Q), P ending at v.

P is nonempty; Q may be empty. These four quantities can be recovered from
two-port future queries. Let K_r=(r,discard newborn), and let E_r be the
command that retains the chosen old root and the newborn, replacing the
other port. Denote counts after these actual commands by H_1,H_RR,H_RN,
where RR means K_r followed by K_r, and RN means E_r followed by a birth at
the retained newborn. Then

    F_r=(H_RN−h)/2,
    q_r=F_r−(H_RR−H_1)/2,
    s_r=h+F_r−H_1.                                  (3.1)

The inherited insertion identities prove these equations at every order.
In encoded words, RN is (1,5) for slot 0 and (3,2) for slot 1. Thus these
are homogeneous linear observations in the seventeen-coordinate space.

Define the nonnegative physical quantities

    δ_r=F_r−h−q_r+s_r,    η_r=q_r−s_r.

The exact endpoint-pair identity in prior ADAPTIVE §3 gives
δ_r=0 iff indeg(v_r)≤1, and for a strong source this means indegree one.
The inequality η_r≥0 follows by deleting the first vertex of a path
starting at v_r.

**Birth-face theorem.** Every C3-seed two-port state satisfies

    δ_0 δ_1 η_0 η_1=0.                              (3.2)

If the last newborn is retained, that port has indegree one, so its δ is
zero. If the newborn is forgotten, the chosen old root has s'=q' by the
fixed-root update, so its η is zero. The initial C3 already has both δ and
both η zero. More specifically, commands 0,4 force δ_0,δ_1 respectively;
commands 2,5 force η_0,η_1; command 1 forces η_0=δ_1=0, and command 3
forces δ_0=η_1=0.

This is a necessary union-of-faces condition, not a realization theorem.
For example, the regular R5 with outgoing rows (6,12,24,17,3), retaining
ports (0,1), has H=15 and at each port

    (h,s,q,F)=(15,3,5,21).

Its four factors are (4,4,2,2), so its observer is outside the two-port
language. The scalar 15 itself is attained there. The stronger exact
row-image tests also reject this source in every last-command branch.

## 4. The strong-source image is not a positive cone

**Lemma.** If a strong tournament of order greater than three has two
indegree-one vertices x→y, then T−y is strong.

Every simple path ending at x avoids y, since entering y would first require
visiting x. Thus every old vertex reaches x in T−y. Let u be x's sole
in-neighbor. Choose z distinct from x,u,y, possible at this order. A path
from z to x avoiding y ends with u→x; before that edge it gives a path
z→u avoiding x. Hence x→z→...→u. Since x points directly to every other
vertex of T−y, it reaches them all. This proves strongness.

In this situation η_x=0 and η_y>0. Deleting x makes y a universal source,
so every Hamilton path of T−x starts at y and can be prefixed by x; hence
q_x=s_x. Strongness of T−y supplies a Hamilton path starting at x, giving
q_y>s_y. Thus when both retained ports have δ=0 in a strong source of
order greater than three, their rooted profiles identify which one is the
legal downstream newborn. The order-three triangle, where both η vanish,
is the necessary boundary case.

**Strong-ray theorem.** Let X_Δ be the seventeen-query vector of the actual
two-marked C3 seed. For a positive integer k, the state kX_Δ has a strong
two-marked tournament representative if and only if k=1.

By (3.1), such a representative would have rooted profiles
(3k,k,k,3k) at both distinct ports. Both δ vanish, so both ports have
indegree one. Both η vanish too. In a strong tournament of order greater
than three, the lemma makes the later indegree-one vertex y legally
deletable. But q_y−s_y counts Hamiltonian paths of T−y starting at x;
there is at least one because T−y is strong. This contradicts η_y=0.
At order three, the strong tournament is C3 and H=3 forces k=1.
The original seed proves the converse.

The actual weak image behaves differently. For any nonempty tournament R,
place R after the seed in an ordinal sum C3→R, retaining the ports inside
C3. Every future newborn dominates all of R, since no vertex of R is live.
After any program w the actual graph is T_w→R. Consequently

    X(C3→R)=H(R) X_Δ.                               (4.1)

This gives actual integral, positive, odd query vectors along the ray,
not merely formal multiples. With R=C3, the source rows are

    (58,60,57,16,32,8), ports (0,1), H=9.

Its complete query vector is exactly 3X_Δ, but the theorem excludes every
strong representative of that observer. A different actual strong H9
source, rows (26,4,1,6,14) from two fixed-root births, has another query
vector. The scaled seed ray passes the rational last-image tests for
commands 0,1,4,5, so those tests do not remove this strong-source obstruction.

Strongness is **not** determined by the observer in general: taking R to
be a singleton in (4.1) gives the same X_Δ from a weak source. The theorem
classifies the strong representatives along this particular ray; it does
not recover graph connectivity from every query vector.

## 5. A strong marked H5 outsider with an actual weak predecessor

Take T4=birth(C3,root 0), with rows (10,4,1,6), but prescribe final ports
(1,2). Its first three observations are

    (H,H_after_root_slot0,H_after_root_slot1)=(5,13,11).

The six one-birth programs from the prescribed oriented C3 seed have only
the following triples:

    (5,9,13), (5,9,9), (5,11,9).

Every later adaptive birth strictly increases H. Thus no C3-seed two-port
program has the prescribed H5 observer. The source graph itself is strong
and is an adaptive descendant; the obstruction is its retained pair.

It nevertheless satisfies (3.2), with factors (2,2,2,0), and lies in the
exact rational image of command 5. This is not solely a formal inverse:
the transitive parent on names (1,2,3), rows (2,0,3) in that order, retaining
its slots (0,1), has H=1. Apply command 5 and identify its newborn with
physical name 0. The name transport (1,2,3,0) gives exactly T4 with ports
(1,2). Thus an actual nonnegative integral graph predecessor and a strong
child still fail to supply a **strong C3-seed predecessor program**.

The lost coordinate is visible before the birth: one would have to retain
root 0 as well as both terminal names 1 and 2. Two live slots cannot retain
those three names simultaneously. This leads to an exact repair.

## 6. Terminal-marked liveness is exactly compilable

Suppose an actual C3 ancestry has been supplied, with chronological roots
r_t for births t=3,...,n−1 and ordered final ports p_0,p_1. Retain the
initial C3 orientation 0→1→2→0. The three allowed initial pair gauges are
(0,1),(1,2),(2,0), with their cyclic name transports recorded.

At each state before birth t, and at the final state t=n, define

    L_t=({r_i:i≥t} ∪ {p_0,p_1}) ∩ {0,...,t−1}.       (6.1)

Each name's mandatory states form an interval: from its availability to
its last root use, or through the final state if it is a final port. Form
the overlap graph on mandatory names, joining two names when they occur
together in some L_t.

**Terminal-port compiler theorem.** The named ancestry can be produced
with the prescribed ordered final ports if and only if:

1. |L_t|≤2 for every t;
2. the overlap graph admits a proper coloring by slots 0,1, with
   color(p_i)=i;
3. the coloring of L_3 extends to one of the three allowed oriented initial
   pair gauges.

The orientation guard matters even for the empty program: the reversed
ordered seed pair (1,0) has width two but is not one of those initial gauges.
For T4 with final (1,2), L_3={0,1,2}, so the first condition already fails.

**Necessity.** Every mandatory name must remain retained from its birth
until its final use, since a discarded name can never be selected again.
Its slot is unchanged during that lifetime. These slots give the coloring,
the width bound, and the initial and terminal guards.

**Sufficiency and construction.** Choose a compatible coloring and initial
pair. Whenever a newborn has a future mandatory lifetime, retain it in its
assigned slot by discarding that slot's old occupant. The old occupant
cannot still be mandatory after the birth: it would overlap the newborn
with the same color. When the newborn is never mandatory, discard it.
Keep arbitrary old fillers in unneeded slots. Induction gives every L_t
in its prescribed slots, so each required root is available and the final
slots are exactly (p_0,p_1). This constructs the command word and all names.

The test and compiler take O(n²) time using the direct construction of the
sets (6.1); no Hamiltonian counts are required. For a supplied tournament,
first use the inherited unique strong ancestry to decide whether its core
is C3 and recover the root word. The theorem then decides the actual
terminal-marked source image. The O(n²) claim concerns the given root word,
not the separate graph-ancestry computation.

This refines the inherited unmarked width theorem. It does not decide
whether a bare seventeen-coordinate vector has *some* compatible actual
source; the graph or ancestry is a necessary input to this compiler.

## 7. The constructive induction reformulated

The inverse problem now has four distinct checks:

* Choose a last-command branch using exact image equations and the
  nonnegative rooted faces; solve its rational affine preimage equations.
* Keep integrality, nonnegative physical query constraints, and a smaller
  parent H when working within strong adaptive growth.
* Supply an actual strong predecessor with C3 ancestry. The ray and H5
  examples show that matrix compatibility and graph existence alone do
  not certify this step.
* Apply the terminal-port compiler to produce the required actual program
  and its named ports, then perform the selected command.

The last check is now an exact constructive theorem, while the predecessor
entry step is still open. For prime/7p induction one may need intermediate
parent counts outside those lanes. No selection of formal nonnegative
query vectors, no inverse of a singular command on its rational image,
and no isolated actual repair proves all-target entry.

## 8. Reproduction and scope

`observer_induction_probe.py` pins the prior 17-query matrix JSON by hash.
Exact fractions verify command ranks, primitive image equations, and both
row-image containments. The graph checks use independently constructed
adjacencies, subset Hamilton DP, and literal vertex-word rooted profiles.
They verify all seventeen queries for the four hostile/control sources,
the actual weak H5 predecessor and its name transport, and the six actual
one-birth response triples.

The liveness algorithm is independently compared with exhaustive retained-slot
transition sets for all root histories of zero through three births and
all ordered final pairs: 2,082 cases, with 549 accepted programs replayed
as exactly the same named graphs and ports. This is a source-compilation
audit, not a target-count coverage scan.

```sh
python3 -I -B 05-knowledge/results/two_port_return_induction_20261103/observer_induction_probe.py
python3 -I -B -O 05-knowledge/results/two_port_return_induction_20261103/observer_induction_probe.py
python3 -I -B 05-knowledge/results/two_port_return_induction_20261103/observer_induction_basis.py
python3 -I -B -O 05-knowledge/results/two_port_return_induction_20261103/observer_induction_basis.py
```

Each normal/optimized pair agrees with its corresponding JSON receipt.
The basis probe verifies 850 full-carrier coefficient identities. All-order claims
use the stated structural proofs and inherited exact observer identities;
the finite controls audit implementation and preserve the small failure
boundaries.

An independent agent proof review passed §§4–6. Using a separate graph/query
implementation, it also verified all seventeen weak C3→C3 ray coordinates,
the six one-birth response triples, the marked strong H5 outsider, and the
literal weak-parent name transport. That review did not independently rerun
the 2,082-case compiler universe; the exhaustive slot-transition comparison
above is the producer's independent algorithmic path.
