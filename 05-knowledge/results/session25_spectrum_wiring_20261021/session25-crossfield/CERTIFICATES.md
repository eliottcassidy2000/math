# Guarded acceleration and the separator-width obstruction

2026-10-21. PROVED elementary scoped results and FINITE-EXACT controls.
No Collatz convergence, uniform153N estimate, six-colour equality gadget,
or value of the plane's chromatic number is established. No priority claim.

## Inheritance

Read-only base f002fc425336. Closest mechanisms: the exact two-port count
algebra and serial/parallel no-go in
`05-knowledge/results/existence_count_amplifiers_20261011e.md`; the true
source-owned clock/credit interface in
`05-knowledge/results/residual_flow_fibre_clocks_20261020k.md` and
`linked_tail_geometry_and_six_center_orders_20261021c.md`.
The earlier `arithmetic_seams_20260921_dynamics.md` constructs an actual
S2 x S3 action on six vectors, with B^3=-I and JBJ=B^-1. Its reflection
reverses time; it is not a commuting symmetry of a directed Collatz map.
Canonical hostiles: serial inequality gates at palette two force equality;
N7 defeats an unrestricted153N claim; a modular clock discards winding.
Least-used sidecar here: the tree decomposition of the *assembly skeleton*.
Board: complete boundary relations; separator width; full integer clock;
macro-state parity; source ownership. Anchor is the guarded Collatz
interface, niche the HN width obstruction, wildcard the paper acceleration.

## 1. Equality forcing requires a wide assembly skeleton

**Theorem.** Fix k>=2. Let H be a finite simple skeleton with distinct
marked vertices u,v. Replace each edge xy by a finite graph F_xy whose
only skeleton vertices are its two distinct ports x,y. All skeleton vertices are retained. There are no precoloured vertices,
colour lists, or restrictions that distinguish colour names. All fragment
interiors are mutually disjoint, and the assembly has exactly the union
of fragment edges (no additional cross-fragment constraints). Suppose
every F_xy has a proper k-colouring in which its ports differ. If H has
treewidth at most k-2, then the assembly has a proper k-colouring with
u and v different. Thus it cannot be a nonvacuous k-colour equality tester.

**Proof.** Take a tree decomposition of H of width t<=k-2. Choose a bag
containing u and a bag containing v. Add u to every bag on the connecting
path. This preserves the connected-bag property for u, changes no other
vertex's bags, and covers the edge uv. Thus H+uv has a decomposition of
width at most t+1<=k-1. Any graph of treewidth at most k-1 is k-colourable:
remove a leaf bag contained in its neighbour, or choose a vertex belonging
only to a leaf bag, of degree at most k-1, delete it and induct. A one-bag
decomposition is the immediate base case. Reinsert greedily. Therefore
choose a proper k-colouring of H+uv. On each fragment, permute a fixed
separating colouring so that the two port colours match this skeleton
colouring. A permutation of the k names maps any ordered pair of distinct
colours to any other. Disjoint interiors make these extensions compatible.
This proves the conclusion. The proof uses existence, not extension counts.

The threshold is sharp for abstract skeletons: H=K_(k+1) minus uv has
treewidth k-1, with bags V\{u}, V\{v}. Its other k-1 vertices form a clique
and each meets both ports. It is k-colourable, but all k-colourings identify
u,v. The K_k clique proves treewidth at least k-1 by the subtree Helly
property (or the elementary clique bag property). Edge fragments themselves
admit separating colourings. This sharpness witness is NOT claimed to have
a unit-distance realization for k=6.

**HN consequence.** At palette six, an assembly of separating two-port
fragments over any skeleton of treewidth at most four always admits a
separating colouring. This extends the inherited series/parallel no-go to
arbitrary bounded-width skeletons, including K5 skeletons. A viable
six-colour equality programme must have skeleton treewidth at least five,
a leaf already forcing equality, or overlaps/additional edges that violate
the two-port assembly interface. This is a necessary obstruction, not a
construction. Exact plane geometry remains a separate requirement.
A supplied unit-distance graph may retain only selected unit edges; extra
unit edges cannot silently be added to the assembly when applying the theorem.

## 1.1 Exact reduction: delete all soft two-port constraints

For each fragment, let e be the number of colouring extensions of one
fixed ordered equal-colour port assignment, and d the number for one
fixed ordered distinct-colour port assignment. Colour-name symmetry makes
these independent of the chosen equal colour or distinct ordered pair.

The width theorem has a stronger support-level form. Under the same
interface hypotheses and d>0 for every leaf, let H_D retain exactly those
skeleton edges whose fragment admits NO equal-port k-colouring (e=0).
Then the assembly's complete relation on the skeleton vertices is exactly
the proper-k-colouring relation of H_D.

Proof: a colour permutation shows that a fragment with e=0 supports every
ordered unequal pair and no equal pair; one with e>0 supports every pair.
Thus restriction of any assembly colouring satisfies H_D, and each
colouring of H_D lifts independently to all fragments. No multiplicity is
preserved or needed. In particular a separating assembly colouring exists
iff H_D+uv is k-colourable; a nonvacuous forced-equality assembly is exactly
a k-colourable H_D whose augmentation H_D+uv is not k-colourable.
This is an arbitrary-skeleton extension of the inherited two-port support
calculus, not a newly discovered formula for serial/parallel counts.

Consequently every nonvacuous equality assembly has chi(H_D)=k and
treewidth(H_D)>=k-1. The upper chromatic bound follows by restricting its
existing colouring. If H_D were (k-1)-colourable, recolour v alone with a
spare kth colour to separate the ports. The width bound then follows from
the elementary treewidth colouring lemma above. Thus at palette six the
*hard-inequality skeleton itself must already be six-chromatic*. In
particular replacing edges by highly equality-favouring but nonzero-d
fragments cannot create a hard constraint that is absent from H_D. Even
an enormous collection of all-positive two-port fragments has the full
boundary support, as long as their interiors really are disjoint.

The sharp K_(k+1)-uv example uses only hard-inequality edge fragments.
This exact reduction identifies the failure boundary: genuine progress
must construct a hard skeleton, use an already equality-forcing fragment,
or retain a coupled interface with more than two shared ports/cross edges.
It is not a restriction on arbitrary unit-distance graphs that lack this
particular decomposition.

## 2. The mod15/mod35 guards cannot reconstruct an arrival clock

For j>=0 put m_j=2+12j and n_j=(2^m_j-1)/3. Each n_j is a positive odd
integer, n_j=1 modulo105, and the accelerated odd Collatz map sends n_j
to1 with exact valuation v2(3n_j+1)=m_j. Indeed 2^12=1 modulo315,
so 2^m_j=4 modulo315 and division by3 gives n_j=1 modulo105.
The conventional unaccelerated map n odd ->3n+1, n even ->n/2 therefore
takes exactly m_j+1 steps to that displayed next odd value. For j>=1 this
is first arrival at1, since the intermediate values are powers of two
strictly above1. For j=0 it is a return to1, not a first arrival.

Consequently identical mod15/mod35 source registers, identical accelerated
successor, and even identical valuation phase modulo6 can coexist with
arbitrarily different physical clocks. No function of those observations
can recover the full clock. This is a proved infinite hostile family, not
just a finite simulation. It does not assert the repository's native P-map
is the ordinary Collatz map or identify the clocks of those different maps.
It explains why the existing authenticated native carrier must retain its
own actual clocks in addition to residue guards.

No credit is changed. The native ledger remains Delta_h=A_h+R'_h+J'_h,
A_h+M_h=153N. A certified upper bound B on Delta_h yields only
R'_h<=B-A_h-J'_h; B<=153N yields R'_h<=M_h-J'_h. Any moved negative
source-pair term must be subtracted from the old remainder exactly once.
Neither the clock hostile nor the paper acceleration pays any part of R'.

## 3. A second exact acceleration of the paper's A.9 recurrence

Source: Dos Santos et al., *At the Limits of Computability: Steps Towards
Autonomous AI Discovery on the Busy Beaver Frontier*, arXiv:2610.22603v2,
19 October 2026, Appendix A.9. The attached paper's rule system is the
input here; its Lean verification is a reported paper claim, not a new
Lean build in this packet. The parameters below are nonnegative integers.
The source configurations C(a,b,c) and rules are transcribed literally in
the checker. Work on the boundary b=0, writing (a,c) for C(a,0,c).

**Proposition.** The first positive return to b=0, or halt, is this exact
piecewise map, with L counting applications of the displayed A.9 macro
rules (NOT raw Turing steps, and NOT native Collatz P-depth):

- c=0: halt, L=1.
- c=1, a=2k+1: (2,k+3), L=1.
- c=1, a=0: (14,0), L=3.
- c=1, a=2k with k>=1: (17,k-1), L=3.
- c>=2, a=2k+e, e in{0,1}, c<=k+4:
  (6c-7+3e, k-c+4), L=c.
- c>=2, a=2k+e, e in{0,1}, c>=k+5:
  (6k+20+3e, c-k-5), L=k+4.

**Proof.** For c>=2 the first source rule gives
C(2+3e,k+3,c-2). Set r=k+3 and s=c-2. While both final coordinates are
positive, exactly min(r,s) rules add6 to the first coordinate and subtract1
from each of the other two. If s>=r this reaches b=0 directly after r+1
rules, giving the last branch. If s<r it reaches C(2+3e+6s,r-s,0);
one more rule gives C(5+3e+6s,0,r-s-1), for s+2=c rules.
All earlier intervening b-values were positive. For c=1 odd a, the first
rule already returns. For even a, the first rule gives C(5,2,k+1).
For k=0 drain once then use the c=0 rule, yielding C(14,0,0); for k>=1
drain twice, yielding C(17,0,k-1). The case c=0 is the halt rule.

The paper's initialized C(5,2,0) reaches C(8,0,1) in one macro rule;
the next boundary states begin (17,3),(14,9),(47,2). This checks that the
actual initialization enters the derived section. It proves no ultimate
halting or nonhalting claim. The parity e is indispensable; erasing it
changes both branches. A clock-weighted lift can attach and sum verified
raw-step weights to the traversed macro rules, but no such weights are
provided here. In particular L may not be substituted for physical time.

## 4. Connection contract and verification

The actual transfer is certificate design: source trace -> boundary return
with exact macro-length, and assembled colouring -> skeleton colouring
with exact endpoint relation. Both preserve a named predicate under an
explicit interface. The first loses intermediate states and raw time;
the second loses geometry and fragment multiplicities. The sidecars are
the full rule trace/true time and exact port overlap/geometry, respectively.
There is no map between the BB machine, ordinary Collatz, and HN that
transports convergence or chromatic conclusions.

`python3 verify_certificates.py` and `python3 -O verify_certificates.py`
run independent stepwise/closed-form recurrence checks on a,c in0..100,
exact Collatz clock controls j=0..20, and exhaustive k-colour counts for
the sharp abstract skeletons at k=2..6. The A.9 test checks first return,
not merely the eventual endpoint. For k=3,4, every absent/edge/two-edge-path replacement of the three
possible edges on three skeleton vertices is checked by independent full
colouring enumeration and restriction to H_D (54 assemblies). All checks
use explicit exceptions.
The theorem proofs cover all stated parameters; finite checks corroborate.
