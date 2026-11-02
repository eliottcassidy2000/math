# Signed incidence carriers for legal path complements

**Status: PROVED** elementary compression and three-window formula below.
**FINITE-EXACT:** exhaustive minimal observers at orders 25, 26, 27; independently
compiled weighted spectra through order 30. These are labelled, guarded
tournament responses, not unlabelled tournament enumeration A000568.

## The object and its operations

Fix a finite labelled graph G with m edges, and a named bank of directed
spanning paths in G. A source is an orientation x of G. Path P is legal when
every one of its directed edges agrees with x. The observation is the entire
set of legal path names. On a legal source, K_P keeps P's edges and reverses
every edge outside P. This is the restriction of the actual tournament
path-complement operation. The binomial(n,2)-m outside-guard bits contribute
a common multiplicity and cannot affect these observations or operations.

Two sources are equivalent iff every legal finite command word gives the
same observations from both. This is the coarsest observation congruence.
There is one dead state when some orientation has no legal bank path:
no command is legal there. This is a deterministic observer with partial
commands, not merely a linear representation of a counting sequence.

Inheritance: [PATH_KERNEL](../joint_incidence_clocked_frames_20261102/PATH_KERNEL.md)
proves the joint signed-incidence expansion; [COMPLEMENT_CARRIER](../lower_frontier_complement_fibonacci_20261102/COMPLEMENT_CARRIER.md)
proves the matching-component formula and the minimum 33 at order 25.
The present extension admits intersecting three-path and larger components.

**Concurrent credit:** Session141 at9e1ae5aaf7 independently supplies the
same n26 polynomial and spectra through30, plus the general common-block
normal form and primitive binary contents. Its
[square-path proof](../session141_metric_aliquot_square_paths_20261102/square_paths/PROOF.md)
and [content proof](../session141_metric_aliquot_square_paths_20261102/square_paths_content/PROOF.md)
are integrated as matching evidence. The observer minima, signed-atom
compiler and command-word formulas are the distinct focus here.

## C1. Signed incidence atoms give an all-size closed finite compiler

For each edge e, record its column across all bank paths: entry 0 or 1 if
that path prescribes that orientation, and * if the path omits e. Discard
all-* columns. Normalize each remaining column by complementing all its
specified entries if its first specified entry is 1. Apply that same
coordinate flip to the source bit. Group edges with identical normalized
columns into atoms.

For an atom of size r retain only its source status:

- E: all normalized bits are 0;
- F: all normalized bits are 1;
- M: the bits include both 0 and 1, available iff r >= 2.

This is an exact observer abstraction for the specified bank and operations.
A path imposing value v on the atom is legal there exactly in state E (v=0)
or F (v=1); a * imposes no condition. A legal K_P fixes constrained atoms
and complements free atoms, exchanging E and F and fixing M. Thus legality
and every successor status factor through the abstraction. Induction on
command length proves all-future faithfulness.

The fibre weights are 1, 1, and 2^r-2, respectively. Different atoms use
independent coordinates, so source weights multiply; all-* coordinates
contribute a factor 2 per coordinate. Consequently the abstraction retains
exact source counts as well as future observations. It need not be minimal.
Ordinary deterministic partition refinement then produces the unique
minimal observation congruence up to state renaming.

With t atoms there are at most 3^t abstract sources, or more precisely
the product of 2+[r>=2] over atoms. Generate only active sources by taking
the union of the cubes allowed by the individual path guards. Both the
construction and its proof apply to every finite bank, including banks
whose compatible paths are not a clique. No asymptotic speedup is claimed
when almost every edge has its own column.

**Failure boundary.** Enlarging the bank can split an atom. For an old
two-bit free atom, sources 01 and 10 are both M, but a new guard asking
for first bit 0 separates them. New attachment operations therefore
require refined incidence data or the original source bits. This is not
a carrier for arbitrary future attachments without that sidecar.

## C2. Compatibility components are independent active controllers

Join two path names when their prescriptions are consistent. Simultaneously
legal paths are a clique in this compatibility graph. K_P preserves P, so
all paths reachable by subsequent legal commands remain in P's connected
component. Different active components are already distinguished by the
observed names. Hence

    minimum = [dead exists] + sum(component active minima).

Apply C1 componentwise, retaining the outside-component edge multiplicity.
An isolated path contributes one active state. Two compatible paths of
equal length, differing in d edges on each side, contribute
3+2[d>=2], as inherited from the matching theorem.

## C3. The three-window law

Consider three compatible paths P,Q,R. Use their common consistent
orientation on the union U as reference. For Z, the set of reversed
union edges, write B_P=U minus E(P), and similarly for Q,R. Then

    P legal iff Z is a subset of B_P;   K_P(Z)=B_P minus Z.

Assume the free sets have the following shape, with four pairwise disjoint
nonempty blocks A,B,C,D:

    B_P=A union B,   B_Q=B union C,   B_R=C union D.

Write m_X=1 if block X has size at least 2 and 0 otherwise. The exact
number of active observer classes is

    K_window = 8 + m_A + m_D + 3m_B + 3m_C + m_B*m_C.       (C3)

In particular all block sizes at least two give **17 states**, independently
of their magnitudes. This statement concerns existence of statuses; their
source weights continue to depend on the actual sizes.

### Proof and minimality

Use E,F,M for each block. The following list partitions all active sources.

| Initial legal set | Classes and count |
|---|---|
| PQR | All empty: 1 |
| PQ | B is F or M: 1+m_B |
| QR | C is F or M: 1+m_C |
| P only | A=F with B=E,F,M: 2+m_B; or A=M with all B states merged: m_A |
| R only | D=F with C=E,F,M: 2+m_C; or D=M with all C states merged: m_D |
| Q only | B and C independently F or M: (1+m_B)(1+m_C) |

There is no PR-only state. Every proposed merger is closed under legal
commands: in particular A=M keeps the controller forever at P only, and
D=M similarly at R only. Summing the table gives (C3).

Different legal sets separate immediately. In a PQ state, K_Q distinguishes
B=F from B=M by producing QR versus Q only; the QR case is symmetric.
In a Q-only state, K_Q distinguishes all four choices by producing,
respectively, PQR, PQ, QR, or Q only.

For P-only states, A=M stays P only under K_P. If A=F, K_P instead
produces PQ or PQR. B=F produces PQR; B=E and B=M produce PQ but are
separated by the next command K_Q. The R-only argument is symmetric.
Thus every two listed classes have a distinguishing legal word of length
at most two. This proves both sufficiency and minimality. QED.

For block sizes (2,2,2,2), the 40 active raw masks have class-weight
histogram {1:8, 2:6, 4:1, 8:2}. The independent finite refinement has
6,15,17,17 classes in successive rounds. All sixteen choices of block
sizes in {1,2}^4 and the asymmetric choice (2,3,4,2) verify (C3).

### Complete weighted command-word formula

There is a stronger structural description. Put
b=2^|B|-2, c=2^|C|-2 and

    J=(2^|A|-2)2^|B| + (2^|D|-2)2^|C| + bc.

The controller's undirected command graph is the disjoint union of:
one K1,3 with source weight 1 at each vertex; one four-vertex path P4
with weight 1; two K1,2 stars with weights b and c; and up to three
single loop states whose total source weight is J. A zero-weight
component is absent. Each graph edge is a named involutive command.
This decomposition follows directly from the six rows of the table:
PQR and its full-window images form K1,3; the full B/C overlap states
and their end states form P4; each mixed middle block forms a K1,2;
the mixed end blocks and the double-mixed middle state give the loops.

Define A_k(t)=k^ceil(t/2)+k*k^floor(t/2), and Fibonacci numbers by
F_0=0,F_1=1. Let W_t count pairs (active source, legal named command
word of length t), with no terminal restriction. Then for every t>=0,

    W_t = J + (b+c) A_2(t) + A_3(t) + 2 F_(t+3).          (C3W)

A star K1,k has exactly A_k(t) total walks; P4 has 2F_(t+3), obtained
by its reflection-even two-state recurrence. Thus Fibonacci enters
through an actual four-vertex path controller, not a numerical match.
The formula applies to arbitrary positive block sizes, including size 1.

Its minimal scalar recurrence order is exactly

    4 + 2[b+c>0] + [J>0],

although its minimum deterministic observer is (C3). Indeed the active
spectral factors are x^2-3, x^2-x-1, optionally x^2-2 and x-1.
Each star has strictly positive coefficients on both square-root
eigenvalues; P4's all-ones vector has nonzero components on the two
reflection-even eigenvalues. The factors are distinct, so no cancellation
reduces the stated order. When all four blocks have size two, the
17-state observer therefore has a **seven-dimensional count recurrence**.

## C4. Actual square-sum guards

The vertices are the integers 1 through n, and an unordered edge ij is in
the guard iff i+j is a square. The bank contains every directed Hamilton
path of this guard, including the two directions as separate names.
Every pair of distinct vertices still has a tournament orientation; only
guard edges enter the observable.

At 26 the compatibility graph is 6 K1 + 6 K2 + 2 K3. Each K3 has precisely
the C3 window arrangement with all four blocks of size two. Of the six
pair components, four have d=2 and two d=5. Thus

    K_26 = 1 + 6 + 6*5 + 2*17 = 71.

The 11,232 active guard orientations refine through 36,66,70,70 classes.
The dead state adds one. Active class weights are
{16:6,32:16,64:12,128:14,256:12,480:4,512:6}.
The named paths, actual edge masks, source partition hash and all weights
are in verification.json. A separate endpoint/subset dynamic program
confirms the 24 directed paths without using the inherited leaf filter.

At 27 the minimum is **637**, by literal exhaustive refinement and,
independently, C1's weighted atom compiler. Its component minima are

    8*1 + 6*5 + 2*16 + 2*17 + 2*39 + 2*52 + 2*175 = 636;

one dead state completes the count. The 62,952 active raw orientations
refine through 150,456,596,632,636,636 classes.

**Hostile to forgetting incidence:** the two three-path components at 27
are also K3, but each has minimum **16**, not 17. One has path IDs
(7,15,18), union size 32, and free-membership atom sizes
{P:4,Q:1,R:1,PQ:1,PR:1,QR:4}; its twenty common edges are fixed.
Thus the same compatibility triangle supports different controllers.
An early probe that tried to reuse C3 solely from the K3 graph was rejected
by its window-shape guard; the signed-incidence compiler handles both.

## Exact path-count polynomials

Let h(x) be the number of legal paths. The following are coefficients of
sum_x z^h(x) on guard orientations. Multiply by 2^(outside free bits) for
all labelled n-vertex tournaments.

| n | Directed paths | Outside free bits | Coefficients from constant term upward |
|---|---:|---:|---|
|25|20|268|4294962432, 4608, 256|
|26|24|291|17179857952, 10240, 928, 64|
|27|70|315|68719413784, 54976, 7360, 488, 120, 8|
|28|104|340|274877723488, 158096, 21840, 2928, 528, 64|
|29|38|366|1099511486464, 129024, 10240, 2048|
|30|40|393|4398046212096, 274432, 20480, 4096|

Two algorithms agree exactly: compatible-family union moments with binomial
inversion, and a separate Shannon recursion retaining residual-pattern
multiplicity. The first moment is independently fixed by
number_of_paths * 2^(m-n+1), and the total by 2^m.

The concurrent primitive normalizations are exactly equivalent:
Z25=2^276(16777197+18z+z^2), and
Z26=2^296(536870561+320z+29z^2+2z^3).
Their binary-content exponents differ by20. This occurrence of20 is
proved from coefficient gcds and joint multiplicities; it is unrelated
to the incorrectly proposed honeycomb deficit20.

The 81/27 pair remains a separate hostile control: it has the same 67
labelled cyclic triples but automorphism groups of orders 81 and 27.
C1 preserves the explicitly supplied operation language; it makes no claim
that a cyclic-triple deck or this path observer reconstructs a tournament.

## All-length recurrences for the actual guards

Write R_n(t) for the number of (active guard orientation, legal command
word of length t) pairs. Component decomposition and C3W give, for all t>=0,

    R_25(t) = 4096 + 256 A_2(t),
    R_26(t) = 8320 + 800 A_2(t) + 64 A_3(t) + 128 F_(t+3).

Multiply by 2^268 or 2^291 to include the outside tournament bits.
These formulas start from active sources. If every source is allowed
with the empty word, add the dead-source mass only at t=0.
They are different observables from the polynomial in the preceding table.

In particular R_25 obeys the minimal recurrence

    R(t+3)=R(t+2)+2R(t+1)-2R(t).

R_26 has minimal annihilator

    (x-1)(x^2-2)(x^2-3)(x^2-x-1)
      = x^7-2x^6-5x^5+11x^4+6x^3-17x^2+6.

Its initial seven values are
11232,12288,14528,16896,21888,27264,38784.
An independent exact rational Krylov elimination annihilates the entire
controller vector, so these are certified recurrences rather than fits
to a sequence prefix. Literal source walks verify the initial counts.

More generally, every fixed-bank controller supplies a finite matrix M,
where M_ij counts the legal names taking class i to class j. Each name
acts as a partial involution, so M is symmetric and neighbouring classes
have equal source weights mu_i. Then R(t)=mu M^t 1. The Hankel matrix
R(i+j) is a Gram matrix of M^i 1 in the positive weighted inner product.
This supplies positivity constraints and distinguishes a small linear
counting recurrence from the larger exact observer needed to answer
named future commands. It does not erase the mixed-bit split required
when a new attachment enlarges the bank.
