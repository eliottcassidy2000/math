# Combining the Moore constraints: frozen spectrum, global matching screens, and anchored action

2026-11-26. **PROVED** rigidity, spectral and collision decompositions;
**FINITE-EXACT** matching and marked-board tests. The degree-57 graph
and the unrestricted marked completion remain open. The current literature also distinguishes [obstructions to a regular derangement construction](https://www.mdpi.com/2227-7390/14/23/4318) from nonexistence. This continues the
[Wilson action](WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md),
[signed residue lift](HEIGHT-RANK-AND-SIGNED-LIFTS.md), and
[three-star fixture](NESTED-CORES-AND-RIGID-FOURTH.md).

## 1. Frozen-star rigidity has a quantitative Wilson version

Let P be a complete zero-action chart, and Q another complete
permutation chart, **not** assumed to be zero action. Freeze a set F
of m branch stars: each link incident to them is the same in P and Q.
At a sheet x of an unfrozen branch i, let s_x be the number of
branch neighbors whose sheets changed. Then the action on the
blocks between frozen and unfrozen branches satisfies

    E_cross(Q) >= 2 Σ_x s_x max(m+1-s_x, 0).                  (1)

For each changed exit j from x, call its old target u in branch j.
For each frozen branch k, there was a route to the uniquely
determined neighbor w_k of u in k. It cannot still be reached via j
because the k-links are frozen and bijective. An unchanged exit or
the direct route to k cannot supply it either, since there were no
collisions in P. A changed exit l other than j can recover w_k only
if its new target y_l and u share w_k as an old common neighbor.
For a fixed l this can happen for at most one k, by the Moore
property of P. Hence at least max(m+1−s_x,0) receivers are missed
for each changed j. The missed routes for different j are distinct
targets within any fixed block row. The number of holes in an
integer row of total 56 equals its total excess; its squared
defect is at least twice the holes. This proves (1).

At action zero it recovers the old minimum of m+1 changed exits
through every affected leaf. It also quantifies near misses:
leaves with only one changed exit each guarantee at least 2m
units on these blocks. From a single permutation change moving r
sheets, take all k−2 other branches frozen. The bound gives exactly
4r(k−2) on the frozen blocks; the direct changed block adds 2r,
matching the single-link formula. The reference P in (1) must be a
complete Moore chart. It is not assumed to exist for k=57.

## 2. Exactly what the three frozen stars already determine spectrally

Let X be the full leaf adjacency restricted to the real subspace
of vectors with sum zero **within each branch**. Permutation links
preserve this subspace. The Moore route condition there is

    X² + X = 56 I.

Zero action would give eigenvalues 7 and -8 on the 3135-dimensional
space W. Its trace is zero, since the diagonal branch blocks vanish;
the multiplicities must be 1672 and 1463.

Separate the three frozen branches (56,0,28) from the other 54.
Write the blocks of X as [C B; B^T D], so C and B are fully known.
The necessary equations include

    BB^T = 56I - C² - C,
    DB^T = -B^T(C+I),
    D² + D = 56I - B^TB.

The first identity is checked directly on the actual three stars:
before restriction to W, its right side has an additional J_56
in each off-diagonal frozen branch block. The link maps from 56
to 0 and 28 are identities, and the 0-to-28 map is the half-turn.
Splitting into the half-turn's even and odd spaces, with the
constant removed in each branch, gives the spectrum of C

    2^27, (-1)^54, 1^56, (-2)^28.

Consequently BB^T is invertible with eigenvalues 50, 56 and 54
of multiplicities 27, 54 and 84. The subspace image(B^T) among
unfrozen fibers has dimension exactly 165. It carries D as a
conjugate of -C-I. Via the isometry B^T(BB^T)^(-1/2),
each eigendirection of C with eigenvalue λ yields a paired block

    [ λ                         sqrt(56-λ²-λ) ]
    [ sqrt(56-λ²-λ)            -λ-1          ] .

The characteristic polynomial of each is t²+t-56, so the
three fixed stars already account for **165 copies of 7 and
165 copies of -8**. Only a 2805-dimensional sector remains,
with multiplicities 1507 and 1298.

In particular define

    D0 = -B^T(C+I)(BB^T)^(-1)B,
    Π  = I - B^T(BB^T)^(-1)B.

Every possible completion must satisfy

    D = D0 - 8Π + 15Q,    Q an orthogonal projection in Π of rank 1507.

This is an exact low-information/high-information split: the 165
frozen modes require no choice at all, while the unknown projection
must still have entries coming from the same discrete permutation
links. An arbitrary real projection does not supply a completion.

## 3. A whole-chart matching screen and a useful negative control

We rebuilt the physical graph on 3250 vertices. For every
unfinished pair of exterior branches, take the legal individual
edges (exclude already occupied branch slots and old routes of
length <=3). A genuine completion requires a perfect matching on
the unmatched rows of each such pair. The replay builds a matching
and removes any edge absent from every perfect matching: relative to
a known perfect matching, a nonmatching edge is usable iff it lies
on an alternating cycle, detected by mutual directed reachability
after contracting the known matching.

Then for every already completed vertex v and every other vertex u
not reached within distance two, require a legal new edge from u
into N(v). This necessary test keeps 240 complete vertices in view
at once, rather than selecting contacts for just the marked source.

| Fixture | Legal optional edges | Failed pair matchings | Empty frozen-radius demands |
|---|---:|---:|---:|
| Fourteen planted neighborhoods, with marker | 3,607,365 | 0 | 0 of 673,728 |
| With the old 51-contact patch | 3,595,917 | 2 | 8 of 673,328 |

The base fixture has no forced edge and no pairwise matching-support
deletions. The pair-domain checks alone have not settled it. In
the old patch the failed pair boards are exactly (1,13) and (1,33);
all eight empty frozen demands involve the source (1,1).
That patch is a good lower-bound witness on how many compatible
contacts can be placed, but **cannot itself extend** while retaining
all those contacts. This is a direct global-screen replay rather
than a failed heuristic run.

## 4. Signed phase and coarse capacity split the same action

For a marked full one-contact-per-branch selection, consider a
double-capacity residue with sheet occupancies a,b. Its physical
energy is exactly

    (a-1)²+(b-1)² = [(a+b-2)²+(a-b)²]/2.

The first term measures the residue-count defect, the second the
phase imbalance. A singleton contributes its whole (occupancy-1)²
to the coarse term. When every residue capacity is met, each double
residue has two incident contacts; half of its physical defect is
one exactly when its signed XOR equation fails. In that case
minimizing the number of squares is exactly minimizing the number
of violated signed edges with the legal phases held fixed. The
all-even long-seam criterion is the zero case of this identity.
For the old least-collision choice (13,16),(33,19) in the two
unfilled branches, the total physical action is 6, with 3 coarse
and 3 phase units. The failure is not solely a phase obstruction.

### Anchored exact search instead of changing arbitrary subsets

The old optimization scanned subsets of reference edits. The
new search starts at the two holes and opens a reference row only
when a selected route displaces one of its occupied slots. It
also allows a **soft** conflict: keep the reference contact
permanently and charge its collision. A collision with another
already inserted choice is also charged; hit reference contacts
are distinct within the three slots by board linearity.

Why this certifies minima: among solutions under any fixed
collision cap, take one with minimum reference edits. Draw directed
arrows from each new chosen contact to each changed reference
contact whose old slot it hits. If any changed reference row
is not reachable from the two holes, revert all such rows to
their old contacts. These are pairwise compatible with retained
references and do not hit any reachable new contact (or would
have been reachable). Reverting cannot add any collision, contrary
to the minimum edit assumption. Thus it suffices to branch only
on reference owners hit by a selected contact. At the first hit
a row is either opened for a change or retained permanently.

Exact results on the marked board:

* **At most nine reference edits:** the minimum is still three
  physical squares. At the ninth budget the cap-two tree exhausts
  after 56,085,562 nodes.
* **At most eleven reference edits:** no full branch selection has
  zero or one square. At the eleventh budget the cap-one tree
  exhausts after 56,071,651 nodes.

Here a selection with nonzero squares is a deliberately penalized
candidate, not a valid new patch. No negative result is claimed
outside either edit budget.

## Replay

From this directory, using Python and a C++17 compiler:

```sh
python nested_core_check.py --board > board.txt
c++ -O3 -std=c++17 frozen_hall.cpp -o /tmp/moore-hall
python export_frozen.py marked | /tmp/moore-hall
python export_frozen.py patch | /tmp/moore-hall
python frozen_spectrum_check.py board.txt
c++ -O3 -std=c++17 anchored_wilson.cpp -o /tmp/anchored-wilson
/tmp/anchored-wilson 2 9 9 120000000 < board.txt
/tmp/anchored-wilson 1 11 11 100000000 < board.txt
```
