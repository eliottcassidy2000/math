# Tiling columns as two-height walks

October 23, 2026. **PROVED:** the exact coordinate change, realizability test,
and reflection/complement orbit formulas below for every n >= 3.
**FINITE-EXACT:** exhaustive replay of every fixed-path tiling for n=3..7.
This is a structural encoding, not a faster general tournament census.

The closest inherited mechanism is the fixed-path tiling cube and its two
commuting involutions. The current [signed-curl reconstruction](../metagraph_port_rigidity_20261023/PROOF.md)
supplies the broader labelled theorem. The canonical hostile here is two
locally permitted triangle curls that demand an impossible arc value.
The corrected near miss is treating allowable local values as a globally
realizable tournament. The least-used coordinate is the cumulative height
inside one tiling column. Live concepts: binary arcs; triangle curls;
cumulative capacity; reflection; complement; actual port ownership.

## Exact coordinates and the capacity constraint

Label the directed Hamiltonian path 0 -> 1 -> ... -> n-1. Write a_ij=1
when i beats j, a_ji=-a_ij. Thus a_i,i+1=1. For j >= i+2 put

    d_ij = c(i,i+1,j) = 1 + a_i+1,j - a_ij.

There are exactly m=(n-1)(n-2)/2 such coordinates, one per free arc. Set
h_ij=(1-a_ij)/2, so h_ij is 0 or 1 and h_j-1,j=0. Then

    d_ij = 1 + 2(h_ij-h_i+1,j),
    sum_{k=i}^{j-2} d_kj = j-i-1 + 2h_ij.                (1)

Consequently an integer array d is the array of a fixed-path tournament
**if and only if**, for every j and i <= j-2,

    sum_{k=i}^{j-2} d_kj is either j-i-1 or j-i+1.        (2)

Necessity is telescoping. For sufficiency define

    a_ij = j-i - sum_{k=i}^{j-2} d_kj.

Condition (2) makes every free arc exactly +1 or -1; set adjacent arcs to
+1 and antisymmetrize. Its consecutive triangle curls are the supplied d.
This proves both existence and uniqueness. The inverse needs O(n^2)
integer additions by carrying suffix sums down the columns.

Thus each column is a walk between heights 0 and 1, starting at height 0
next to the forced path. A step of +1 consumes its sole free height; another
+1 must wait for a -1. Zero steps leave the height unchanged. Every sequence
of the binary heights is allowed, so this describes exactly 2^m tilings.
It does not reduce their number or their information content.

At n=4, d_13=3 and d_03=3 are each in the local alphabet {-1,1,3}, but
together yield a_03=3-(3+3)=-3. They cannot coexist. Setting d_02=1
completes a smallest cross-coordinate obstruction with every individual
triangle legal for its position.
This is a precise capacity obstruction that bucket totals do not detect.

## The two reflections and the blue/black totals

Let tau retain path arcs and reverse every free arc. In height coordinates,
it toggles every free bit. In curl coordinates it is the affine reflection

    d'_ij = 4-d_ij  if j=i+2,
            2-d_ij  if j>i+2.                           (3)

The extra boundary term occurs because the adjacent height is fixed at 0,
not toggled. It is the same boundary issue that the rooted curl formula
records by its path-flow term.

Let sigma reverse path positions and take tournament converse. On forward
arc/height coordinates it permutes (i,j) to (n-1-j,n-1-i). It commutes with
tau. A tiling is blue exactly when sigma fixes it, equivalently when path
reversal is an anti-automorphism. The f=floor((n-1)/2) free positions on the
antidiagonal are sigma-fixed; the other m-f positions occur in pairs.
Hence, with k=(m+f)/2,

    blue tilings = 2^k,
    literal blue tau-lines = 2^(k-1),
    literal black tau-lines = (2^m-2^k)/2.              (4)

Because f>=1, tau*sigma cannot fix any tiling: a fixed antidiagonal bit
would have to equal its own complement. Every nonblue orbit of the group
{1,tau,sigma,tau*sigma} therefore has four elements; every blue orbit has
two. After pairing the two reflection-related black lines, the counts are

    fundamental black lines = (2^m-2^k)/4,
    fundamental blue lines  = 2^(k-1).                 (5)

These are line **multiplicities**, including parallel lines and loops;
they are not counts of distinct endpoint supports. Nor is tau full converse.
Converse merging of tournament nodes and the pairing of black lines must
retain these distinctions. This recovers inherited orbit-count mechanisms
in explicit height coordinates, with no claim of new global count formulas.

| n | tilings | blue tilings | fundamental blue lines | fundamental black lines |
|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 1 | 0 |
| 4 | 8 | 4 | 2 | 1 |
| 5 | 64 | 16 | 8 | 12 |
| 6 | 1024 | 64 | 32 | 240 |
| 7 | 32768 | 512 | 256 | 8064 |

## What this says about counterfeit metagraphs

A proposed wiring can satisfy all totals (4)-(5) and still fail to be the
tournament metagraph. The missing information is which **actual** binary
tiling, or equivalently which legal height array, occupies each node port.
Once these arrays are attached, (3) fixes every partner; sigma fixes color
and the paired-black quotient. A different attachment fails this exact
local transformation even if all aggregate degrees and parities agree.

The source is a path-labelled tournament; the target is its column walks.
The map preserves all adjacency, path, color and complement data. It loses
no information. Dropping the path labels or cumulative heights is a further
quotient and needs a separate restoration coordinate. The cheapest decisive
test is (2), followed by the complement formula (3).

The incoming [all-flips memory theorem](../session33_continuation_states_20261023/metagraph-memory/ADDENDUM_FULL_FLIPS.md)
shows why this should not be advertised as a smaller exact state space:
with every labelled free-arc flip allowed, all 2^m tilings are distinguishable.
Here the gain is an explicit integrability constraint and a decoder, rather
than fewer states. The [incidence-moment packet](../metagraph_bucket_obstruction_20261023/README.md)
instead tests how far a genuinely lossy summary can force small-order wiring.

Run `python -B 05-knowledge/results/metagraph_two_height_20261023/check.py`
and the same command with `-O`. The consumer exhausts every tiling at n=3..7,
checks both inverses, direct complement and reflection against the formulas,
and independently constructs the involution orbits. It includes malformed
and impossible arrays. Larger n are covered by the elementary proofs, not
by extrapolation of these finite checks.
