# Source-preserving inverse growth obstructs a bounded-length planar lattice clock

**Status: PROVED elementary consequence of an inherited inverse tree;
FINITE-EXACT clock controls. Collatz remains OPEN.**

Fix the native positive-integer map C(n)=n/2 for even n and3n+1 for odd n.
Every node has exactly one outgoing edge. Its incoming nodes are2n and,
when n=4 mod6, the additional odd node(n-1)/3. Thus the requested indegree
bound is correct. It does not establish a unique orbit through all integers.

## 1. Recover the older binary tree instead of rediscovering it

[Source refinement, Proposition SF5](../collatz_source_refinement_floor_20261005.md)
already constructs a proper inverse binary tree rooted at5. For an odd unit
n>1 modulo3, choose

    k in {2,4,6} if n=1 mod3;
    k in {3,5,7} if n=2 mod3,
    child=(2^k n-1)/3,

and discard the one choice divisible by3. Exactly two remain: 2^k cycles
through the three appropriate units modulo9. Both children exceed n, and
their accelerated odd-parent image is n. Unique forward parents and strict
growth give an embedded binary tree. This is the inherited construction.

The refinement needed for the lattice question is its literal clock:
each child reaches n in k+1<=8 NATIVE steps of C. The exact exponent table is

| n mod9 | Two allowed inverse exponents |
|---:|---|
|1|2,4|
|2|3,7|
|4|2,6|
|5|3,5|
|7|4,6|
|8|5,7|

The verifier checks all8191 distinct sources through depth12, the two-child
guard, strict growth, and every native return word. The all-depth argument
is the preceding arithmetic proof, not extrapolation from those8191 nodes.

## 2. A metric obstruction, not a nonplanarity assertion

**Proposition.** There is no injective map of all positive Collatz nodes
into a fixed finite-dimensional lattice of polynomial ball growth for
which the endpoints of every native Collatz edge have distance at most a
single uniform constant L.

**Proof.** The depth-d inverse binary tree has2^d distinct nodes, each at
native undirected distance at most8d from5. Under the supposed map they
all lie in a lattice ball of radius8Ld, which has O(d^D) points in a
fixed D-dimensional polynomial-growth lattice. This contradicts2^d>O(d^D)
for sufficiently large d. QED.

For a finite planar realization through inverse depth d, if the largest
native edge length used is L_d and balls obey |B(r)|<=K(r+1)^2, the proof
already forces

    L_d >= [sqrt(2^d/K)-1]/(8d).

Thus a planar lattice picture can retain source names only by allowing
unbounded edge lengths, collisions/multiple sheets, or a growing hidden
address. The proposition applies both to vertices of the honeycomb and
to centers of hexagonal tiles (the triangular lattice). The additional
1->4->2->1 cycle rules out a unit-edge honeycomb embedding immediately,
but not a triangular-lattice one; the growth argument handles both.

If uniform L is kept but collisions are resolved by a binary sidecar of
at most b_d bits through depth d, source recovery requires
2^b_d K(8Ld+1)^D>=2^d. Hence b_d>=d-D log2(d)-O(1). This is a precise
growing-memory lower bound under the stated geometric clock assumptions;
it does not apply to a representation allowed to forget source identity.

This says nothing against abstract planar drawings with nonuniform curves.
Nor does it exclude a visual ordering that is not the actual Collatz step.
Forward merging alone already prevents an injective forward orbit from
being a space-filling enumeration of every positive integer. Collatz's
termination/root-cycle conjecture is neither assumed nor established.

## 3. Which hexagons and 3,5,16+/-1 clocks genuinely survive

The old [coherence note, THM-4521(2)](../collatz_coherence_20260930_trichotomy_hexagon_motifs.md)
has a literal six-state cycle: the units modulo9 under multiplication by2,
1,2,4,8,7,5. Its half-turn is negation. That finite residue/Fourier clock
is not an injective embedding of the integer graph.

The [Eisenstein-torus note, Theorem B2](../collatz_pentagon_fixed_points_argument_styles_20260930.md)
has a different map, multiplication by1-omega, under its local-hexagon and
no-extra-hexagon hypotheses. Retaining those incidence hypotheses is
essential. It supplies no identification with native Collatz edges.

For the nearby tuple, 3*5=16-1 and3*5+1=16 give the actual native trajectory
5->16->8->4->2->1. The inherited [Fermat/Pepin lane](../s58_fermat_pepin_collatz_20261026/README.md)
contains the complementary modulo17 clock. For the unhalved affine map
f(x)=3x+1, the change y=2x+1 gives y(f(x))=3y(x) modulo every odd modulus.
Modulo17,3 has order16, giving one fixed point and one16-cycle. Modulo15,
3 is not a unit; the image has only five residues. The verifier checks both.

These are operations with different clocks. Reinstating the native oddness
guard and exact halving valuations is necessary to return from an affine
residue permutation to integer Collatz. In particular, the17-cycle does
not prove integer termination, and the nonunit modulo15 chart cannot be
treated as reversible. The division guard, not the neighboring integers
alone, explains the difference.

## 4. The attachment-carrier connection

The same obstruction test is useful in metagraph design: if a proposed
low-dimensional geometric carrier admits only polynomially many states
within a bounded operation radius, a source-preserving binary family with
uniformly bounded legal return words defeats it. One may instead count
anonymous responses on a smaller quotient, just as a fixed marker in the
Conway example needs fewer states than the complete frame. The discarded
source address must stay discarded by every later allowed operation.

Applying this general test to actual strong-atom complement families
requires proving their own embedded binary source family and bounded legal
clock. That tournament application is OPEN; the Collatz tree by itself
does not provide those hypotheses.
