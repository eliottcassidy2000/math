# Four vertices, six edges, and four port observations

Status: **PROVED** for the stated marked linear model; **FINITE-EXACT** field realization and exhaustive small controls. This is a coordinate bridge, not an identification of fields with tournaments or a tournament-class invariant.

## The inherited question

[Session86](../session86_relations_attachment_carriers_20261028/NORM_PORTS.md) identifies the nine norm-one elements of F64, under its binary alternating pairing, with the nonzero vectors in three mutually orthogonal binary two-dimensional subspaces. Write

\[
V=V_1\oplus V_2\oplus V_3,\qquad
U=(V_1\setminus0)\sqcup(V_2\setminus0)\sqcup(V_3\setminus0).
\]

Fixed block-local partner queries require six bits; adaptive partner queries require five; arbitrary linear queries require four. The four-query lower bound is simply that nine possible ports cannot be separated by three binary answers. The question here is whether the six coordinates and four answers have an actual K4 realization.

## The boundary-plus-parity map

Give K4 its six binary edge coordinates. Let B be its unsigned vertex-edge incidence matrix over F2. Define

\[
M(x)=((Bx)_0,(Bx)_1,(Bx)_2,\sum_e x_e).
\]

The omitted fourth boundary coordinate is determined by the other three. The map has rank four. Its kernel has dimension two and consists of zero and the three four-cycles of K4. Indeed Bx=0 describes the cycle space; its even-weight subspace consists precisely of these four vectors. Equivalently this is K4's bicycle space, the intersection of its cut and cycle spaces.

Partition the six edges into three pairs and identify the three V_i with those coordinate pairs. Then M separates U **if and only if all three pairs consist of adjacent edges**.

Proof: a collision between two U vectors is a nonzero kernel vector supported in at most two pairs. A four-cycle is the complement of an opposite pair. Thus it misses an entire coordinate pair exactly when that pair is opposite. Within a single pair the kernel has no nonzero vector. This proves both directions.

There are fifteen pair partitions in total and exactly eight successful ones. Inclusion-exclusion on the three possible opposite pairs gives 15−9+3−1=8. A successful partition has three distinct pair-centres. The omitted vertex is a root; the remaining three edges form a directed triangle when each is assigned to its pair-centre. Thus the eight choices are exactly four roots times two cyclic directions. Conversely, pair each root spoke rv with the outgoing triangle edge vv+.

Declaring the root a source identifies these eight scaffolds with the eight labelled tournaments 1→C3; declaring it a sink gives their converses. The source/sink choice is additional data, not determined by the wedge partition. Each scaffold has stabilizer C3 under vertex relabeling.

## An explicit F64 realization

Use F64=F2[a]/(a^6+a+1), with integer masks encoding coefficients in the polynomial basis. Take ordered binary bases

\[
(1,58),\quad(6,26),\quad(11,20)
\]

for the three orthogonal planes. Map these basis elements respectively to the edge pairs

\[
(03,01),\quad(13,12),\quad(23,02).
\]

The nine ports are {1,58,59,6,26,28,11,20,31}. With B1(x,y)=Tr(xy^8), the four rows of M become the four field probes

\[
48,29,25,56.
\]

The kernel in field coordinates is {0,9,35,42}. It misses every difference of two distinct ports, so the four answers separate all nine. The field-to-edge map is linear on all 64 vectors, not merely a lookup table on the nine ports.

The producer and an independent checker verified the field arithmetic, all 64 coordinate vectors, all 15 pair partitions, and the 36 avoiding two-dimensional kernels among the 651 binary planes. The independent checker also tested all 1,536 pairs of a K4 orientation and a vertex permutation. Reproduction scripts will accompany the complete session packet.

## Hostile controls and retained sidecars

1. **The natural opposite-edge pairing fails.** In the tested opposite-pair assignment, ports 28 and 31 have the same four-bit answer. Merely saying “six edges and four vertices” does not identify the successful map.
2. **The four-bit map does not classify tournaments.** Number edges 01,02,03,12,13,23 and use the transitive orientation as zero. Flip masks 0 and 45 have the same response but different isomorphism types. Across all 64 orientations, twelve output fibres each contain two transitive tournaments (H=1) and two strong tournaments (H=5). The other four fibres each contain four source-over-C3 or sink-under-C3 tournaments (H=3).
3. **Relabeling is affine in flip coordinates.** The vertex-boundary bits are changes in outdegree parity. Total flip parity changes by the permutation-sign cocycle when the reference transitive orientation is relabeled. A fixed orientation baseline must not be silently erased.
4. **This is not a multiplicative field quotient.** The nonzero proper kernel cannot be an ideal of F64. The faithful order-nine multiplication action on norm-one ports cannot be realized by K4 vertex relabelings, since S4 has no element of order nine.
5. **The selected alphabet is not closed under addition.** The map is lossless on U but not on arbitrary sums of its elements. Full additive closure is all six-dimensional V and needs all six bits for exact reconstruction.

## What the connection buys

The six-versus-four gap has a genuine marked tetrahedral realization: a rooted cyclic scaffold selects a successful encoding, and boundary plus edge parity supplies four linear tests. The information saving comes from restricting the source to nine one-block states. It does not recover arbitrary joint attachment information for free, identify tournament classes, or replace the exact transporters required by the 81/27 hostile pair.
