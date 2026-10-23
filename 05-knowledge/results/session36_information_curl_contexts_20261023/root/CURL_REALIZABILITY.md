# Rooted triangle-curl realizability and an exact finite decoder

Status: PROVED by the elementary argument below; FINITE-EXACT tests are separate. Session36 research, not a complete metagraph-wiring realization theorem.

## Inheritance and meaning

THM-4739, `01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`, distinguishes the score-gradient from full signed triangle curl. Its low-order curl-profile summaries do not determine all wiring. Session35's augmentation form retained the full restriction to the sum-zero lattice but discarded its metric placement in some comparisons. The companion Session36 lattice note classifies the fibres of that full restriction. Here we solve the inverse problem for an arbitrary proposed rooted curl array, including arrays which no tournament realizes.

Fix labelled vertices 0,...,n-1 and root r=n-1. Let M be skew with M_ij=+1 precisely for i→j. Define

K_ij = M_ij − M_ir + M_jr  (i,j<r).

These are the signed triangle circulations around (i,j,r), or equivalently the matrix of the alternating form on e_i−e_r. All signs and the common root are retained. This is not the scalar count of cyclic triangles.

## The sign-constraint theorem

For n≥1, a supplied skew array K is realizable if and only if its upper entries lie in {−3,−1,1,3} and the following Boolean implication/pin system is consistent. Write a_i=M_ir∈{−1,+1}:

- K_ij=3 forces a_i=−1, a_j=+1.
- K_ij=−3 forces a_i=+1, a_j=−1.
- K_ij=1 requires a_i≤a_j.
- K_ij=−1 requires a_j≤a_i.

Every satisfying assignment gives exactly one realization by M_ir=a_i and M_ij=K_ij+a_i−a_j; every realization arises once.

Proof: substitute the possible differences a_i−a_j∈{−2,0,2} and require M_ij∈{−1,+1}. These four rows list all possibilities, independently for every pair. The resulting M is skew and has the requested K. No additional tournament axiom exists: exactly one oriented arc per pair suffices.

## A threshold decoder, not exponential assignment enumeration

Represent a_i≤a_j by a directed implication i→j. Reachability propagates a forced plus sign forward, and a forced minus sign backward. If any vertex is forced both signs, reject. Otherwise remove all forced vertices.

Every pair of remaining free vertices has K=±1: an entry ±3 would have pinned both endpoints. Thus their implication digraph is a tournament. Its strongly connected components have a unique linear condensation order B_1→...→B_b. Any satisfying assignment is constant on each component and consists of a minus prefix followed by a plus suffix. Conversely every one of these b+1 thresholds satisfies all implications, including those through pinned vertices: any prohibited connection to a fixed sign would already have propagated a pin.

Therefore a realizable array has exactly b+1≤n labelled tournament realizations. Reachability takes O(n³) by elementary transitive closure; reconstructing all at most n matrices also takes O(n³). Linear-time SCC algorithms improve the implication phase, but the stated bound needs no specialized implementation.

For n=1 the empty array has one realization; n=2 has two. The companion cyclic-SCC classification identifies b+1 with the number of strong components of any realized tournament. A strong tournament is uniquely reconstructed on its given labelled vertex set.

## Full triangle data and a minimal hostile

If instead all oriented triangle values C_ijk are supplied, first require alternation and the tetrahedron cocycle identity

C_jkl − C_ikl + C_ijl − C_ijk = 0.

Set K_ij=C_ijr and apply the decoder. The cocycle condition ensures every other triangle is K_ij+K_jk+K_ki. These tests are jointly necessary and sufficient.

Even alternation, the cocycle identity and the correct local alphabet {±1,±3} are NOT sufficient alone. On four vertices set K_01=K_02=K_12=−3. The remaining triangle also has circulation −3, so every triangle has an individually allowed cyclic value and the cocycle identity holds. But K_01 forces a_1=−1 while K_12 forces a_1=+1. No realization exists. This is minimal: for at most three vertices every proposed rooted array with entries in the alphabet is realizable.

Thus local circulation magnitudes need a common sign-potential sidecar. The obstruction is global compatibility, not a missing count formula.

## Scope and tests

`check_curl_decoder.py` uses explicit exceptions rather than optimization-removable assertions. It independently enumerates every labelled tournament through n=6, computes H by subset dynamic programming, groups by K, and compares each complete group with decoder output. It also tests every alphabet-valid proposed K through n=5, including impossible inputs. Results: distinct curl arrays for n=1..6 are 1,1,4,38,728,26704; every fibre has constant H and size equal to the number of SCCs.

These finite checks corroborate, not replace, the all-order proof. The decoder reconstructs a tournament from a quadratic-size labelled invariant. It does not provide a polynomial-time algorithm for H, compress arbitrary metagraph wiring, or solve the Hamiltonian-count spectrum problem.
