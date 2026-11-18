# A two-coordinate contraction for the hypercube collision sum

Let U_d be the fair-landmark union bound obtained by taking, for
each unordered edge-pair orbit, the maximum-weight spanning forest
of its distance-level graph and multiplying the central binomial
atoms beta(N)=2^(-N) binom(N,floor(N/2)) across that forest.
The orbit formulas and cell convolution used below are from
[Allikvere](https://arxiv.org/html/2608.09983v1), Sections 6–7.

**Claim.** U_(d+2) ≤ (3/4) U_d for every d≥11. Together with
the exact base values U_11<4/25 and U_12<11/1000 this gives a
single estimate from d=11:
  U_d < (4/25)(3/4)^floor((d−11)/2).
In particular a resolving set exists. The two rational base checks
are finite calculations, so this is a uniform analytic recurrence
with finite base certificates, not a completely computation-free proof.

## All pair types and the lift

There are d−1 parallel orbits P_h (1≤h≤d−1), and d−1 cross-direction
orbits C_h (0≤h≤d−2). Write g for the number of coordinates common
to the two distance readings: d−1−h for parallel pairs and d−2−h
for crossed pairs. Each type in dimension d+2 is assigned a parent
in dimension d: keep h if g≥h; otherwise decrease h by 2.
The latter is legal, including both new high-offset types.
The cell polynomial is multiplied by (1+z)^2 on the chosen group.
Consequently the old level graph, shifted by +1 on both indices,
is a subgraph of the new graph and every corresponding cell has
at least twice its old weight.

The parallel level graph has a−b≡h mod 2, so it has two parity
components for even h and one for odd h, except P_(d−1), where
it is an antipodal matching. Connectivity follows because the
level graph contains all edges of a fixed positive step, and, if
h is odd and h<d−1, also an edge joining the two parity classes.
The crossed graph contains all consecutive-level edges: using
the independent two distinguished bits, fix p so h−2p is zero
or one and use their different values when necessary. Varying
the common-bit count fills the entire chain. Thus the respective
forest ranks are d−2, d−1, and floor(d/2) in the parallel even,
parallel odd/cross, and antipodal cases. The replay checks these
closed classifications independently through d=24.

For every old non-antipodal forest the lifted tree(s) omit both
extreme levels 0 and d+1. Both new extremes have a level edge:
choose all relevant common and differing bits uniformly at either
end, with the two distinguished bits chosen oppositely for crossed
pairs. Each new edge weight is at least 2 and beta(2)=1/2. Attaching
both extremes is acyclic. For antipodal P_(d+1), the lift adds one
disjoint extreme matching edge of weight 4, giving beta(4)=3/8.

For N≥2 the exact binomial atom obeys beta(2N)≤(35/48)beta(N):
the quotient is maximal at N=3 (N=2 gives 3/4? For N=2 it equals
3/4, but all parallel weights are multiples of 4 and cross
weights multiples of 2; see correction below). In fact cross
weights may equal 2, so use beta(2N)≤(3/4)beta(N) universally,
while parallel weights N≥4 and divisible by 4 yield the stronger
35/48 bound. Parallel non-antipodal forests have at least d−2≥9
edges, crossed forests at least d−1≥10; hence each child's
collision factor is at most (35/48)^9/4 or (3/4)^10/4
times its parent's, respectively. Antipodal matching has at least
five edges, giving (35/48)^5·3/8.

The total orbit multiplicity mapped to a parent is at most
42 times the parent's for non-antipodal parallel types and
52 for crossed types: the binomial coefficients satisfy
C(n+2,h)≤4 C(n,h) for h≤(n+2)/2, and
C(n+2,h+2)≤4 C(n,h) for the complementary side.
There are at most two children, and the additional factors from
coordinate choices and 2^d are respectively at most
4(d+2)/d and 4(d+2)(d+1)/(d(d−1)).
The antipodal orbit ratio is 4(d+2)/d≤52/11<5.
The resulting three factors are
  42(35/48)^9/4 < 3/4,
  52(3/4)^10/4 < 3/4,
  5(35/48)^5(3/8) < 3/4.
Summing over parents proves the contraction including both
new boundary offsets. The replay evaluates the exact orbit sums
and checks the lift for finite sample dimensions; the proof is
the general binomial convolution and graph classification above.

## The number five and the order-six metagraph

At Q6 there are exactly five parallel and five crossed pair orbits,
indexed by distances 1..5 and 0..4. They are ten different
obligations, not six parallel types. The earlier six-family bound
is for d≥11, where offsets 1..5 and the antipodal offset are distinct.

The literal merged tournament metagraph at order six has 56
oriented classes and 34 converse-merged nodes. Its pure-blue
counts are 1 and 9; its transitive node has one blue incidence
to count 17. The replay enumerates every fixed-frame tiling and
recounts the Hamiltonian paths. Prior [THM-4730](../../../01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md)
states that the only strong pure-blue atoms are 1, C3, T5.
Using its corrected all-order classification, the bivariate
enumerator for pure-blue classes by order z and path count H is
  (1+z+z^3 H^3+z^5 H^15) / (1−z²−z^6 H^9).
A reflected pair of C3 atoms contributes order six and count nine,
explaining the order-six point. Marking a C3 pair by x and a
singleton pair by y gives multiplicities binom(i+j,j); their
parity is one exactly when the binary additions of i and j have
no carries. This is a rigorously periodic, dyadically self-similar
substructure, not a classification of the entire metagraph.

A testable broader hypothesis is that after retaining the strong
atom and endpoint-interface data, larger pure-blue fibres admit
the same finite-state binary carry description. The enumerator
above proves it for this stratum, but says nothing by itself about
mixed/black fibres or Chan–Pak coverage. The earlier fivefold
archive explicitly distinguishes the fivefold lattice character,
five Platonic solids, and the list of known Fermat primes; the
common cardinality alone does not identify them with Q6 or T5.

Replays: python3 cube_transfer.py and python3 blue_period.py.
