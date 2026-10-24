# Equal symmetry data, different path-fibre masses

Status: FINITE-EXACT explicit counterexample, proved by complete integer
Hamiltonian dynamic programs and full score-preserving permutation searches.
No minimal-order claim or novelty priority claim.

The two tournaments, with bit i marking an outneighbor i in each integer row,
are:

A = (102,108,408,657,771,220,412,787,555,103)
B = (102,108,664,401,771,220,412,787,551,107).

Both have the same two stable cells {0,1,2,3,4} and {5,6,7,8,9}, inducing
consecutive cyclic five-tournaments. Every vertex in the first cell has two
outneighbors in the second. Their score multisets are both five 4s and five
5s. Both crossing translation stabilizers are trivial. Consequently the
prime cyclic-cell decoder gives exactly the same graph, a single rigid edge,
and the same automorphism group order1, ten singleton vertex orbits, and
refinement error40.

Nevertheless H(A)=15465 and H(B)=15549. Their ordinary path-fibre masses
H/|Aut| differ by84. Thus they are nonisomorphic and the complete gain/pin
symmetry data, even augmented by these cell sizes, internal types, and all
crossing degrees, does not determine path-fibre mass, and therefore cannot reconstruct a full
metagraph object that retains those masses. This witness alone does not prove
different bare unweighted metagraph adjacency.
No explanation of the forbidden Hamiltonian counts7,21 follows from the
accidental divisibility of the difference84.

The mechanism is precise: trivial crossing stabilizer says which rotations
preserve a crossing, not which crossing matrix is present. The decoder
correctly preserves the automorphism action while discarding data needed
for Hamiltonian counts. Retaining all actual crossing matrices repairs the
information loss but sacrifices that compression.

Run python3 -B check.py and python3 -B -O check.py. Independent forward subset
DP and backward memoized deletion DP agree. Full ordinary automorphisms are
checked over all score-preserving permutations (5! squared per tournament),
not inferred solely from the decoder. The copied decoder source is pinned
in SOURCE_PINS.json; no network or sibling packet is required.
