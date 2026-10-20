# SCC basin normalization, finite deficit centres, and golden words

20 October 2026, continuation d. This is the maintained synthesis of the
packet. The graph is the converse-merged Hamiltonian-path-complement
**line multigraph**, with colours, loops, and multiplicities retained.
It is not the single-arc-reversal graph.

**PROVED:** arbitrary-context R4 dissolution; ambient black connectivity
of the entire `{1,C3,R4}` SCC language at fixed vertex count and feedback
number; finite strong-centre inventories at every fixed positive black
deficit; rationality of every fixed-deficit self-converse stratum.
The incoming, independently audited theorem gives exactly `D_2(z)=z^4`.
**OPEN:** black-core connectivity for arbitrary tournaments, and the
existence or exclusion of a genuine mixed-to-mixed blue bridge.

**Subsequent continuation:** [spectrum_basins_20261020e](../spectrum_basins_20261020e/README.md)
connects all feedback inventories in this golden language, normalizes every
strong five-vertex end atom with a source context, and supplies black
detours around every blue line at a nonstrong mixed source. The remaining
blue-bridge question has two strong endpoints. The fixed-inventory results
below remain valid components of the stronger proof.

## 1. Connecting terminal basins by actual black paths

Write `A -> B` for an ordered SCC join and R4 for the unique strong
four-vertex tournament. A tournament is black-active if it has a black
Hamiltonian-path orbit. The [interior theorem](arithmetic_interior.md)
proves, for arbitrary nonempty tournaments A and W,

```text
A -> R4 -> W  <->  A -> C3 -> 1 -> W
```

by four black moves whenever the right endpoint is active. Converse gives
the alternative `R4 -> 1 -> C3`. Every internal context arc and every
final cross arc is restored literally. The proof chooses a starting
vertex in W's first SCC to exclude the only remaining potentially blue
intermediate step, then checks compatibility with self-converse endpoints.
It has no two-singleton guard and no bound on the outside SCC sizes.

The [end dissolution](arithmetic.md), [empty-context repair](moves.md),
Section10, and [prefix interchange](prefix.md) settle the boundaries.
Together they prove the following all-order result.

**Fixed-inventory normalization theorem.** All active ordered SCC words
over `{1,C3,R4}` with n vertices and feedback-arc number f lie in one
component of the ambient converse-merged black graph. They reach

```text
1^(n-3f) -> C3^f,       1 <= f <= floor((n-1)/3).
```

Here `f=#C3+#R4`. Dissolve the R4 atoms, then sort the rigid singleton/C3
word by adjacent swaps. Rigid pure-blue words are exactly palindromes;
every such word has inversion number at the midpoint. Sorting away from
that midpoint avoids them, and converse identifies the two sorted ends.
This supplies paths between distinct strict triangle minima, whose
intermediate tournaments can leave the SCC language.

This does not assert that f is preserved along every intermediate move,
that different f values are different black components, or that every
terminal tournament has only these three SCC types. Those would be
additional statements. The inherited
[all-order descent theorem](../black_descent_20261020c/cycles.md) still
reduces the full connectivity problem to arbitrary nonstrong terminals.

## 2. Blue bridges: new exclusions and the precise leaf convention

[Black bypasses](moves.md) now include:

* every blue line whose paired outer strong SCC has order four or five;
* every blue line at a mixed `1 -> C -> 1`, at **every order**, by at most
  five black moves;
* Hamiltonian-connected outside components, with the cited 4-strong
  sufficient condition.

The singleton-padding theorem combines an elementary regular-tournament
cut bound, score-changing Hamiltonian cycles, and a two-vertex-cut
argument excluding the even alternating-score obstruction. It uses two
explicitly CITED theorems of Thomassen; the original full text could not
be retrieved, and no audit of that original proof is claimed. The exact
imported statements and access limit are stated in moves.md, Sections9,11.
The [independent residual audit](moves_residual_audit.md) checks the
unbounded cut proof and reconstructs the complete small blue universe
without importing the producer code.
Complete fixed-blue-path checks and the uniform reductions cover
**every nonstrong mixed source through total order eleven**.
This is not a full order-eleven metagraph census.

There is also an all-order exclusion at every mixed degree-one vertex
of the **paired-black multigraph**, obtained by dividing each literal
black-line multiplicity by two. Such vertices have beta2 and are exactly

```text
M_W = W -> R4 -> W^op,
```

for rigid singleton/C3 words W. Every incident mixed-to-mixed blue line
has a black bypass. The [centre analysis](centers.md), Section6, gives
their unique black neighbours Q_W explicitly. For nonempty W, Q_W is
strong, non-self-converse, and has trivial automorphism group. An
intrinsic score-shell decoder recovers W, so different mixed leaves
have different pure-black attachment nodes. At empty W the neighbour
is the merged type `{1 -> C3, C3 -> 1}`.

Degree one in the paired multigraph is stronger than degree one in its
simple support: the latter can have several parallel black pairs.
No classification of all simple-support leaves is claimed. The
correction is recorded in [MISTAKES.md](../../../01-canon/MISTAKES.md).

## 3. D(z)=z^4 and what higher powers actually mean

For a self-converse T set `h=H(T)/|Aut(T)|`, let b be its blue path-orbit
mass, and put `beta=h-b`. Define D_beta to count **strong self-converse
types** of that deficit by vertex count. The incoming
[complete deficit-two classification](../deficit_two_centers_20261020/README.md)
proves

```text
D_2(z)=z^4.
```

The coefficient at degree four is one; all others vanish. We read its
structural cutoff and independent audits, then replayed the independent
affine reconstruction and every one of its **2,121,984 residual exclusion
certificates**, with zero survivors. These are necessary periodic
presentations, not 2,121,984 isomorphism classes. The retained
[receipt](incoming_replay.json) identifies the exact source hashes.
This is a rerun of the independent verifier, not a third implementation.

The new [general cutoff](general.md) proves that every fixed positive
deficit has a finite strong-centre polynomial. For even beta>=2, let L
be the least power of two strictly greater than beta+2. Then

```text
deg D_beta <= (beta+2)(2L-1) < 4(beta+2)^2.
```

Odd beta is impossible. The mechanism compares a power-of-two
cycle-rerouting clique against the connectivity forced by automorphic
rotation orbits. This gives a terminating coefficient algorithm, not an
evaluation of all higher-deficit coefficients.

The [SCC stratum theorem](strata.md) goes further. If B_beta counts all
self-converse types at fixed positive beta, its generating function is
rational, with denominator dividing a power of `1-z^2-z^6`. The first
deficit at which a pole of order j+1 occurs is exactly

```text
beta_j=3^j(3^j-1).
```

| deficit | pole order | coefficient prefactor on each parity |
|---:|---:|---|
|6|2|linear|
|72|3|quadratic|
|702|4|cubic|
|6480|5|quartic|

These powers count independent rigid gaps around nonrigid strong SCCs.
They are not additive powers of beta. For example, one through four R4
atoms have deficits `2,20,110,600`. The exact SCC law is

```text
beta(A -> Z -> A^op) = q^2 c - q d,
q=h(A), c=h(Z), d=b(Z).
```

The complete symmetry-slack-three centre sector has explicit Fibonacci
polynomial coefficients and degrees in [centers.md](centers.md),
Sections3–4. This is a second, separately graded family.
The previously defined metagraph degree polynomial D_n and maximum M_n
are different objects; their inherited status is in the
[preceding classification ledger](../black_descent_20261020c/README.md).

## 4. The golden ratio and 105 have actual maps

For ordinary SCC words over `{1,C3,R4}`, the
[arithmetic packet](arithmetic.md) proves

```text
F(z)=1/(1-z-z^3-z^4)=1/[(1+z^2)(1-z-z^2)],
[z^n]F(z) ~ phi^(n+2)/5.
```

The Fibonacci prefix code is `1 <-> 1`, `21 <-> C3`, `22 <-> R4`,
with a possible final unpaired2 retained as a flag. It preserves vertex
weight and, after retaining that flag, the directed-triangle grading.
Thus the golden growth belongs to an explicit tournament language on
which the new normalization theorem operates.

The recovered 105 character atlas identifies the golden character sqrt5
on `F2^3`. Its parallel matching on radical labels is

```text
{1,5}, {3,15}, {7,35}, {21,105}.
```

Under the previously proved marked matching-to-junction incidence, it
gives the three Fano lines through that character:
`{3,5,15}`, `{5,7,35}`, `{5,21,105}`. Across all seven directions this
selects a 21-dimensional bank of actual junction probes. The bank still
annihilates all seven Kempe blind coordinates. It therefore does not
complete the obstruction-preserving Petersen reduction.

The older golden-unit clocks 80,105,210 were replayed on their distinct
quotient rings. They are not interchangeable periods. Fibonacci Euclid
parameters `(13,8)` give the primitive triple `(105,208,233)`, but their
limiting leg ratio is1:2. The isosceles limit instead requires Euclid
parameter ratio `1+sqrt2`. These are explicit surviving connections and
explicit boundaries of the original analogy.

## 5. Incoming work integrated after the first proofs

The checkpoint also incorporates commits 760a29ae33 and a852aca417.
We read the
[endpoint-ear proof](../deficiency_three_20261020/black.md),
[larger-ear repair](../deficiency_three_20261020/ear_probe.md), and
[physical junction certificate](../deficiency_three_20261020/boundary.md)
before using them.

The incoming R4 construction improves every paired-R4 blue bypass to
**three** black path reversals. Consequently the six-move end-dissolution
bound in our original proof improves to four. To transport its exact
labelled certificate to our K convention, replace three reversal words
`(Q1,Q2,Q3)` by `(Q1,reverse(Q2),Q3)` as K words. Converse after the odd
steps gives exactly the original K target; path validity and blackness
are preserved. The existing six-step receipts remain valid witnesses.
The larger-ear theorem adds arbitrary-order outer SCCs with specified
closing/chord arcs and a nonautomorphic cyclic rotation; these conditions
retain the boundary data lost by a bare SCC-size invariant.

The physical junction theorem uses fourteen zero closure measurements to
force the target support, then Kempe switching forces a uniform tensor
`mJ`; one normalization fixes m. It certifies this particular tensor
without reconstructing the entire 91-coordinate space. A separate rooted
minor certificate remains necessary. This provides a concrete positive
direction beyond the seven-dimensional measurement blindspace.

We then tested whether our golden 21-probe bank supplies that certificate.
The [complete comparison](arithmetic_junction_certificate.md) separates
two predicates. Its full measurement vectors distinguish all 105
junction labels, but its zero rows leave 24 unwanted colour orbits for
each of the 21 bank targets and 30 for each of the other 84 targets.
Thus this particular zero-support proof fails. The incoming fourteen
target-adapted probes pass independently for all 105 targets. We have
not inferred a physical counterexample or failure of all twenty-one
measurement equalities from the weaker support-cover failure.

The incoming synchronized-triangle comparison also keeps an important
hostile: the two centres have the same feedback count, moved support,
automorphism order and Hamiltonian count, but different blue counts.
The [retained relative-orientation coordinate](../prime_batch_triangle_20261020/CONNECTIONS.md)
must survive any transfer to the fixed-deficit series.

## 6. Exact controls, reproduction, and next proof obligations

All new finite conclusions retain their source, universe, hostile
controls, and JSON receipt. Checkers use raising conditions that remain
active under `python -O`.

| checker | exact universe or consequence |
|---|---|
|`prefix_macro.py`|All 33,868 fixed-path prefixes through order7; ordinary-class controls through8|
|`moves.py`, `moves_five.py`|Uniform boundary ledgers, complete small outer presentations, literal black detours|
|`moves_endpoint_probe.py`|Every fixed-blue-path centre through order9; identifies fixed-start failures|
|`moves_residual.py`|All 3,21,537 strong two-score blue presentations at orders4,6,8 have explicit unbalanced cycle edges|
|`moves_residual_audit.py`|Independent universe reconstruction, exact certificate-set equality, and all 561 cycle replays|
|`moves_leaf.py`, `moves_end.py`|All central R4 leaf cases; 1,101 sink contexts through order7|
|`arithmetic_dissolution.py`|1,100 end-prefix presentations; 1,096 active targets|
|`arithmetic_interior.py`|4,400 context pairs; 4,391 active targets, including blue-intermediate hostiles|
|`arithmetic_wordgraph.py`|All 7,920 golden words through order18; 3,957 active merged nodes, 45 fixed-inventory strata|
|`arithmetic_normalization_audit.py`|Independent implementation: all 142,129 words through order24; 84 strata; 2,038,158 active rigid swaps|
|`centers.py`|Complete inherited low-deficit census through order9; 57 slack-sector words; 1,872 leaf decodings through order40|
|`strata.py`|Two independent word-series enumerations through size30 and mass81 on a declared restricted alphabet|
|`arithmetic_probe.py`|Fibonacci codec, actual tournament controls, exact closure-bank ranks and arithmetic clocks|
|`arithmetic_junction_certificate.py`|All 105 target labels and 2,187 colour words; golden-bank support failures and successful target-adapted covers|
|`incoming_replay.py`|Independent reconstruction plus all 2,121,984 incoming residual certificates|

Run a listed checker from the repository root, for example:

```text
python3 05-knowledge/results/scc_basins_20261020d/arithmetic_interior.py
python3 -O 05-knowledge/results/scc_basins_20261020d/arithmetic_interior.py
```

The macrograph checker tests edges represented by proved black paths,
not direct metagraph adjacency. The
[independent normalization audit](arithmetic_normalization_audit.md)
checks the proof's simultaneous path choices and all boundary cases;
its finite code imports no producer implementation.
The optional `interior_blue_probe.py`
is an affine discovery survey for nonempty context sizes `1<=m<=6`,
`1<=k<=14`; its finite absences are not the universal proof.

The continuation linked above resolves different feedback inventories and
arbitrary outer-SCC blue boundary cases. Remaining connectivity obligations
are normalization of arbitrary terminal SCCs and blue lines with two strong
endpoints. The higher-deficit
counting problem now has finite unknown strong-atom polynomials, but
counting alone still loses destinations of Hamiltonian-path complements.
No genuine mixed-to-mixed blue bridge or global connectivity proof has
been obtained.
