# Independent audit of the marked native cut presentation

**Audit verdict: PASS for the exact reviewed proof and code below.** The
producer remains RESERVED pending root integration; this review does not
promote it or edit its mathematics. All finite claims retain their explicit
order limits. No source or artifact outside this audit file was edited by the
reviewer in the native lane.

| Reviewed file | SHA256 |
|---|---|
| `RESULT.md` | `344d0dec8efdce732c51a243fab11ccdbdb22905447cfffc08f5661077d6f2ba` |
| `check.py` | `c639decc73052fdb238d63911db2105c749476544af26bbc50fc7c0a73687550` |
| `normal.json` | `bd096be5452c49a286538112819913cd0e20278dc8c7e5d249b55b2e986cf1a3` |
| `optimized.json` | `bd096be5452c49a286538112819913cd0e20278dc8c7e5d249b55b2e986cf1a3` |

The immutable prior matrix input is
`incidence_closure_20261113/closure/data.json`, SHA256
`c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7`.
Base commit for this packet: `68d0d86bef480b6a36efc3d9864fd01ff3244bf5`.

## All-order proof review

**Affine action — PASS.** The proof's reversal mask and the code's
increasing-arc mask differ by the all-one vector. Its line-graph coboundary
vanishes, so both conventions give the same beta relation. The native move
complements nonpath arcs while retaining the ordered Hamilton path; its mask
change is `1+chi_P`. Applying delta therefore gives `delta chi_P` exactly.
The cut count is `sum_v d_P(v)(n−1−d_P(v))`. Substituting the two degree-one
endpoints and n−2 degree-two internal vertices gives `2(n−2)^2`. This is
positive at every n>=3, hence labelled beta states have no native loops.
Relabelling can identify the two endpoints later and create quotient loops.

**Exact adjacency guard — PASS.** For two realizable beta fields, their
difference is a cut and has precisely the two preimages s and its complement.
A native jump exists iff at least one is the support of a spanning path
consistently oriented in the reconstructed tournament. Connectivity, degrees,
and consistent direction are all needed. Reversing the tournament reverses
the path, preserving the predicate. The proof supplies necessity and
sufficiency, not just an invariant necessary condition.

**Order-four exception — PASS.** Equality of two different support jumps
forces their symmetric difference to be the whole complete graph, hence
complementarity. Equal path sizes then force n=4. The displayed two ordered
paths are complementary and assign every edge a consistent orientation.
They prove attainment of the double-support case. The statement is made at
fixed labels; quotienting targets by relabelling can combine additional
jumps at every order and is not covered by that uniqueness assertion.

**Marked carrier and source weight — PASS.** For n>=2 an ordered path
selects exactly one of a tournament and its converse. Thus `(beta,p)` is
bijective to the ordinary labelled pair `(T,p)`, equivariantly under vertex
permutations. That group action is free because a full ordered vertex list
has no nontrivial stabilizer. Normalizing p to the identity path leaves
exactly the free nonpath coordinates, giving `2^((n−1)(n−2)/2)` fixed-spine
states. The native involution retains p, and converse reverses p; these
operations commute. Reflection has no fixed ordered list before quotienting.
The revised text correctly locates its fixed words on `D_n/S_n`, the
fixed-spine carrier. These are precisely the inherited blue darts.

The number of reflection-invariant coordinate assignments is the power of
two in the proof: among the `(n−1)(n−2)/2` free coordinates, exactly
`floor((n−1)/2)` are individually fixed. Taking paired-coordinate orbits
therefore gives `floor((n−1)^2/4)` free reflection-orbit bits. Prematurely
forgetting path direction identifies every nonfixed reflection pair and
changes integer dart mass. In particular the order-four and order-five
losses 8 to 6 and 64 to 40 are correct.

On an ordinary class the path-marking action of its automorphism group is
free, giving H/a darts. Combining one or two ordinary classes in a converse
fibre gives exactly mH/a. The native complement has no fixed fixed-spine
word for n>=3, so quotient loops still have two native darts. Blue and black
weights are retained separately; they are not counts of unordered marked
paths. This is the crucial distinction from merely specifying which merged
nodes are adjacent.

**Scope — PASS.** The presentation gives exact labelled states, legal moves,
markings and quotient multiplicities. It does not assert a smaller
changing-order enumerator, connectivity of the native graph, or a Hamiltonian
spectrum realization theorem. The full wedge code's minimum distance and the
native jump weight concern different allowed cut supports. Their common
linear carrier does not make their transition graphs identical.

## Code audit and exact replay

The implementation discovers every ordinary class by the least-unused
labelled tournament and its complete permutation orbit; it does not read a
prior orbit map. It separately enumerates fixed-spine states and all physical
Hamilton paths of each independently found representative. The four ordinary
and merged uncolored/blue matrices agree entrywise with the immutable prior
census. Source automorphism division is tested separately on physical path
counts, including blue counts. Thus matching the prior totals is not its only
normalization check.

`cut_inverse` propagates edge bits and checks every constraint; the path
support test checks degree, a connected traversal of all vertices, and either
possible path direction. Every actual Hamilton path through six passes the
exact target guard, and every possible target beta is tested through five.
The source clearly distinguishes those two finite universes. Counterexamples
check a disconnected triangle-plus-edge degree pattern and a spanning path
with inconsistent orientation. There are no removable assertions.

The reviewer independently reran both normal and optimized executions. Both
produce the frozen receipt above, with **51,504 explicit checks**. This is
execution verification of the reviewed code, not an additional independent
algorithm for every matrix.

A separate small calculation, without importing the producer, enumerated
unordered spanning path edge sets directly from permutations and computed
their line-graph cuts as sets. The exact results are:

| n | Distinct undirected path supports | Distinct jumps | Every jump weight |
|---:|---:|---:|---:|
| 3 | 3 | 3 | 2 |
| 4 | 12 | 6 | 8 |
| 5 | 60 | 60 | 18 |
| 6 | 360 | 360 | 32 |

This independently checks the exceptional support multiplicity and the
constant weight, without relying on orientation or any prior metagraph file.
The algebraic proof is still what establishes the all-order statement.

During review the producer repaired the companion filename from
`reconstruction/RESULT.md` to `reconstruction/PROOF.md` and made the quotient
location of reflection fixed points explicit. These were reference/type
clarifications; no formula, matrix, or computed receipt changed. The exact
hash above includes both repairs.
