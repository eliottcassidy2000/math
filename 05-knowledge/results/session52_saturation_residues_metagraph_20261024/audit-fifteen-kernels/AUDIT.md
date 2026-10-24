# Independent audit: the fifteen-kernel channel

**Verdict: PASS within the stated finite and cited geometric scopes.**
Session 52, 2026-10-24. No publication or new research lane was undertaken.

The producer's elementary proofs, exact actions, inheritance credits, and
arithmetic stopping boundary are sound. Its geometric quotient descriptions
remain cited inputs, not consequences of the finite checker. No claim about
the MRB constant or Ogg primes is promoted by this audit.

## 1. Exact frozen candidate

The audit consumed precisely the 13 files in the producer's
`delivery-manifest.json`, plus that manifest itself. The frozen copy is in
`frozen_candidate/`. Every size and SHA-256 was verified before copying;
there are no symlinks and no undeclared source files in the copied tree.

- Candidate manifest SHA-256:
  `af8794a44f63d6251679a0b89d3228cb4c7d7407d9597421cd367a5fa24573c7`
- Candidate README SHA-256:
  `747ab7952914b4a9002ecc33ac98d6714596f9f1eb7cfc1cc9ace240930a8c07`
- Candidate checker SHA-256:
  `b2aa92ae33c2351b2fae36e24f291edfdbb7f5e0a8a4d672d8e6eac25a85b248`
- Replayed receipt SHA-256:
  `079757c20937b18f90d36f9917f88f310e3dce93307f4e195270787c5480dc73`

The producer's 28,604 checks pass from an actual manifest-only relocated
copy, under both normal and optimized isolated Python (`-I` and `-I -O`).
The working directory is `/tmp`, not the source directory. Both stdout
streams and both receipt files equal the frozen original byte for byte.
The original producer tree and frozen candidate remain unchanged.

## 2. Independent reconstruction, not an import of the producer

`independent_check.py` uses only the Python standard library. It does not
load producer code or generated producer data. Its 545,153 active gates
give byte-identical normal and optimized results.

The independent coordinate system differs from the producer's standard
symplectic basis. An even six-bit subset is replaced, if necessary, by its
complement so that bit 5 is zero. Its first four bits are coordinates;
bit 4 is their parity. The induced form is

    B(x,y) = parity(x & y) + parity(x)*parity(y), over F2.

The checker verifies bilinearity, alternation, and nondegeneracy. It
enumerates isotropic planes as additive three-element subsets of the
fifteen nonzero vectors. Independently it enumerates matchings as
three-element subsets of the fifteen edges of K6 with vertex degrees one.
The exact map from duads to binary vectors identifies these universes.
There are fifteen Lagrangians, three points on each line, three lines
through each point, and 45 flags.

For the matching graph, adjacency means sharing a duad. All vertex-pair
tests give strong-regularity parameters (15,6,1,3). This is also justified
by the producer's proof: adjacent matchings have the third matching on
their common edge as their sole common neighbor; disjoint matchings form
a six-cycle and have three common neighbors. The matrix identity and
spectrum in the candidate follow.

All six one-factorizations are enumerated independently. Every matching
is contained in two; sending it to that pair gives every two-subset of
the pentad set exactly once. Sharing an original duad is equivalent to
disjointness of the pentad pairs. Thus the displayed KG(6,2) isomorphism is
an actual relation-preserving map, not a comparison of parameter lists.
The producer's triangular-prism proof of unique extension of two
edge-disjoint matchings is valid and proves completeness without relying
on enumeration.

## 3. Action and stabilizer boundary

The independent checker verifies every one of the 720 relabellings and
all 518,400 multiplication pairs in the pentad action. It is faithful
and onto the permutations of the six pentads. A transposition has pentad
cycle type (2,2,2), so after choosing labels for the pentads the resulting
automorphism is outer.

The ordinary actions on duads and synthemes do not become equivalent.
A transposition fixes seven duads and three synthemes; a triple
transposition fixes three duads and seven synthemes. Independently,
both stabilizers have order 48, but the duad stabilizer has original-point
orbits of sizes 2 and 4, while the syntheme stabilizer is transitive.
They cannot be conjugate in the natural S6. Equality of stabilizer orders
does not repair the obstruction.

Consequently the graph isomorphism is equivariant for the induced
pentad action, not for the ordinary duad action on the original six
labels. The candidate states this correctly. This is not an isomorphism
from the syntheme graph to a graph of abelian-surface isogenies.

## 4. Source and characteristic audit

Primary sources were read directly, not inferred from search snippets:

- [Florit--Smith, Atlas, introduction, §§2.1, 3.1, 4.14 and Table 12](https://arxiv.org/html/2101.00917v2).
- [Florit--Smith, Automorphisms and isogeny graphs, Definitions 2.1--2.4 and Lemma 3.2](https://arxiv.org/pdf/2101.00919).
- [Howard--Millson--Snowden--Vakil, §1.4](https://arxiv.org/html/0710.5916v1).

The Atlas excludes characteristics 2,3,5. Its three displayed permutations
match the frozen transcription exactly. Independent closure gives S4,
order 24, orbits {1,4,8,9,12,13}, {2,7,14,15}, {3,6,10,11}, {5}, and
stabilizer orders 4,6,6,24. Table 12 supplies the codomain information:
the six-orbit reaches the elliptic square with j=8000; the singleton is
a loop; generically the four-orbits reach Type IV. At p=7 the four-orbits
are loops; at p=13,29 one target becomes Type V. Lemma 3.2 supports the
stated automorphism/dual-edge ratio law. HMSV supplies the classical
pentad construction, which is not new here.

Full geometric automorphism claims are interpreted over an algebraic
closure, or a sufficiently split source field, not an arbitrary unsplit
ground field. A subgroup action alone cannot certify the full group.

The following additional audit calculations require no geometric citation:

    8000 = 2^6 * 5^3,
    8000 - 1728 = 2^7 * 7^2.

Thus in odd characteristic the collisions are exactly p=5 with j=0 and
p=7 with j=1728. The first lies outside the cited Atlas scope. Neither
the abstract permutation closure nor this arithmetic extends the full
geometric classification to p=5.

At p=7, the cited table and orbit partition give a particularly decisive
endpoint counterexample: the four edge classes have only two distinct
destination vertices, the source itself and the product vertex. Three
distinct loop classes have weights 4,4,1. Their combined endpoint weight
is nine; this does not merge them into one isogeny class. Over an unsplit
field, kernel rationality is an additional restriction. Weighted degree
fifteen is not an unconditional simple degree over every base field.

## 5. Inheritance and mathematical scope

The frozen `cardinality_bridges_20261018m.md`, §3, already proves the
matching/Lagrangian carrier, incidence counts, and fixed-point obstruction.
The candidate credits these as inheritance. The older six-point observer
note already uses pentads and the outer action. The current packet's
contribution is a self-contained, executable reconstruction, explicit
relation transport, and a correctly scoped use of the published Bolza
kernel orbit data; it must not be advertised as discovery of classical
finite geometry.

THM-4012's frozen theorem retains its characteristic-zero setting and
the conditional face-stable specialization hypothesis. The occurrence
of the same Bolza Jacobian and j=8000 factor is a genuine identified
object shared by the two discussions. It does not remove the stable-model
hypothesis or prove a new obstruction for an MRB quantity. The packet
supplies no map from Ogg primes to kernels, no preserved predicate on such
a hypothetical map, and no privileged cyclic ordering or residue modulus
fifteen. These negative boundaries are explicit in both prose and receipt.

## 6. Hostile controls and reproducibility

`verify_replay.py` rejects all four mutations in both normal and optimized
Python, without emitting a success receipt:

1. Append a byte to the cited action input without changing its digest.
2. Remove that required input.
3. Change the claimed orbit partition and reseal its input digest.
4. Introduce overlapping cycles in a generator and reseal its input digest.

The last two distinguish substantive validation from a mere hash check.
The independent mathematical hostiles are the natural-action fixed-point
mismatch, equal-order nonconjugate stabilizers, and the separation of
kernel count, edge class, and endpoint count.

From this audit directory:

    python3 independent_check.py > independent_normal.json
    python3 -O independent_check.py > independent_optimized.json
    cmp independent_normal.json independent_optimized.json
    python3 verify_replay.py > replay_audit.json

On subsequent runs, `verify_replay.py` uses the frozen candidate and still
creates a new actual temporary relocation. It checks the exact manifest
file set before executing. The frozen source files are never mutated.

No historical checker absent from the sparse checkout was claimed to have
been rerun. No quotient curve, finite-field isogeny, or full automorphism
group was recomputed geometrically by this audit.
