# Independent audit of the H=4151 equality classification

**Verdict:** PASS, subject to the source snapshot recorded below.

## All-order proof review

1. **Strict strong growth is valid.** A Hamiltonian cycle in a strong induced
   subtournament `U` has a `Y -> X` transition for every mixed exterior vertex
   `z`, where `z -> X` and `Y -> z`. Cutting the cycle there produces one path
   with two distinct insertions of `z`. Different old paths cannot collide
   after insertion because deletion recovers them. If there is no mixed
   vertex, `U` is a module and strongness supplies a two-vertex cyclic ear,
   yielding at least `3H(U)` paths. Thus proper strong enlargement is strict.
2. **Equality forces order ten.** The one-seven-orbit reduction gives a module.
   A rooted strong four-set in the quotient lifts to a strong ten-vertex
   subtournament with at least 4151 paths. At equality, strict strong growth
   excludes every larger strong ambient tournament. The exact substitution
   table leaves only the carousel and roots zero or two.
3. **Exactly two ordinary core classes.** Each core has automorphism group of
   order seven and a unique nontrivial vertex orbit, the seven-module.
   Isomorphisms must preserve that intrinsic orbit and therefore induce marked
   quotient isomorphisms. The quotient has trivial automorphism group, so
   roots zero and two cannot be identified.
4. **Converse pairing and no self-converse core.** The quotient's unique
   anti-automorphism exchanges roots zero and two. Combined with negation in
   the carousel, it pairs the cores. An anti-isomorphism fixing either core
   would induce a rooted quotient anti-isomorphism, which does not exist.
5. **SCC padding and counts.** The proved absence of H=7 and factorization of
   `4151=7*593` force exactly one nontrivial strong component. The others are
   singletons. There are `n-9` choices of the number of singleton components
   before a ten-vertex core, and two core choices, giving `2(n-9)`. Converse
   switches the core type and swaps the padding sides, so no padded class is
   self-converse, and there are `n-9` converse-merged classes.
6. **Automorphism and weighted corrections.** The padded classes all have
   automorphism order seven. All other H=4151 classes are rigid by the same
   equality classification and prime-order symmetry exclusion. Therefore
   `L_n/n! = N_n - 12(n-9)/7` and the path-pointed count is
   `4151 N_n - 7116(n-9)`.
7. **Minimum order for any symmetric 7p.** Other prime symmetry divisors are
   excluded independently of the count threshold. Below order ten a seven-
   automorphism has one nontrivial orbit, and its strong quotient has at most
   three vertices. The complete values `175,189,2163,2205` are not seven times
   a prime. This establishes the separate minimum-order assertion.

No inference above treats an induced Hamiltonian count as a divisor of the
ambient count. Divisibility is used only for whole strong components.

## Independent finite controls

`verify_equality_classification.py` was written separately from the source
verifier. It exhausts all 33,867 labelled tournaments of orders one through
six. On the strong tournaments it checks every proper induced strong subset
of size at least three, using a forward subset DP, and every possible rooted
strong-four containment.

- 398,488 strict strong inclusions checked.
- 136,736 rooted strong-four vertex instances checked.
- Strong labelled counts by order: `1,0,2,24,544,22320`.
- Both equality cores independently have `H=4151`, automorphism order seven,
  and strong connectivity.
- Complete score-preserving bijection enumeration proves the cores distinct,
  neither self-converse, and converse to each other.

The first core is

`[142,156,184,240,225,195,135,512,255,383]`.

The root-two core is

`[2,512,59,115,227,451,391,271,31,509]`.

An independently found isomorphism from the first core to the converse of
the second is `[2,8,7,6,5,4,3,9,0,1]`.

Reproduce with `python3 verify_equality_classification.py`. The matching
receipt is `equality_independent_receipt.json`. Normal and `-O` runs of this
independent verifier also produced byte-identical receipts. These finite tests are hostile
controls; the all-order classification is justified by the structural proof.

## Snapshot and source replay

The audited immutable source snapshot is:

- `session196-orbit-spectrum/PROOF.md`:
  `fafa5fb9224ae6fee2fcb91a65684717447c691ba6546c7180475b6ec2d25d8c`
- `session196-orbit-spectrum/verify.py`:
  `28aff6dc254415a5c79c9aada1041be6bb4b37c22daf46bbc9aa72f2cf2439c8`
- Replayed source receipt:
  `2c080e6f5142cdda2d36a2dfbddac95a37c7f4001a9ee36202781c2f68d7f3e3`

The frozen verifier was copied locally and run both normally and with `-O`.
Both executions passed 196,123 explicit checks, with byte-identical stdout
and receipt files. The matching receipt is preserved as
`orbit_equality_replay_receipt.json`.

The final appendix was reviewed in this snapshot. It removes the prior
spectrum-omission dependency by proving H7 absence and hypothetical H21
rigidity internally. No mathematical correction is outstanding.

Publication may change the single relative cross-reference from
`../prime-audit/PROOF.md` to `../symmetry-audit/PROOF.md` to match the parent's
actual directory layout. That link-only change does not alter this verdict.
