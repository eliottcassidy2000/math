# Independent order-9 quotient census and parity-filter audit

**Status: FINITE-EXACT / VERIFIED.** This packet proves the supplied order-9
converse-merged census complete and independently verifies its prime/parity
filter. It does not, by itself, prove a response-polynomial exclusion or a
minimum-size theorem.

## Result

The independent audit establishes:

- 97,144 distinct converse-merged tournament classes;
- 191,536 distinct ordinary isomorphism classes, including 2,752 self-converse
  classes;
- exact labeled orbit mass **68,719,476,736 = 2^36**;
- 69,445 prime merged classes, representing 136,655 ordinary prime classes;
- 67,097 prime merged classes with some nonzero endpoint parity
  `p_i = (S_i + E_i) mod 2`;
- **2,348 prime merged classes with every p_i = 0**, representing 4,293
  ordinary classes, including 403 self-converse classes;
- exact agreement with all 2,348 candidate IDs, source codes, nine adjacency
  rows, nine starting counts and nine ending counts. This checks 42,264
  candidate endpoint entries. In total, literal path enumeration computes
  1,250,010 endpoint entries on all 69,445 prime representatives.

The independent program enumerates 33,867,002 degree-preserving vertex orders
across representatives and converses at order 9. Literal endpoint enumeration
visits 778,459,859 directed path prefixes on order-9 prime representatives.
The program runs only through order 9.

`summary.tsv` gives corresponding coverage checks for orders 2 through 9.
For n=2, the `prime_*` columns mean only absence of a proper nontrivial module;
this is the module-free convention, not a claim that the order-2 strong-tree
node is a prime node. The order-9 statement has no such convention issue.

## Independent completeness certificate

The certificate does not infer completeness from familiar class counts and does
not rely on the producer's color refinement or inductive generator.

For a labeled tournament T, partition vertices by outdegree, with cells sorted
by increasing degree. The auditor enumerates **every permutation within every
cell**, with no color refinement or prefix pruning. It encodes each resulting
adjacency matrix in row-major pair order `(i,j)` with `i<j`, unlike the source's
column-major pair order. Its canonical key is the minimum resulting integer.
An explicit check verifies that the number of enumerated orders equals the
product of the factorials of the degree-cell sizes.

Why this is a canonical form:

1. Every isomorphism preserves outdegree. It therefore maps each sorted degree
   cell onto the corresponding cell, and induces a bijection between all
   degree-respecting vertex orders. Isomorphic tournaments have identical key
   sets and identical minima.
2. Equality of two resulting full adjacency codes gives an explicit
   isomorphism through their vertex orders. Thus equal minimum keys cannot
   conflate two nonisomorphic tournaments.
3. Fix one vertex order attaining the minimum. Composition with each
   automorphism gives another order attaining the same code. Conversely, two
   orders with the same code differ by an automorphism. Since automorphisms
   preserve the degree cells, **minimum multiplicity is exactly |Aut(T)|**.

The auditor canonicalizes each input T and its converse T*, checks their
independently computed automorphism counts agree, and inserts the keys into an
exact ordered map. Every insertion is required to be new, except that T=T*
up to isomorphism is intentionally inserted once. Thus all 191,536 ordinary
classes are pairwise distinct, and every merged class is exactly one ordinary
class or one converse pair.

By orbit-stabilizer, a class contributes exactly `9! / |Aut(T)|` labeled
9-vertex tournaments. The auditor checks integral divisibility and sums these
orbit sizes over the distinct ordinary classes. The sum is exactly `2^36`,
which is the number of all labeled tournaments on nine vertices. These are
pairwise disjoint subsets of that finite universe whose total cardinality is
the cardinality of the universe. Therefore **none is missing**.

This establishes coverage without a separate re-generation of every extension.
Repeating the source's generator is unnecessary for this proof. The same
orbit-mass argument succeeds independently at every order from 2 through 8.

## Independent prime filter

A module M is a vertex subset such that every vertex outside M has a uniform
orientation toward all vertices of M. A nontrivial proper module has size
between 2 and n-1.

Instead of using the producer's exhaustive-subset loop, the main auditor starts
with each unordered pair and repeatedly adjoins any outside vertex that has
both an in-neighbor and an out-neighbor in the current set. Every module
containing the initial pair must contain each adjoined distinguisher, by
induction. If the closure stabilizes before reaching the whole vertex set, it
is itself a nontrivial proper module. If a nontrivial proper module exists,
starting with any pair inside it can never escape it. Consequently T is prime
exactly when every pair closure reaches the full vertex set.

Auxiliary controls also compare this closure algorithm against an independently
written literal subset/uniform-orientation predicate on all 97,144 order-9
representatives and all 33,866 labeled tournaments of orders 2 through 6.
There are no discrepancies. Another 1,002 deterministic shuffled-label tests
verify invariance of the independent canonical key and automorphism count.

## Independent endpoint filter

The main auditor uses no subset dynamic program for endpoint counts. Starting
at each vertex, it recursively extends every directed simple path along each
unused out-neighbor. At a path covering all vertices it increments the start
and end counters. Each directed Hamiltonian path has a unique sequence of
extensions, so it is counted exactly once. The nine starting and ending
counts yield every p_i directly.

All prime input representatives are checked, including ones discarded by the
producer's p-nonzero filter. The sum of their starting counts equals the sum of
their ending counts and agrees with the producer's H column. Producer tiling
counts are also checked on prime representatives. No approximation, modular
reconstruction, random endpoint sampling, or floating-point arithmetic occurs.

Primality is invariant under isomorphism and converse. Isomorphism permutes
endpoint coordinates; converse swaps starting and ending counts. Therefore the
predicate `p_i = 0 for every i` is invariant under both equivalences. The 2,348
merged candidates are a complete classification of prime order-9 quotients
satisfying this predicate, up to isomorphism and converse.

## Audit of the original source

The source `session189-deletion-incoming/source.c` has Git blob hash
`98953760287aa68650fb91728adea9686ea03e6d`, matching the recorded pin at
commit `e58a09c18bc7536e2388c85df239afd9d10a9c1d`.

Its coverage argument is mathematically sound in the audited n<=9 range:

- The initial vertex colors are the invariant pair (outdegree, number of
  directed triangles through the vertex), encoded as `1000*score+c3`.
  The c3 term is far below 1000 at these orders, so keys do not conflate
  distinct pairs.
- Each refinement signature consists of the prior color and the number of
  out-neighbors in each prior color class. Signatures and their deterministic
  ordering are isomorphism-invariant; inclusion of the prior color ensures
  cells only split. At most n-1 strict splits can occur. The use of `memcmp`
  on native integer signatures is not a portable claim of numerical
  lexicographic order, but on the execution platform it is a deterministic
  invariant order, which is sufficient for canonicalization. All signature
  entries are small nonnegative integers.
- The recursion enumerates every vertex order respecting the resulting
  ordered cells. Its prefix pruning compares already fixed bits against the
  corresponding prefix of the current best code; a strictly larger prefix
  cannot produce a smaller or equal full code. Equal prefixes are retained,
  so the minimum multiplicity still counts all automorphisms.
- Starting from the one-vertex tournament, `build_order(n)` uses **all ordinary
  classes**, not merely converse-merged classes, at order n-1, and all
  `2^(n-1)` new-vertex orientations for each. Deleting any chosen vertex from
  any n-vertex tournament gives one of those predecessor classes. Relabeling
  that predecessor to its stored representative determines one of the
  enumerated masks. Thus induction covers every isomorphism class.
- The hash map stores full canonical codes and compares exact keys. Linear
  probing resolves collisions; a hash collision never merges distinct codes.
  Valid codes have at most 36 bits, so they cannot equal the all-ones 64-bit
  empty sentinel. The observed ordinary class counts are below capacities
  2,048 (n<=7), 32,768 (n=8), and 1,048,576 (n=9). Thus probe tables have empty
  slots throughout these successful runs.
- Codes use at most 36 bits, masks use at most 9 bits, automorphism counts are
  at most 9!, and Hamiltonian counts are at most 9!. The relevant arithmetic
  is within the declared integer widths in this bounded run.
- Converse codes are looked up among ordinary classes and merged using their
  sorted positions. The independent audit verifies each self-converse flag
  and automorphism order, every merged class, and the final full orbit mass.

This is a mathematical/source audit plus independently replayable exact
computation, not a claim of formal verification of a C compiler or runtime.
Even if the source generator were distrusted, the independent disjoint-orbit
mass certificate is sufficient to establish the completeness of its data.

## Reproduction and scope

From the workspace root:

```
g++ -std=c++17 -O3 -Wall -Wextra -pedantic session189-c3-census-audit/audit.cpp -o /tmp/c3-census-audit
/tmp/c3-census-audit session189-deletion-incoming/data session189-c3-census-audit
g++ -std=c++17 -O3 -Wall -Wextra -pedantic session189-c3-census-audit/controls.cpp -o /tmp/c3-census-controls
/tmp/c3-census-controls session189-deletion-incoming/data
python3 session189-c3-census-audit/verify.py
python3 -O session189-c3-census-audit/verify.py
```

`verify.py` uses explicit exceptions, not assertions. Normal and optimized
Python produce identical receipts and output. `all_rows.tsv` retains every
independent canonical key, converse key, automorphism count, module witness or
prime classification, and all prime endpoint arrays. `candidates.tsv` retains
the independently reconstructed 2,348-candidate set. `receipt.json` records
input and output SHA-256 hashes. A self-contained compact replay needs only `audit.cpp`, `controls.cpp`,
`verify_compact.py`, the producer's `census.tsv` (node/code rows), and its
`candidates.json`:

```
mkdir -p /tmp/c3-compact-replay
/tmp/c3-census-audit --compact session189-c3-frontier/census.tsv /tmp/c3-compact-replay
/tmp/c3-census-controls --compact session189-c3-frontier/census.tsv
python3 session189-c3-census-audit/verify_compact.py session189-c3-frontier/census.tsv session189-c3-frontier/candidates.json /tmp/c3-compact-replay
python3 -O session189-c3-census-audit/verify_compact.py session189-c3-frontier/census.tsv session189-c3-frontier/candidates.json /tmp/c3-compact-replay
```

The compact mode computes its own self-converse flags and automorphism orders;
it has no source metadata dependency. Every order-9 output row, endpoint and
summary field agrees with the full tower-input mode. `verify_compact.py` also
independently recomputes the exact ordinary-key distinctness and orbit-mass
sum from the retained audit rows, and uses explicit exceptions in optimized
Python as well as normal Python.

The release allowlist explicitly names files;
there is no reliance on a broad workspace scan.

The certified universe is all prime order-9 quotient tournaments modulo
isomorphism and converse, followed by the all-p-zero predicate. Applying a
separate verified response exclusion to these candidates can justify a
finite exclusion for that quotient frontier. The derivation connecting this
frontier to a global minimum, the response exclusion, and all smaller
strong-tree cases must be supplied and audited separately.
