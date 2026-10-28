# Complete parent-only construction through order ten

**PROVED:** the inherited augmentation mechanism, with the elementary shortcuts
below. **FINITE-EXACT:** this implementation and its complete saved order-ten
bank. **VERIFIED local timing:** not a matched comparison with another machine
or another output contract.

## Canonical root mechanism

Start with the empty parent. For every representative A of order n, compute a
canonical representative and **all** attaining transports. Conjugate these to
the concrete group G=Aut(A) in the saved representative. Enumerate one mask m
per G-orbit; a set bit means the new vertex points to that parent vertex. This
enumerates pointed (n+1)-tournaments exactly once.

Initialize ordered vertex cells by ascending outdegree. Repeatedly split each
cell by the vector of outgoing counts into the current ordered cells, keeping
the old cell order. This is an isomorphism-invariant ordered partition. Select
its last cell C. A unique maximum-degree vertex is automatically the sole
member of C: no subsequent split can change that singleton or its last position.
This permits the implemented early acceptance without completing refinement.

For every v in C define the **exact** pointed key

    K(T,v) = (canonical(T-v), minimum transported root mask).

The mask minimum uses every canonical transport of the deleted parent. Keep a
pointed extension exactly when its root belongs to C and has minimum K among C.
Two roots give the same key iff an automorphism of T carries one to the other:
an isomorphism of the deleted parents preserving the transported mask extends
over the roots, and restriction proves the converse. Therefore precisely one
pointed orbit survives per ordinary class. Ties count vertices in that orbit,
not extra ordinary classes. If k vertices share the winning key and the parent
mask stabilizer has order h, orbit-stabilizer gives |Aut(T)|=kh.

The selected cell is only a candidate filter. Distinct roots of a regular
tournament can share all refinement colours. No equivalence is inferred from
colours, degree sequences, cyclic triples, or group orders alone. The mechanism
is inherited from [root ports, CARRIERS](../root_ports_closed_carriers_20261028/CARRIERS.md),
which already reached order nine with 1,286 competing-root queries. This C++
implementation extends that exact mechanism to ten.

## Actual output and costs

| Order | Pointed candidates | Ordinary classes | Strong classes | Competing root queries |
|---:|---:|---:|---:|---:|
| 8 | 54,256 | 6,880 | 6,008 | 146 |
| 9 | 1,716,608 | 191,536 | 178,133 | 1,286 |
| 10 | 97,213,472 | 9,733,056 | 9,355,949 | 15,107 |

At ten, the maximum-degree filter rejects 81,528,886 pointed candidates;
stable-cell membership rejects another 5,951,080; exact key comparison rejects
450. There are 191,536 parent preparation searches and 15,107 competitor
searches, totaling 206,643. A runtime guard forbids production canonicalization
above parent order nine. This is **not** the full root-switch matrix Q: emitting
all its entries has a different output cost.

The accepted-order-ten automorphism histogram is:

    1:9677605, 3:54370, 5:274, 7:10, 9:738, 15:32,
    21:12, 25:1, 27:10, 63:2, 81:2.

The production count uses no supplied A000568 values. The publisher compares
with the classical counts only after generation and verification. The strong
count agrees with the already-existing [nauty order-ten census](../order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md).
That packet predates this complete bank; do not describe this as the first
order-ten census.

## Complete alternate audit

`verify_bank.cpp` decodes **every** saved order-ten representative. It independently
implements stable colour refinement and then enumerates every permutation
within the resulting colour cells. It does not use the producer's recursive
individualize/refine canonical search or its selected-root rule. Counting
attaining minimum orders recovers each full automorphism order. Sorting all
9,733,056 resulting canonical codes establishes uniqueness.

Finally it checks

    sum over saved T of 10!/|Aut(T)| = 35,184,372,088,832 = 2^45.

Every labelled tournament belongs to an isomorphism class of that size.
Uniqueness plus this exact mass proves completeness; matching the known scalar
count alone would not. The verifier shares primitive row encoding with the
producer, so independence is at the canonical-search level. A separate Python
literal-permutation audit checks those encodings and rooted fibres through five.
The original 81/27 examples also pass the inherited direct-isomorphism search,
Hamilton dynamic program, and complete labelled cyclic-triple comparison.

## Bank format and consecutive reuse

`parents9.u64.gz` contains 191,536 sorted canonical nine-vertex codes, eight
little-endian bytes each. Edge bits use lexicographic i<j order and bit 1 means
i points to j. An order-ten representative is reconstructed from a parent and
the nine-bit root-outmask. Its address is

    address = (parent_index << 9) | mask.

The address fits 27 bits. `flags10.varint.gz` stores sorted accepted addresses
as unsigned base-128 variable-length **differences**, first address measured
from zero. Gzip has deterministic mtime=0. `bank_manifest.json` records the raw
and compressed checksums, lengths and conventions. Total compressed size is
829,899 bytes, compared with 77,864,448 bytes for one u64 per representative.
The file is a sequential bank; random access needs a delta-stream index.

This saves repeated parent storage and makes parent/child pairs available
together. It is compression of a particular structured census, not an entropy
bound for arbitrary independently chosen tournaments. A tiny program can also
describe a huge census, while expanding the census still costs output time.
The bank does not persist all parent automorphisms; `generate.cpp` reconstructs
them exactly when needed. An attachment consumer must do the same or cache them.

For **independent aligned** consecutive pairs, use the concurrent
[Session89 double-coset decomposition](../session89_consecutive_tournament_pairs_20261028/generation/ALIGNED_GENERATION.md).
For fixed cores T,S it enumerates Aut(T)\S_n/Aut(S), then mask orbits under the
intersection stabilizer. For an overlapping stream with common middle object M,
the gluing set is A\Aut(M)/B, where A and B are the two embedded residual
stabilizers. This is the correct place to reuse a middle dictionary. Replacing
A or B by the full Aut(M) can erase earlier attachment constraints.

**Next experiment, OPEN:** augment the order-ten dictionary with lazily cached
full parent groups and selected-root transporters, then measure actual stream
queries. The present receipts establish the bank and its output costs, not a
measured order-ten independent-pair or arbitrary-stream implementation.
