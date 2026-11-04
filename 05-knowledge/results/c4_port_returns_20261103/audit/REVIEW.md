# Independent whole-attachment and native-transport audit

2026-11-03. **FINITE-EXACT independent reproduction** of the fixed marginal
63 attachment classification through suffix order seven. **CITED** Busch
bound completes the all-order reduction. **PROOF REVIEW PASS** for the
separately pinned native padded-guard and fixed-prefix arguments. These
reviews do not promote arbitrary two-digit source closure or the global
Hamilton-path spectrum conjecture.

## 1. What was independently enumerated

The reviewed producer is `regeneration/small_bank.cpp`, with its `full` and
`complete7` receipts. It enumerates all labeled tournaments at orders one
through six, then takes all64 births from every one of56 order-six
isomorphism representatives. Its degree-respecting canonical key is exact:
an isomorphism preserves degrees, and every permutation within each ordered
degree block is tested. Every order-seven tournament is covered by deleting
one vertex, identifying the remaining tournament with a bank representative,
and recording the deleted vertex's incident cut. Every port subset is then
tested, so this quotient loses no marked-subset choice.

The independent `whole_attachment.cpp` deliberately uses different paths:

1. For orders one through six, scan all labeled codes and mark their entire
   permutation orbits. This gives 1,1,2,4,12,56 classes without any
   degree-partition canonicalization. Use all56*64 births at order seven.
2. Enumerate every vertex permutation of the suffix. It is a two-slot cover
   exactly when all consecutive arcs are directed except possibly the one
   split arc. A directed Hamilton word contributes every split, including
   the two empty-slot cases. A word with one bad consecutive arc contributes
   only its split there; two bad arcs contribute none. This independently
   recovers h, every first-slot endpoint count, and every ordered endpoint
   joint. No subset endpoint convolution is used.
3. For every actual port set, construct the full clock-one attachment and
   test forward and reverse reachability. This replaces the producer's
   terminal-SCC filter with actual full-graph strongness.
4. For every strong candidate with f=63, count Hamilton paths in its actual
   full attachments at clocks one and two. Check H_2-H_1=126 and
   2H_1-H_2=h-j. This is a direct full-graph check of both coefficients.

The complete finite universe and results are:

| Suffix order | Graphs tested | All port sets | Strong attachments | f=63 candidates | 22,24,42,44 hits |
|---:|---:|---:|---:|---:|---|
| 1 | 1 | 2 | 1 | 0 | 0,0,0,0 |
| 2 | 1 | 4 | 2 | 0 | 0,0,0,0 |
| 3 | 2 | 16 | 11 | 0 | 0,0,0,0 |
| 4 | 4 | 64 | 45 | 0 | 0,0,0,0 |
| 5 | 12 | 384 | 308 | 10 | 0,0,0,0 |
| 6 | 56 | 3,584 | 3,133 | 22 | 3,0,3,0 |
| 7 | 3,584 | 458,752 | 426,587 | 153 | 7,0,0,0 |

At order seven the extension bank retains isomorphic duplicates. The
185 fixed-slope candidates and their hit multiplicities are not counts of
distinct rooted isomorphism classes. Duplicates cannot invalidate the
absence result; they are retained deliberately to avoid relying on an
order-seven canonicalizer. The 22 and42 hits are positive controls.

## 2. Why the finite result excludes every fixed common-clock suffix

The exact attachment formula is

    H(M(TT_t;B,W)) = c+2^t f,  c=H(B)-J_W(B).

The two proposed second-letter +4 targets require (f,c)=(63,-35) or
(63,-29). At clock one they have H=91 or97. The actual full attachment is
strong exactly when W meets the terminal SCC of B; suffix B need not itself
be strong. This condition makes the same suffix attachment strong at every
nonempty prefix length, including one.

Busch's [official paper, Theorem1 and Corollary1, pp3-4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
was independently checked. Its minimum at strong order ten is125, and its
bound is at least125 at every larger order. Therefore either proposed
clock-one full attachment has order at most nine and suffix order at most
seven. The exhaustive computation above contains every possibility and
finds none. The bound is **CITED**, while the enumeration is **FINITE-EXACT**.

This argument also excludes the identity if requested only at sufficiently
large clocks: two distinct clocks force the same coefficients, and the
fixed actual suffix still defines the forbidden clock-one attachment. No
native operation at clock one is assumed. Conversely this does not exclude
other slopes yielding the isolated count91 or97, a fixed output-clock
shift, clock-dependent suffixes, or a different attachment model.

## 3. Native padded-guard proof review

Reviewed `native/PROOF.md` and the pairwise native interpreter in
`native/check_native.py`; this audit does not repeat its finite replay.
The padded gadget is sound. For an actual Hamilton guard P containing
a->b, terminal z and helper h with sole old out-neighbor a make

    P1=(P,z,h),  P2=P with h inserted between a,b, followed by z

legal successive guards. The two new adjacencies a->h->b are initially
oppositely oriented and are not consecutive in P1, so the first complement
reverses them. All other adjacencies of P2 were retained. The pad handles
the case where b was the old last vertex. On old pairs, the guards' edge
sets differ only at {a,b}; their two native complements therefore reverse
exactly that old pair. Every other old vertex has the same predecessor and
successor in both guards, including old endpoints.

Homogeneous transitive-prefix inflation is also valid. Replace the quotient
marker by its full transitive path. A first native move may reverse its
nonconsecutive internal arcs, but the explicitly lifted second guard uses
the retained consecutive arcs, and the second complement restores every
old prefix pair. The three displayed quotient guards were checked in
their actual successive incidence/port states. Their completed count and
root-defect formulas agree. Intermediate helper graphs need not be strong,
as the producer explicitly states.

## 4. Fixed-prefix arithmetic proof review

Reviewed `prefix/PROOF.md`, `prefix/probe.cpp`, and `prefix/check.py`.
For a fixed Q ordered before or after TT_(t-k), the necessary coefficients
are ac=-35 or-29 and pf=63*2^k, where a=H(Q),p=P_Q(2). Negative k is
impossible by integrality. SCC multiplicativity and the inherited all-order
absence of Hamilton count7 reduce each possible nontransitive Q to one
strong core of count5,29,or35 and singleton components. Busch's order-eight
minimum45 bounds these cores to order seven. The recorded complete core
inventory and its independent literal counts supply odd P2 parts
1; 5,145; and169, respectively.

Thus every nontransitive fixed Q is excluded for44. For24 only a=5 remains,
forcing c=-7 and f=63*2^(k-5-z), k>=5+z, where z counts singleton components.
This is a residual target, not a realization. The same scalar proof works
if the variable transitive block is inserted at a fixed ordered cut
Q1->TT_(t-k)->Q2, by taking Q=Q1->Q2. It does not cover non-ordered
incidences, clock-dependent Q, or changing suffixes.

This is a logical/code review of the prefix producer, not an independent
replay of its nine selected-core counts. Its own checker supplies those
permutation controls. The distinction is retained in the audit receipt.

## 5. Reproduction and pins

    python 05-knowledge/results/c4_port_returns_20261103/audit/check_audit.py
    python -O 05-knowledge/results/c4_port_returns_20261103/audit/check_audit.py

The wrapper compiles the independent C++ implementation in a temporary
directory. All guards use explicit failures, and both modes give identical
stdout and JSON. `receipt.json` records the exact universe and hashes of
reviewed source files. The previously pinned language packet is unchanged.
