# Session67 incoming prime-family / C3 review

## Pins, credit, and boundary

- Reviewed integration commit **24f8e62c8f4cfbb54b80639b8bbd801ded9a5948** of `eliottcassidy2000/math`.
- Producer commit **62430ddb33e685494ae50f00d05e24c837bb3e02**, credited by the source to **klein, 2026-10-26**, worktree `codex/session-tournament-twist-unify`.
- Producer packet: `05-knowledge/results/a000568_intransitivity_families_20261026/`.
- Primary canonical records inspected and copied unchanged to `source-canonical/`: THM-4803 and HYP-9413. THM-4803 declares `CITED + PROVED + FINITE-EXACT`, explicitly **not independently audited**; HYP-9413 declares its score/no-cancellation statement OPEN beyond finite n<=7.
- THM-4803 SHA256: `fbf58f1d9702f625039643aad917ca4e423aed588afef4835722304c547495b5`.
- HYP-9413 SHA256: `92e1640941d3d47e1483409693980909a0f7fb712d4e325bafbd8fa5e3373e2b`.
- All **48/48** entries in the producer's SHA256SUMS match the pinned bytes (`source-integrity.json`). Integrity is not correctness.
- This review writes only `session67-incoming/`. Publisher-owned repository was read-only. No push, downloaded external binary, or large class census was run. No blocked action or permission uncertainty occurred.

## What can be integrated

**Exact source-data arithmetic independently replayed.** A fresh standard-library Fraction implementation evaluates the 23 aggregated prime cycle-index rows, taking source order-10 totals P10=7,906,161 and self-converse P10=5,525 as inputs. It recovers, for n=1..10:

- A000568: `1,1,2,4,12,56,456,6880,191536,9733056`.
- Strong tournaments: `1,0,1,1,6,35,353,6008,178133,9355949`.
- Unlabelled C3 structures: `1,1,2,3,9,32,247,3328,93612,4774417`.
- Strong C3 structures: `1,0,1,1,6,22,214,3074,90255,4680547`.

The two ordinary species equations are A=1/(1-S), S=x+sum_Q Z_Aut(Q)[A-1]. Their C3 analogues replace ordered strong components by a multiset and prime quotients by duality classes acted on by Aut±. This is an effective supplied-input evaluation of the Session56 identity, **not an independent prime census** and not an all-orders closed enumeration. Attribution for the cycle-index census and supplied n=10 totals stays with klein. The source's nauty census, n<=7 brute-force fibre checks, and n<=9 H/Aut finite agreement receipts were read, not independently rerun. The existing packet's normal/optimized arithmetic receipts are not independent implementations.

`independent_review.py` and `independent-review.{normal,optimized}.json` contain the independent arithmetic and explicit counterexamples; normal/-O agree byte for byte.

**Cited structural theorem verified at primary source.** [Boussaïri–Chergui–Ille–Zaidi (2018)](https://arxiv.org/html/1805.04607v1), Theorem 10, gives shared strong modules for a realizable labelled C3 hypergraph and every realization. Proposition 51 makes quotient choices at those modules independent; Theorem 24 restates the earlier Boussaïri–Ille–Lopez–Thomassé prime uniqueness-up-to-duality result. Corollary 57 gives N(T)=2^p product(k!). These are cited inputs, not Session67 discoveries. The paper's special hypergraph-module definition must not be confused with another hypergraph modular decomposition.

## Material corrections to incoming statements

### 1. Automorphism order is NOT an all-orders C3 invariant; first failure is order 12

Let A be the order-four tournament `singleton -> directed triangle`. Its converse A* is not isomorphic to A: their sorted scores are `(1,1,1,3)` and `(0,2,2,2)`. Set

    T = C3[A,A,A],    U = C3[A*,A,A]

on the same twelve labelled vertices, changing only the first block's internal arcs. Their C3 hypergraphs coincide: internal cyclic triples survive reversal; triples with two vertices in a block are transitive; transversal triples are unchanged. Both have **67** cyclic triples. Nevertheless

    |Aut(T)| = 3^3 * 3 = 81,    |Aut(U)| = 3^3 = 27.

The internal factor is three rotations per triangle. T admits the quotient's three rotations, whereas U's distinguished A* block kills the two nonidentity rotations. Complete degree-pruned arc-preserving permutation backtracking independently verifies 81 and 27. Exact Hamiltonian-path DP gives **H(T)=H(U)=462915**, so this counterexample does not challenge H-invariance. The full adjacency matrices are in `independent-review.json`; the ordinary labelled C3 fibre has N=128.

**Why the incoming proof fails:** checking a quotient dualization while its children are held fixed misses changes to a parent's stabilizer when a descendant's ordinary isomorphism type changes. Finite agreement through n=9/10 cannot establish the proposed all-order invariant.

**Sharp minimality proof (Session67 root/review synthesis; conditional on the cited structural theorem).** Induct on subtree order. At a prime node, Aut(Q) has odd order: an involution would exchange an arc's endpoints, which is impossible. Restrict to the subgroup preserving child sizes. For total order <=11, there are at most two children of size >=4; this odd-order subgroup fixes both individually. Children of sizes <=3 retain their ordinary tournament isomorphism type throughout their labelled C3 fibre. Hence the child-decoration stabilizer is unchanged. Child internal automorphism orders are unchanged inductively. A quotient dualization leaves the literal group Aut(Q) unchanged; at a linear node, reordering children preserves the product of child automorphism orders. Therefore |Aut| is C3-invariant for every tournament of order <=11, and the order-12 example is the first possible failure. This is an all-orders structural proof of the cutoff, **not exhaustive enumeration through 11**.

### 2. Strong-component ordinary isomorphism types are not C3-invariant

The invariant is the component **vertex partition / sizes**, and recursively the modular data modulo permitted moves, not the ordinary tournament isomorphism types. Already V=C3[A,1,1] and V* are strong order-six tournaments with equal labelled C3 but sorted scores `(1,2,2,2,4,4)` versus `(1,1,3,3,3,4)`. Independent isomorphism backtracking returns zero maps. The source analysis script only tested the **size multiset** (`compsizes`), not the stronger prose claim.

### 3. Additional local corrections / scope cautions

- Claimed monotone level trends are contradicted by FAMILIES.md's own n=9 table: mean log2 N rises 1.21->1.23 at c3=23->24 and 1.14->1.20 at 29->30; prime fraction drops .82->.81 at 23->24. Say broad tendency, not monotonicity.
- HYP-9413(d): order-10 cycle indices alone do not predict A000568(11). The coefficient also contains the unknown additive **P11**; C3 counting needs `(P11+P11,self-converse)/2`. A known A11 can instead infer P11 after lower-order contributions are computed.
- The C3 dualization bit is not generally the A_n relabelling-parity bit. For the regular C5 with i->j iff j-i mod5 is 1 or 2, the even reflection i->-i, `(1 4)(2 3)`, maps T to T*. Thus its two distinct labelled C3 realizations lie in the **same** A5 orbit. Independently checked in `independent_review.py`.
- INFORMATION.md's `L=id: binom(n,2) bits` mixes ambient capacity with fibre loss: H0(T|id)=0, whereas the labelled tournament universe has capacity binom(n,2).
- The binary, unsigned C3 hypergraph is **not** the full signed addressed triangle-curl field. Session67's scores-plus-signed-curl exact recovery is not contradicted by residual C3-plus-scores ambiguity.
- Do not promote claims about wild-pair completeness or the OEIS correction from this review: the full underlying independent census was not rerun here.

## A genuine strengthening: settle the score-fibre no-cancellation question

**Theorem (Session67 proof, with cited structural input).** Let b(T) be the number of prime nodes in the shared modular tree whose quotient is balanced with respect to child sizes: for every child i,

    sum_{i->j} |M_j| = sum_{j->i} |M_j|.

On the fixed labelled vertex set, the number of tournaments with both the same C3 hypergraph and the same **vertex-addressed score function** as T is exactly **2^b(T)**. Consequently the pair determines T precisely when b(T)=0. This upgrades HYP-9413(b)'s finite n<=7 evidence to a proof; the structural C3 theorem remains explicitly cited.

**Proof.** For any shared strong module M of size m, let d(v) be the full vertex score. Every exterior vertex sees M uniformly, so its external outdegree is

    D(M) = (sum_{v in M} d(v) - binom(m,2)) / m.

This is determined by the addressed scores, independent of orientations internal to M. If M_i is a child of M, then

    D(M_i) - D(M) = sum_{M_i -> M_j} |M_j|.

Thus the score function separately determines the weighted quotient outdegree at **every** node: descendant or ancestor moves cannot cancel one another. At a linear quotient these weighted outdegrees are distinct positive-weight suffix sums, recovering the order of all children. At a prime quotient the only choices are Q and Q*, and their weighted outdegrees coincide exactly in the size-balanced case. Each balanced reversal preserves every vertex's score. The cited independent node-choice bijection therefore leaves exactly two choices at each of b nodes and one at every other node. QED.

**Independent bounded checks.** `score_fibre_check.py` extracts strong modules directly by checking every subset. On **all 1,099 labelled tournaments of orders 1..5**, it groups by exact `(C3, score vector)` and verifies fibre size 2^b; zero mismatches. It also checks the ordinary C3 fibre formula. Additional generated tree-fibre controls:

| Tournament | order | b | full C3 fibre | C3+score fibre |
|---|---:|---:|---:|---:|
| C3[TT2,TT2,TT2] | 6 | 1 | 16 | 2 |
| C3[TT1,TT2,TT3] | 6 | 0 | 24 | 1 |
| C3[C3,C3,C3] | 9 | 4 | 16 | 16 |
| C3[A,A,A] | 12 | 4 | 128 | 16 |

Normal and -O receipts are byte-identical (`score-fibre.normal.json`, `score-fibre.optimized.json`). These four controls enumerate their explicit tree fibres, **not** all tournaments at orders 6, 9 or 12.

## Session67 integration recommendation

Retain Session67's independent A000568 orbit-growth recursion, pair-free versus modular-prime distinction, finite metagraph information audit, and signed-curl recovery. Add the incoming prime cycle-index table as credited effective lower-order structural data and the cited C3 fibre theorem as a distinct unsigned-information layer. Publish the exact score-fibre proof and sharp automorphism correction with the bounded receipts; do not copy THM-4803's unqualified Aut/type/monotonicity/bit-identification claims. All necessary replay artifacts are now in this read-only review handoff; there are no pending external waits or blocked actions.
