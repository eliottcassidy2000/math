# Session67 — A000568, metagraph families and information in the relation

**Audited synthesis, 2026-10-26.** A000568 counts ordinary unlabelled
tournaments; converse is not identified unless expressly stated below.

The repository already contains complete **all-order structural and arithmetic
recurrences** for these classes. This session recovers and credits them, adds a
count-free local recognition certificate, independently enumerates the complete
order-eight metagraph carrier, and reconnects the construction to the older
relation-first information work. A substantially smaller closed family grammar
or scalar recurrence is still **OPEN**, not supplied by these exact algorithms.

## 1. Complete families were already available — and what is added here

**Inherited PROVED classification.** Every nonempty tournament has a unique
partition into maximal transitive modules. Contract them to a **pair-free core**
`Q`; retain their positive sizes `l`. Then

```
T = Q[TT_(l_1), ..., TT_(l_q)],
classes over Q = {positive size vectors}/Aut(Q),
Aut(T) = Stab_Aut(Q)(l).
```

Pair-free means no homogeneous pair, not modular indecomposable or strongly
connected. Homogeneous-pair adjacency itself is not transitive: its connected
components are the canonical chains. Arbitrary unreduced weighted quotients
can duplicate classes; the canonical quotient prevents that loss.

The October 20 `weighted_core_transport`, `defect_boundary`, and
`tournament_future_boundary` packets already provide alternating augmentation
of the complete core bank, signed cycle-index inversion, and full odd-cycle
operator recurrences. Their earlier ordinary-count evaluation through 200 and
core-count evaluation through 100 are **inherited**, not rerun census claims.
[Recovered sources, exact equations and scope](session67-metagraph/REPORT.md).

For each core the complete ordinary family generating function is

```
F_Q(x) = (1/|Aut(Q)|) sum_(g in Aut(Q))
             product_(cycles c of g) x^|c|/(1-x^|c|),
A(x) = 1 + sum_(nonempty pair-free Q) F_Q(x).
```

This is local Burnside on the full automorphism **action**, not merely its order.
The inherited structural recurrence closes the core bank as a set recurrence;
it is not a short recurrence on the scalar `a_n=A000568(n)`.

**New checked refinement.** `A[n,q]` counts classes jointly by total order and
canonical core order. Its full integer triangle is evaluated through **40**;
actual duplicate-free core families and automorphism actions are independently
verified through **7**. In particular, `A[8,8]=4602` out of `a_8=6880`.
Order-40 formula evaluation is not generation of every order-40 tournament.

The all-order core-deficit bound is

```
sum_(q<=n-d) A[n,q] <= binom(n-1,d) a_(n-d),  1<=d<n.
```

Together with an explicit ordinary-class Burnside estimate it proves
`A[n,n]/a_n -> 1`. More generally, as `q -> infinity`, **uniformly for every
n>=q**, `A[n,q] ~ a_q binom(n-1,q-1)`. This is not a fixed-q asymptotic.
A fixed finite core bank therefore misses almost all classes; this does not
rule out an unbounded symbolic grammar with a small description.
[Independent proof audit and hostile limits](session67-cycle-audit/CORE_TRIANGLE_AUDIT.md).
No priority claim is made for standard species identities or asymmetry results.

## 2. A count-free local enumerator, including the literal metagraph lines

**PROVED corollary of inherited all-order relocation connectivity.** Start with
all `2^binom(n-1,2)` binary relations having one fixed forward Hamiltonian path.
Union legal single-vertex path relocations. The components are exactly ordinary
tournament classes, without a graph-isomorphism call or a supplied A000568 value.
A proposed partition is correct iff each block is connected by legal moves and
**every** legal edge stays in its block. Neither condition alone suffices.

Add reflection to merge converse classes. Pair each cube word with its bitwise
complement to recover the literal coloured lines, retaining loops and repeated
lines. Relocation edges are auxiliary witnesses, not those complement lines.

**FINITE-EXACT, complete through order 8:** ordinary counts at n=0,...,8 are
`1,1,1,2,4,12,56,456,6880`. Order zero is the separate empty-class convention.
At order 8 there are **3,528 converse-merged nodes**, **176 self-converse
ordinary classes**, and **1,048,576 literal lines**, stored as 438,374 weighted
coloured-endpoint rows. Every class representative and normalized multiplicity is exported; all
lines are retained in the weighted endpoint ledger. A loop contributes twice to degree; ordinary multiplicity is
`tc(T)=H(T)/|Aut(T)|`, not raw Hamiltonian-path count `H(T)`.

Over any field, the relocation incidence quotient has dimension `a_n`.
The invariant-function algebra is a separate, equivalent class-coordinate
object, not a pointwise-algebra quotient by the edge span. A mod-two Laplacian
cannot replace incidence: the smallest hostile is an order-four `K_(2,3)`
component, with incidence rank 4 but Laplacian rank 2.
[Proof and algorithm](session67-root/COUNT_FREE_ENUMERATION.md) ·
[Reproduction](session67-root/REPRODUCE.md) · [Independent audit](session67-local-audit/AUDIT.md).
The implementation is bounded to n<=8 and uses exponential space.

## 3. Why present family symmetry does not determine future counting

The exact rooted-incidence contribution of parent `T` to child family `G` is

```
K_n(T,G) = sum_(attachments z yielding G) |Aut(T+z)|
              / ((n+1)|Aut(T)|).
```

This is a rational counting mass, **not a transition probability**. Correct
canonical augmentation selects a complete **root orbit**, not a minimum unrooted
deletion card; the independent implementation accepts each child exactly once
through order eight, before any global deduplication.

For every `k>=2`, the rigid strong parents
`A_k=C3[3k-3,1,1]` and `B_k=C3[k-1,k,k]` have the same actual core, order,
SCC-size sequence and full current automorphism cycle index. Nevertheless:

- For the **single child** `U_k=C3[k,k,k]`, their masses are **0 and 1**.
- For **all children with core order three**, their masses differ by **2/3**.

The missing coordinate is the weight vector's alignment with **future** symmetry.
This refutes universal transfer from the stated coarse family key, valid for
arbitrary assigned parent masses. It does not forbid scalar identities for the
actual global sequence or other family recurrences. The complete weighted core
remains sufficient. An SCC-interval law also gives a smaller exact observer for
strong **rooted** extensions, not strong unrooted counts.
[Recurrence and controls](session67-recursion/REPORT.md) ·
[Independent cross-audit](session67-cycle-audit/RECURRENCE_CROSS_AUDIT.md).

## 4. Recovering “tournaments as information” without losing the relation

The March-labelled information scripts and June relation-first hypotheses are
recovered with exact source pins, original snapshots and date qualifications.
Their useful viewpoint survives, but several literal identifications do not.
[Historical lineage, repairs and proofs](session67-information/INFORMATION_INTRANSITIVITY_RECOVERY.md).

For signs `S_ij` in `{+1,-1}`, put `c_ijk=S_ij+S_jk+S_ki`. The binary cyclic-triple
indicator is exactly `(c_ijk^2-1)/8`. Transitivity means a **sign ranking**
`S_ij=sign(p_i-p_j)`, not an additive gradient: a transitive triple already has
curl ±1. Scores plus the full signed, addressed curl recover the relation;
triangle count, curl energy, switching shadows and unsigned curls generally do not.

A second exact relation-first enumeration groups classes by their cyclic-triple
hypergraph `H`, then counts its **full realization orbits** `R(H)/Aut(H)`.
Positive-triple XOR constraints need negative-triple NAE constraints and the
signed group action. The hypergraph alone is not a complete class invariant.
The classical C3/interval-reversal theorem is **CITED**, not claimed as new.

Information requires a specified source. A uniform class source has entropy
`log2(a_n)`; a uniform labelled source weights class `T` by
`n!/(2^binom(n,2)|Aut(T)|)`. The canonical weighted-core datum is lossless for
**every** source, giving a conditional entropy chain rule, not independent
“order × shape” factors. The old path-entropy surrogate is not noisy-ranking
mutual information: at 50% comparison noise the latter is exactly zero.
Current-scope corrections retain the historical bodies and original snapshots.

## 5. Incoming prime families: credited arithmetic and a sharp correction

Klein's concurrent prime-cycle-index census supplies 23 aggregated rows through
order 9 and order-10 prime totals. A new independent Fraction implementation
reproduces the supplied-data ordinary counts through **10**, ending
`191536, 9733056`, and the C3-structure counts ending `93612, 4774417`.
This **does not regenerate the prime census**, and credits its producer.
Order-11 prediction still needs P11 (and self-converse P11 for C3), not merely
order-10 cycle indices. [Incoming review and pinned source](session67-incoming/INCOMING_INTEGRATION.md).

That review corrects an all-order claim without retracting the valid cited fibre
theorem or species work. Let `B=singleton -> C3`. On the same 12 labelled
vertices, `C3[B,B,B]` and `C3[B*,B,B]` have the **same 67 cyclic triples** and
**H=462915**, but automorphism orders **81 versus 27**. Reversing a descendant
changes its parent's decoration stabilizer. A proof using the cited shared
modular tree and odd automorphism orbits shows that order 12 is the **first
possible failure**; this is not an exhaustive census through order 11.

The same cited structural input yields a new exact refinement: for the full
**vertex-addressed score function** together with C3, the labelled fibre has
size **`2^b`**, where b counts size-balanced prime nodes. For every shared module M,

```
D(M) = (sum_(v in M) d(v) - binom(|M|,2))/|M|.
```

Child-minus-parent D isolates each weighted quotient score, so changes at
separate tree levels cannot cancel. This closes exactly HYP-9413(b), not its
other open clauses. Independent checks exhaust all 1,099 labelled tournaments
through 5 plus four explicit larger tree fibres; the order-12 automorphism
witness is checked by two independently written searches. The dualization bit
is not generally an alternating-group parity bit, and identity-observer fibre
loss is zero. Current source corrections retain pinned pre-correction originals.

## Evidence and remaining boundary

All-order assertions above have explicit proofs or identified inherited inputs;
finite checks are not promoted into all-order proofs. Local independent class
oracles reach 6, complete line/table agreement reaches 7, and the production
LUT/arc audit and metagraph mass accounting reach 8. The core-triangle and
recurrence cross-audits add 10,664 and 3,195 checks respectively. Information
and realization certificates exhaust their stated universes through order 5.

[Publication/replay gates](PUBLICATION_CHECKS.md) and [source provenance](PROVENANCE.json)
record the frozen packet. `python3 verify_packet.py` (also with `-O`) verifies
its exact file set and hashes; see the lane reports for reproduction commands.

**OPEN:** a substantially smaller closed all-order family grammar or scalar
recurrence, its closure under the requested line/extension operations, a proved
complexity advantage, and black-core connectivity. No LRC(14) result follows.
