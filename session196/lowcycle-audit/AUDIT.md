# Session 196 two-apex automaton: independent audit

## Verdict

**VERIFIED, with the stated class restrictions.** The frozen packet gives a valid all-order cycle-language theorem and a valid finite certificate of the bounded scalar image at arbitrary order. No unresolved mathematical or implementation defect was found. A certificate-loading gap found during review was repaired before this final snapshot: `verify_language.py` now reads and checks the stored `cycle_dfa.json`.

The audited proof SHA-256 is `702f9ec3a6f784f71ae1fb6b3cbf02ed29c71fd1f953adf2fe2a386531121bd6`. The complete eight-file manifest is `source_hashes.json`; `source/` preserves the exact audited bytes. No repository checkout was modified by this audit.

## Independent checks

`independent_audit.py` reconstructs tournaments with the roots numbered first, counts paths/cycles by a last-vertex subset recurrence, independently implements sparse slot transfers, independently refines the output partition, and explores the bounded products. It does not import the proposed algorithms for these computations; its import of `apex_transfer.py` is used only as a third comparison target in the graph checks.

- All **5,461** words of base lengths 0 through 6 agree between independent adjacency counts, independent sparse transfers, and the proposed transfer.
- An additional **100** deterministic random words of base lengths 7 through 10 agree by all three routes.
- Independent cycle automaton: **251 reachable saturated vectors**, **32 minimal Moore states**. The independently rebuilt machine is equivalent to the stored JSON over its entire reachable product.
- Independent bounded product at `B=67`, cap 68: **3,388 states**, **3,999 retained transitions**, discovery depth **15**, missing strong odd values **{7,21,67}**.
- Independent bounded product at `B=501`, cap 502: **50,749 states**, **58,358 retained transitions**, discovery depth **20**, missing strong odd values **{7,21,67}**.
- All stored positive witnesses for both bounds have the claimed exact path count and cycle count 1 or 2.
- **1,200** additional length-40 word checks verify coordinatewise saturation against exact integers for all three transfers and caps 1, 2, 3, and 68.

`language_audit.py` uses Brzozowski derivatives rather than the packet's Thompson NFA construction. The complete product with the independent cycle machine has **32 reachable states**, all agreeing on whether `C=2`. This is a second all-order regular-expression equivalence proof, not a finite-length sample.

The independent machine inventory agrees with the proposed formula through length 99. More importantly, the residual recurrence

`A_n - 5 A_(n-1) + 8 A_(n-2) - 4 A_(n-3) = 0`

vanishes for 93 consecutive terms starting at `n=7`. Since this residual is a scalar sequence from the same 32-dimensional transition matrix, 32 consecutive zeros imply all subsequent zeros by Cayley–Hamilton. Together with the checked initial coefficients, this independently certifies the displayed rational generating function and all-order closed formula, without depending on the prose unambiguity argument.

The producer's transfer test, stored-certificate verifier, and both bounded-image runs were also rerun against the final snapshot. All eight source hashes still matched afterward.

## Proof and scope review

1. **Cycles and rotation.** Starting each directed Hamiltonian cycle at the unique distinguished root `a` chooses exactly one representative modulo cyclic rotation. The partition into the `a`-to-`b` and `b`-to-`a` base slots is unique. Reversal is not an additional identification. The empty base is a two-vertex tournament with `H=1,C=0`; the one-letter word `2` gives `H=3,C=1`.
2. **Nine cycle coordinates.** The first inserted vertex must satisfy its root-to-base boundary; subsequent inserted vertices automatically follow the existing last vertex because the base is transitive. The last-base-to-root condition is stored exactly by G/B. Empty-slot permission is asymmetric because `a->b` but not `b->a`. These facts prove the transfer bijection.
3. **54 path coordinates.** A path has exactly one root order and three increasing slots. The two path-order transfers are disjoint and complete. A terminal slot cannot enter B, so some of the declared 54 coordinates are unreachable; this is harmless redundancy, not an omitted case.
4. **Saturation.** For nonnegative integral transfer matrices, clipping at `K` is a homomorphism into the capped semiring: `min(K, sum x_i) = min(K, sum min(K,x_i))`. This is the precise invariant behind the cycle cap 3 and path cap `B+1`. No subtraction or cancellation is used.
5. **Why `B+1`, and why prefix pruning is sound.** Counts at most `B` are recovered exactly from capped accepting-coordinate sums. Larger counts remain larger than `B`. The proof's insertion map is injective because removing the newly inserted vertex recovers the original Hamiltonian word. Thus an over-bound induced prefix can never return below the bound. Cycle coaccessibility is an independent exact rejection condition. Every surviving word prefix is retained, and a finite closed reachability set therefore proves absence at every length. Discovery depth is not an order cutoff.
6. **Minimality and expression.** Output-preserving partition refinement and pairwise distinguishing continuations justify the 32-state Moore minimality claim. This is minimal for the four outputs `0,1,2,>=3`, not a claim about the minimal Boolean acceptor for `C=2`. The expression's first-exit branches are unambiguous as stated; its alphabet and concatenations match the implemented expression.
7. **Scope.** Orientation words include distinguished roots and a fixed transitive base order. Their inventory is not an unlabeled tournament enumeration and contributes no new A000568 coefficient. The bounded image concerns **two-apex tournaments with `C` equal to 1 or 2**, not all strong two-apex tournaments, not all low-cycle tournaments, and not the global `{7,21}` conjecture. The finite bound is on the scalar `H`, not on tournament order.

## Hostile controls

- Word `210` has `H=15,C=3`. If the cycle cap were incorrectly changed from 3 to 2, it would be a false `C=2` positive. The implemented cap 3 rejects it correctly.
- Word `12301` has `H=67,C=4`, confirmed independently from adjacency and slots. It is a concrete counterexample to extending the 67 exclusion to all two-apex tournaments. The implemented cycle filter rejects it correctly.
- The empty word and word `2` test the two-vertex and directed-triangle boundaries.

## Reproduction

Run from this directory:

```
python3 independent_audit.py
python3 language_audit.py
python3 source/apex_transfer.py
python3 source/verify_language.py
python3 source/bounded_image.py 67
python3 source/bounded_image.py 501
```

The machine-readable independent results are `audit_results.json` and `language_results.json`. The six logs record the final-snapshot runs. The audit script's elapsed-time field is diagnostic and is not expected to reproduce byte-for-byte.
