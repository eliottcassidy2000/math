# HYP-2463 - LRC14 Q27 hard resources are non-stacking

**Status:** OPEN synthesis with exact finite support lemma.

**Source:** codex-2026-06-13.  Extends HYP-2459's parity-typed ledger, HYP-2444's Pisano quotient residual, HYP-2443's marked ladder support gate, and the HYP-2460/HYP-2461/HYP-2462 resonance-bonus lesson from unit distances.

## Claim

For the LRC(14) proof route, the one-stranger obstruction from HYP-2444 behaves like a resource packet, not like an indefinitely stackable bad direction.

Let

```text
CORE = 7*{1,...,12}
HARD = {260,351,442,611,702,793,962,1053}
Q27 = {d*m : d in {1,2,7,14}, 1 <= m <= 27} \ {1}.
```

For every `R subset HARD` with `|R|=k` and every `D subset CORE` with `|D|=k-1`, the 13-speed row

```text
S(D,R) = (CORE \ D) union R
```

has a strict LRC witness at some `q in Q27`.  Equivalently: even after all eight hard shell-27/13-clock residues are stacked into the 7-core frame, the fibered Q27 ladder always opens.

The proof target is stronger than this finite lemma:

> Any primitive LRC14 row that spends the shell-27 missing/zero class and the 13-clock resource must either descend toward the AP/Vstar/2AP wall atoms, open a low clock, open a divisor-fiber witness in Q27, or expose an odd marked owner/Bprime deletion.  These resources should be independent enough that they cannot all be spent indefinitely.

## Exact Finite Evidence

Script: `04-computation/lrc14_parity_typed_q27_ledger_codex.py`

Stored output: `05-knowledge/results/lrc14_parity_typed_q27_ledger_codex.out`

The script scans the complete hard replacement hull:

```text
sum_{k=1}^8 binom(8,k)*binom(12,k-1) = 77520 rows.
```

Results:

- `77520/77520` rows have a Q27 witness.
- `0` rows miss Q27.
- Only `10` rows miss the plain shell `q<=27`.
- The eight original one-stranger rows are the only plain-shell misses at `k=1`.
- The only non-original plain-shell misses are:
  - `delete (28), add (351,1053)`, caught by `q=30`.
  - `delete (28,63), add (351,962,1053)`, caught by `q=34`.
- The only rows whose first Q27 witness is as late as `q=91` are the two original one-stranger rows `r=611` and `r=702`.
- Once `k>=4`, every replacement row already has a plain `q<=27` witness.

This is not merely a larger brute-force check.  The strengthened script uses a bitset safe-twist mask for each `(q,speed)`, so the replacement hull can be audited as an exact finite support lemma rather than a fragile long search.

## Interpretation

HYP-2444 found the packet

```text
shell27 missing/zero class + 13-clock
```

as the one-stranger plain-shell obstruction.  HYP-2463 says that packet has poor tensor power: when multiple copies are forced into the 7-core frame, deleting core runners exposes other clocks faster than the shell-27 obstruction can accumulate.

The hard packet has two escapes:

1. **Divisor-fiber escape.**  The original `q=91=7*13` rows are not new shell-27 walls; they are fiber-address rows.  Their missing coordinate is the `7 x 13` address.
2. **Low-clock collapse.**  Replacing enough core speeds opens `q<=13` in most rows, and deleting `56`, `77`, or `84` opens a low Q27 witness in every checked replacement row containing that deletion.

Thus the obstruction appears stack-hostile, analogous to the unit-distance result that `27=3^3` is resonance-bonus hostile: a bad local packet can exist, but it cannot carry the crossing resource at the next multiplicative layer.

## Proof Route

The next proof target is a resource-independence lemma.

1. **Normalize the packet.**  Prove that any row blocking all plain `q<=27` shells has a marked subpacket projecting to the HYP-2444 types: shell-27 class `0` or `+-10`, together with a 13-clock constraint.
2. **Localize replacement cost.**  Show that adding each hard packet forces deletion, domination, or owner transfer in the 7-core quotient.  The finite hull suggests deletion of core speeds is not neutral: it opens low clocks or a fiber address.
3. **Promote from finite hull to arbitrary rows.**  Replace "row is exactly CORE with replacements" by a compression/descent statement: if a primitive row has no Q27 witness, its parity-typed ledger compresses to the hard replacement hull without losing blockedness.
4. **Use marked odd channels.**  HYP-2459 says invariant scalars live in even channels, while owner/carry/deletion data live in marked odd channels.  A would-be wall must keep the even scalar obstruction and suppress the odd opening simultaneously; HYP-2463 predicts that this cannot persist through Q27.

## Assumption Challenge / Tournament Analysis

This session did not use "runner" as the default tournament vertex.  Candidate vertices considered:

- runners,
- hard residues,
- deleted core speeds,
- denominators,
- divisor addresses,
- safe twists,
- owner/Bprime targets,
- proof obligations.

The selected Tournament Analysis uses proof obligations as vertices.  It preserves which quotient is needed for the proof and destroys raw time geometry.  The stored proof-obligation tournament is transitive, with leader `typed_Q27_ledger`, followed by `hard_resource_replacement_hull`, `deleted_core_low_clock`, `shell27_x_13_packet`, `Bprime_target_transport`, and `resonance_bonus_analogy`.

The challenged assumption is that a hard one-stranger obstruction should remain hard when copied.  The complete replacement hull says the opposite: the packet is fragile under stacking.

## Relation To Incoming Work

The unit-distance `n=27/28` work (HYP-2460, HYP-2461, HYP-2462) changed the reading of `27`: the number can be a bonus-hostile carrier even when nearby structure is rich.  Here, `C=2n-1=27` plays the same cautionary role.  The LRC shell-27 obstruction exists, but the hard packet is hostile to self-stacking once the missing divisor/fiber address is retained.

The analogy is not an equivalence.  Unit distance uses resonant products and Euclidean distance `sqrt(t)` packets; LRC uses torsion clocks, antipodal quotient classes, and marked owner/carry support.  The shared principle is that scalar product-count shadows become misleading unless the missing address coordinate is reattached.

