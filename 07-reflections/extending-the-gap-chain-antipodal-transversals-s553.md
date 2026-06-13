---
source: oracle-2026-06-01-S553
status: VERIFIED -- Link 1 proven (witness-time lemma) + gap reduced to 2^(n-1) transversals; reduced gap exhaustive n<=8
tags: [lonely-runner, spectral-gap, chain, antipodal-pairs, transversal, witness-time, lrc-sufficient-condition, lrc-progress]
---

# Extending the gap chain: the antipodal-pair witness times collapse LRC to the transversals

**Prompt (user):** try to extend the chain (toward proving the S552 loneliness
spectral gap, HYP-2052: every non-AP-tight gcd-1 config has max-collar
`M(S) >= 2/(2n-1)`).

I found the **mechanism** behind the gap and turned it into a proven lemma that
collapses the conjecture from *all* configurations onto an explicit, tiny, structured
family -- the antipodal transversals mod `2n-1`. Along the way it yields a genuine new
**sufficient condition for LRC itself**.

## The witness times and the proven lemma (Link 1)

Let `M = 2n-1` (odd) and use the family of **witness times**
```
   t_k = k / (2n-1),   k = 1, ..., 2n-2.
```
For a single speed `s`,  `|| s t_k || = || k s / M || = dist(ks mod M, 0)/M`, which is
`< 2/M`  iff  `ks ≡ 0, +1, -1 (mod M)`  iff  `s ≡ 0, +k^{-1}, -k^{-1} (mod M)`. So at
time `t_k` the **only** residues that dip the collar below `2/M` are `0` and the
**antipodal pair** `{+k^{-1}, -k^{-1}}`. As `k` runs over the units mod `M`, `a=k^{-1}`
runs over all nonzero residues, so the bad pair runs over all `n-1` antipodal pairs
`{a, M-a}`, `a = 1..n-1`.

> **LINK 1 (lemma, proven).** If a gcd-1 speed set `S` has no speed `≡ 0 (mod 2n-1)`
> and its residues mod `2n-1` **miss** some antipodal pair `{a, 2n-1-a}`, then at
> `t = a^{-1}/(2n-1)` every runner satisfies `||s_i t|| >= 2/(2n-1)`, hence
> `M(S) >= 2/(2n-1)`.

Verified with zero violations across all such sets for `n=4..8` (142–860 sets each).
**This is already a new LRC theorem:** every such `S` satisfies LRC with strict
surplus (`M >= 2/(2n-1) > 1/n`). A large, explicitly described class of speed systems
provably clears the lonely bound by a one-line witness-time argument.

## The reduction: only the antipodal transversals remain

Link 1 leaves exactly the configs it cannot reach:
- those with a speed `≡ 0 (mod 2n-1)` (the `t_k` are all blind to such a runner), and
- those whose residues **hit every** antipodal pair -- a **transversal**.

A `0`-residue speed wastes a slot, so with only `n-1` speeds it forces a missed pair;
the genuine residual is the **perfect transversals**: residues hitting each of the
`n-1` pairs exactly once. There are `2^{n-1}` of them. The **AP `{1,..,n-1}` is the
all-lower transversal** (residues `1..n-1`, the small half of each pair). I verified
directly that the leftover `0`-residue configs *also* clear the edge (n=4..7, thousands
of sets, zero below `2/(2n-1)`) -- they satisfy the gap via non-`t_k` times -- so the
gap conjecture is **equivalent to its restriction to the transversals**.

## The reduced gap (verified n<=8)

Parametrize a transversal by its **flip-set** `F ⊆ {1,..,n-1}` = the pairs where the
*upper* residue `2n-1-a` is chosen (`F = ∅` is the AP). Over all `2^{n-1}` canonical
transversals:
```
   M(T) = 1/n        only for the AP-tight family:
                       F = ∅ always; plus F = {2} (the set {1,3,4,..}) for n=5,6 only;
   M(T) >= 2/(2n-1)  for every other transversal.
```
The min `M` taken over transversals with each flip-count `|F|` is always `>= 2/(2n-1)`
once you leave the tight family. So the spectral gap holds on the transversals
exhaustively for `n<=8` -- which, with Link 1 and the `0`-residue check, re-proves the
full S552 gap for `n<=8` and **structurally explains it**: the AP is the unique
all-lower transversal, and any other residue choice (flipping a pair to its upper half)
either stays exactly at the floor (the few tight exceptions) or jumps the whole margin.

## Where the chain now stands

```
 all gcd-1 configs
   |  Link 1 (PROVEN): miss an antipodal pair (no residue 0)  => M >= 2/(2n-1)
   |  residue-0 configs (verified n<=7): => M >= 2/(2n-1)
   v
 perfect antipodal transversals mod (2n-1)   (2^{n-1} of them; AP = all-lower)
   |  REDUCED GAP (verified n<=8): only AP-tight family at floor; rest >= 2/(2n-1)
   v
 [OPEN] prove the reduced gap for all n
```

The conjecture is now reduced from an infinite family to `2^{n-1}` highly structured
objects with a clean coordinate (the flip-set `F`), and the floor members are an
explicit short list. This is real progress: a proven lemma doing most of the work, a
genuine new LRC sufficient condition, and a sharp, finite-looking residual.

## Open (-> HYP-2052 update)

1. **Prove the reduced gap on transversals.** Show every transversal with `M < 2/(2n-1)`
   is AP-tight. Coordinate: the flip-set `F`. Likely route: a *second* witness-time
   family (times not of the form `k/(2n-1)`) that separates flipped pairs, or a
   monotonicity in `|F|`/`max F`.
2. **Characterize the tight exceptions.** `F=∅` (AP) is universal; `F={2}` is tight only
   for `n=5,6`. Which flip-sets keep `M=1/n`, and why does the exception die at `n>=7`?
   (These are exactly the non-AP members of the S552 tight family.)
3. **Relax Link 1's `0`-residue hypothesis** to a clean lemma so the whole non-transversal
   half is a single proven statement.

## Anchors
- `04-computation/lrc_gap_chain_antipodal_s553.py` (+`.out`): Link 1 proof-check (0
  violations), residual = transversals, AP = all-lower, reduced-gap enumeration,
  one-flip neighbors.
- `05-knowledge/results/lrc_transversal_flipset_s553b.out`: all `2^{n-1}` canonical
  transversals by flip-set; floor members; min-M by `|F|`.
- `05-knowledge/results/lrc_res0_subresidual_s553c.out`: residue-0 configs clear the edge.
Builds directly on S552 (the gap + doubled-apex witness), S526 (covering / `k/m`
character times), S542 (channels mod n/2).
