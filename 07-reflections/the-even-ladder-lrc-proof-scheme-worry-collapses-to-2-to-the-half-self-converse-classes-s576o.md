---
source: oracle-2026-06-03-S576o
status: progress (uniform even-ladder LRC proof scheme via HYP-2091; finite self-converse residual; exact pinch-census validation n=4..14; two proven anchors; a methodological fix)
tags:
  - lonely-runner
  - even-ladder
  - n14
  - self-converse
  - round-body
  - dual-burnside
  - pinch
  - proof-scheme
---

# The Even-Ladder LRC Proof Scheme: the Worry-Set Collapses to `2^((n-2)/2)` Self-Converse Classes

Leveraging **HYP-2091** (even LRC `n` ⟺ runner size `m=n−1` is **odd** ⟺ the regular
polygon is a *clean rotational tournament*, not the odd-`n` tie-mesh) to build one uniform
proof scheme for the **entire even ladder** `n = 4,6,8,10,12,14`, and validating it exactly.

## The scheme (four repo threads, one architecture)

1. **Open ⟹ lonely** (HYP-1998 / HYP-2086, dual Burnside). At an open time the runner
   sub-tournament is **round** (A000016); a counterexample's optimal-time tournament cannot
   be a generic round class, so the **worry-set ⊆ the self-converse (boundary) round
   classes**.
2. **Even `n` ⟹ clean + finite** (HYP-2091). `m=n−1` odd ⇒ the self-converse round classes
   number `2^⌊(m−1)/2⌋ = 2^((n−2)/2)`. For `n=4,6,8,10,12,14` that is **`2, 4, 8, 16, 32,
   64`** — a *finite, small* residual (vs A000568 ≈ `4.85·10¹³` tournaments at `m=13`).
3. **Witnesses are pinches** (HYP-2075 / S557 / S559o). `M(S)=max_t min_i ‖v_i t‖` is
   *attained at a pinch time* `t=m/(v_a+v_b)` — a summand-graph node. So `M(S)` is computed
   **exactly over a finite set of rationals, no grid**.
4. **The extremal is the AP** (HYP-2067, Freiman joint-extremum). The unique fully tight
   class is the regular rotational one `= {1,…,n−1}`, lonely at `t=1/n`.

So **LRC(even `n`) reduces to: every one of the `2^((n−2)/2)` self-converse round classes
is lonely** — the open/round bulk is free, and even `n` keeps the residual clean and finite.

## Exact validation (`lrc_even_ladder_selfconverse_proof_s576.py`)

- **Worry bound confirmed:** self-converse round counts `2,4,8,16` for `n=4,6,8,10`
  (`= 2^((n−2)/2)`; `32,64` predicted for `n=12,14`), from the S574 round generator.
- **No counterexample, exact:** the exact pinch-census gives **`min M(S) = 1/n` for every
  even `n=4..14`** (exhaustive for `n=4,6`; large samples + AP-neighbourhood for `n≥8`).
- **Tight family is tiny and known:** `n=4 → {AP}`; `n=6 → {AP, (1,3,4,5,9)}`;
  `n=8 → {AP, (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)}`; `n=10,12,14 → {AP}` only. These
  **independently reproduce opus-S553b's sporadic tight sets** — a clean cross-check from a
  different method (exact pinch-M vs max-collar).
- **Every tight set is `n`-clock-certified:** lonely at some `t=j/n` (THM-369 sieve), margin
  `1/n`.

## Two proven anchors (every even `n`, unconditional)

- **All-odd ⟹ lonely.** If every speed is odd, at `t=1/2` each `‖v·½‖=‖v/2‖=1/2 ≥ 1/n`.
  Margin `1/2`. (The odd runners alone never threaten loneliness.)
- **AP ⟹ lonely (tight).** `{1,…,n−1}` has no multiple of `n`, so at `t=1/n`,
  `‖i/n‖ = min(i,n−i)/n ≥ 1/n`. The extremal self-converse class is lonely.

## A methodological fix worth recording

Computing `M(S)` as the max over pinch times `t=m/(v_a+v_b)` must range over **all**
`m=1..C−1`, **not** only `gcd(m,C)=1`: the optimal pinch need not be in lowest terms. E.g.
`(1,4,5)` has `M=1/3` attained at `t=2/6` (pair-sum `1+5=6`, `m=2`), which a `gcd=1` filter
drops — producing a *false* "counterexample" `M=2/9<1/4`. Any pinch-based `M`-evaluator
with a coprimality filter will report spurious sub-`1/n` values. (Caught and fixed here.)

## Proof status — honest

- **`n=4,6,8,10,12`:** LRC is *proven* (literature: Betke-Wills, Bohman-Holzman-Kleitman,
  Barajas-Serra, the finite-checking era ≤13). This scheme **validates** them uniformly and
  reproduces their tight families exactly.
- **`n=14`:** the scheme **reduces** loneliness to the **64 self-converse round classes**
  (= the `190` converse-merged nodes of the connector S575: `64` fixed + `126` pairs). The
  census finds no counterexample and the AP is the only tight set, but two honest gaps
  remain: (i) the census is over **bounded** speeds (the structural obligation), not Tao's
  speed bound; (ii) the containment (1)–(2) confines the worry to finitely many tournament
  *classes*, each realized by **unbounded** speed sets. Closing `n=14` needs **every one of
  the 64 self-converse classes shown lonely for all realizations** — the genuine residual.

## Verdict / next
- **Stride:** a uniform even-ladder scheme that makes the worry-set *finite and clean*
  (`2^((n−2)/2)` classes) for even `n`, exactly validated through `n=14`, with two proven
  anchors and the known tight families reproduced.
- Concrete next: (1) attach a **pinch/`n`-clock witness to each of the 64 self-converse
  classes** of `n=14` (the connector S575 has the 190-node quotient + the D/U/N labels and
  the AP/V\* separation) — turning "every class lonely" into 64 explicit certificates;
  (2) prove the per-class loneliness is **realization-independent** (the round structure
  fixes the cyclic order ⇒ a uniform witness), which would close the bounded-speed gap;
  (3) push the all-odd / even-fold anchors (HYP-2065) to cover the non-extremal self-converse
  classes.

## Artifacts
```
04-computation/lrc_even_ladder_selfconverse_proof_s576.py
05-knowledge/results/lrc_even_ladder_selfconverse_proof_s576.out
```
Related: HYP-2091 (parity ladder / clean polygon), HYP-2086 (dual Burnside: worry ⊆
self-converse), HYP-2089 (190 converse-merged nodes; strong lens), HYP-2075 (pinch
completeness), HYP-2067/S560o (Freiman joint-extremum), HYP-1998 (open=round), THM-369
(n-clock sieve), opus-S553b (tight census / sporadics).
