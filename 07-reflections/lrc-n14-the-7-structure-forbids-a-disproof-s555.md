---
source: opus-2026-06-01-S555 (remote-control)
status: honest NEGATIVE — a 7-steered disproof hunt finds no counterexample; 7 is a PROOF lever, not a disproof one
tags: [LRC, n14, n7, disproof-hunt, counterexample, even-fold, negative-result, S554, S555]
---

# LRC@14: the 7-structure forbids a disproof (it is a proof lever, not a disproof one)

**Prompt (user):** consider the 7-impossibility the key to the *disproof* of LRC n=14.

**Honest framing.** LRC is believed true and no counterexample is known at ANY n;
a disproof of LRC(14) would be historic. I therefore treated this as a rigorous
hunt for a config with max-collar `M(S)=max_t min_i ||v_i t|| < 1/14`, verified
every candidate with exact arithmetic, and report the result truthfully. **No
counterexample was found**, and the 7-structure is exactly why one is not expected.

## Where 7 says a counterexample could hide

The even-fold bridge (S554): even `v=2u` gives `||vt||=||u·2t||`, so with
`fold=halved evens`, `M14(S) ≤ M(fold(S))`. For `e=|fold| ≤ 6`, **LRC(7) forces
`M(fold) ≥ 1/7`** and the preimage construction produced witnesses 127/127. Hence
any counterexample must live in the **unprotected `e ≥ 7` regime**, where
`M(fold)` is an unproven LRC(8+) quantity. So the *sharpest* disproof probe is to
make the even part a **doubled tight AP** (`M(fold)=1/(e+1)`, smallest the
unproven cases allow) and stack the odd coupling.

## The hunt (`lrc_n14_disproof_via_7_s555.py`, exact)

- **A. AP with any subset of its 6 even runners doubled** (63 configs): min
  `M=1/14`, attained *only* by `V*={1..11,13,24}` (the S553 sporadic). No `M<1/14`.
- **B. doubled-AP_k even part `{2,…,2k}` + (13−k) odd speeds, k=7..12**
  (`M(fold)=1/(k+1)`): every config is **loose** (`M>1/14`); none even tight.
- **C. mod-7 / 7-resonant random** (4000): all loose; min `M>1/14`.
- **D. measure-minimising hill-climb forced `e≥7`, seeded by doubled-AP_12 and a
  7-split** (speeds ≤120): **global min `M=1/14`, attained at the AP `{1,…,13}`** —
  the search, even forced into `e≥7`, migrates back to the AP as the minimiser.

> **Verdict: 0 verified counterexamples.** The minimum collar over the entire
> 7-steered / e≥7 / doubled-AP / hill-climb search is exactly `1/14`, achieved by
> the tight family (AP and V*).

## Why this is the expected answer

1. **The fold bound cannot produce `M<1/14`.** `M14(S) ≤ M(fold(S))` and
   `|fold| ≤ 13`, so `M(fold) ≥ 1/(|fold|+1) ≥ 1/14` whenever LRC holds for
   `≤13` speeds — and the only *unproven* slack (`|fold| ≥ 7`) was searched and is
   empty. The bound is an *upper* bound, so it never forces non-loneliness; it
   only caps how lonely the even part can be.
2. **The protected half is genuinely protected.** For `e ≤ 6` (≈51% of configs)
   LRC(7) + the antipodal-preimage construction yields a witness (S554, 127/127),
   so no counterexample can have `e ≤ 6`.
3. The `e ≥ 7` half, where the 7-fold gives no guarantee, is exactly where a
   disproof would have to live — and the targeted exact search there is empty.

So the "7 impossibility" is the key to **believing / proving** LRC@14, not to
disproving it. It localises the only conceivable counterexample to `e ≥ 7` and
then finds nothing there; combined with the S553 census (no `M<1/14` among
transversals or non-unit-pair non-transversals) and S551 (no random failures),
the minimum collar at n=14 is `1/14` across everything tested.

## Honest limitation

This is not a proof that LRC(14) holds: the searches are over bounded speeds
(≤120) and heuristic hill-climbs, so a counterexample with very large speeds is
not excluded by computation (only by the conjecture's truth and Tao's
bounded-minimal-counterexample reduction, whose bound is astronomically large).
But there is **no evidence whatsoever** for a disproof, and the structural
argument points the other way.

## Convergence with oracle-S552o (independent, same prompt)

A concurrent session (oracle-2026-06-01-S552o, HYP-2057) attacked the identical
"7 impossibility, for a disproof" prompt via the **complementary prime factor**:
the **mod-7 CRT decomposition** of the 13 runners into six pair-classes
`{i, i+7}` and the **singleton** `{multiples of 7}`. It reached the **same honest
verdict**: 7 is the right structural lens but yields **no impossibility** — its
natural 7-gon-window construction *provably fails* on the singleton class (the
window half-width `1/(14V)` is too small to clear even one generic multiple of 7).

So both factorisations of `14 = 2·7` tell the same story:

| lens | clean part | residual coupling |
|------|-----------|-------------------|
| **mod 2** (this work, S554/S555) | even-fold settled by LRC(7) (e≤6) | the odd antipodal split |
| **mod 7** (oracle S552o) | six pair-classes near-independent | the singleton {7·w} ↔ pairs |

Neither lens produces a counterexample; both *localise* the residual and both
point to the proof side. My S555 adds the direct exact-search confirmation (min
collar `= 1/14`, attained only by the tight family) to oracle's constructive/
measure analysis.

**Conclusion for the owner:** a 7-based disproof of LRC(14) does not appear
possible — by **either** factor of 14, the 7-structure is a *lower-bound* (proof)
mechanism, not a disproof one. Independently confirmed two ways. The productive
use of "7 is the key" remains the proof direction: close the `e≤6` no-odd-split
(S554) and bound the 7-way correlation on the singleton class (S552o).

**Artifacts:** `04-computation/lrc_n14_disproof_via_7_s555.py` (+`.out`).
Builds on S554 (even-fold lever), S553 (tight census, V*), S551 (sieve).
