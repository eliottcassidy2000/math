# Necessary anatomy of an LRC(14) disproof — and a Fibonacci foil heuristic

> **SCOPE CORRECTION — MISTAKE-221 (2026-07-21).** The rigorous necessary
> conditions recovered here are: after gcd normalization, a counterexample may be taken
> primitive; it must be Cover14 with `M<1/14`; its maximum-deletion core is
> non-AP by THM-1017, hence has a strict Schur-triple deficit by THM-730. The
> stronger “near-AP, anti-golden, higher-order-autocorrelation” description is a
> search heuristic, not a characterization. THM-731's `disc_v` is peel-specific
> and **smaller** discrepancy strengthens its sufficient safety certificate.
> The corrected script uses the exact THM-1002 pair-sum engine on every row,
> after displaying a strict failure of the old `q<=Qmax` scan. All fifteen
> listed rows are exact and safe, but a finite bank gives no global ranking.
> Read
> [`CURRENT-FRONTIER.md`](../00-navigation/CURRENT-FRONTIER.md) first.

*boxeph-2026-07-21-S206. Owner: mine the repo for connections to 12 (esp. Fibonacci); think about a
disproof construction for LRC(14) and refine what it must be. Builds on klein-S124's Fibonacci-foil
proposal, boxeph-S103 (CF descent and `lcm(13,14)=182`), the stratified
THM-724/726 covering results, THM-730/731/732, THM-1017, and HYP-7310.*

## The `12`, the AP, and the anti-golden Eisenstein benchmark

The exact **deep-well benchmark** `{1,…,12,182}` has value `14/183`:
the `12` non-max speeds are the arithmetic progression `{1,…,12}`, and Cover14
forces the added speed to be divisible by both `13` and `14`, hence at least
`182=lcm(13,14)`. This proves the minimum in the AP-deletion stratum, not the
global LRC(14) covering minimum. Its maximizing time is
`t* = 14/183 = [0;13,14]`. The
arithmetic of `t*` is **Eisenstein**: `Φ₆(14)=14²−14+1=183`, `14` is a primitive 6th root of unity mod
183, the powers cycle with period 6 under `s_{k+2}=s_{k+1}−s_k` (the `x²−x+1` recurrence, `ℤ[ω]`, Heegner
`−3`). And `[0;13,14]` has **large** partial quotients, the fast-converging
pole of this comparison.

**Fibonacci is an exact sign-flip and an opposite continued-fraction pole:**
golden `x²−x−1`, `[0;1,1,1,…]`, and slow convergence contrast with this
Eisenstein/large-partial-quotient benchmark. That makes Fibonacci-structured
packets useful hostile controls. The contrast alone gives no uniform theorem
that all golden packets are safer than all non-golden packets.

## Proved necessary kernel

LRC(14) holds for a 13-speed set `S` iff
`M(S)=max_t min_{v∈S}‖vt‖≥1/14`. The AP deep-well benchmark is
`14/183≈0.07650`; a disproof must therefore lie below it by more than
`14/183-1/14=13/2562`. Precisely, a counterexample can be normalized and must
satisfy:

1. **Primitive after gcd normalization.** Dividing all speeds by their gcd
   preserves `M`.
2. **Cover14 and `M(S)<1/14`.** If no speed is divisible by some
   `q∈{2,…,14}`, then `t=1/q` is already a lonely witness.
3. **Non-AP maximum-deletion core.** THM-1017 excludes
   `S\{v_max}=d{1,…,12}`. THM-730 then gives its Schur-triple count
   `T(S\{v_max})≤65`, since the AP alone attains `C(12,2)=66`.
4. **Large relative peel discrepancy.** A counterexample has lonely-time
   measure `L=0`. THM-731 therefore forces, for every peel `v`,
   `disc_v≥6|G'_{~v}|²`; otherwise its rigorous lower certificate for `L`
   would be positive.

These facts do not determine whether the core is near an AP, the continued-
fraction type of a maximizing time, or which runner carries each modulus.

## The original joint-order idea, repaired as a question

> Can a non-AP twelve-set with `T≤65` nevertheless force
> `disc_v/|G'_{~v}|²≥6` for every legal peel, and can such a packet also satisfy
> the Cover14 and pointwise constraints?

This is a legitimate joint-order search question, not an equivalence with Wall
A. THM-730 controls the triple count; THM-731/732 control a peel-specific
autocorrelation discrepancy. No theorem currently resums one into the other.

## Heuristic search directions

- **Fibonacci/golden controls.** Test these as the opposite CF pole; do not
  assume they are uniformly safe.
- **Near-AP, anti-golden candidates.** A useful experimental family is a
  generalized AP or Bohr set that trades a little triple count for different
  higher-order autocorrelation. In the THM-1979/2013 continuum language the AP
  is the cold transitive point (`τ=0`); nearby hotter points are test cases, not
  a complete search universe.
- **Alternate modulus ownership.** Do not assume the maximum speed alone covers
  `13` and `14`; multi-owner and multi-defect packets remain live.

## Corrected finite audit

The companion script first preserves the one-sided bounded-scan logic, then
uses THM-1002-pair-sum-denominator-bound §1 to compute every row exactly by
enumerating every numerator on every pair-sum ruler `q=v_i+v_j`. It normalizes
gcds and explicitly validates that each input contains thirteen distinct
positive speeds. Non-coprime numerators are retained;
dropping them would miss maxima whose reduced denominator merely divides a
pair sum.

The row `{1,...,12,5460}` is a concrete regression test for the old bug:

`L_1200=92/1197 < M=420/5461`.

More generally, if `Q>=14`, `d=lcm(2,...,Q)`, and the deep well is dilated by
`d`, then every old sampled phase has value zero while dilation invariance gives
the true value `14/183`. Completeness and gcd normalization are separate issues.

Selected rows from the frozen output are:

| normalized family or control | result | logical status |
|---|---:|---|
| `{1,...,12,182}` | `14/183` | exact and safe |
| `{1,...,12,364}` | `28/365` | exact and safe |
| normalized `2*AP` row | `7/92` | exact and safe; misses Cover14 after normalization |
| `{1,...,12,5460}` | `420/5461` | exact and safe; old cutoff was strictly smaller |
| Fibonacci-12 plus `720720` | `5/29` | exact and safe; `720720=lcm(1,...,16)` |
| six original valid random completions | exact values `>=4/31` | all safe |

All fifteen valid tested rows have exact `M>=1/14`; no row is a disproof. This
is a rigorous finite safety result, not a global extremal theorem. The earlier
random generator lacked a distinctness assertion, although its six displayed
seed-2026 rows happen to be valid; the corrected generator preserves those
rows and validates the condition explicitly. In particular, its old fourth
row improves from the truncated `69/467` witness to the exact `56/379`.

## Connection to the parity–Hasse audit

THM-2043 supplies a sharper warning about coarse search statistics. The AP and
the infinite family

`{1,...,11,13,96+3444n}`

have identical owner-residue, raw phase, complete characteristic-seven Hasse,
`q<=13` blockedness, and threshold packets at period 14, yet the latter has the
strict resolved witness `(q,a,integer slack)=(41,17,1)`. Even adjoining any
fixed finite number of 7-adic lift-height digits does not separate them.

That makes the exact resolved-phase certificate

`C_{q,a}(S)=min_v(14 dist(av,q)-q)`

a more theorem-facing search coordinate than an unlabelled autocorrelation
scalar. Higher-order autocorrelation may still help *supply* a phase, but the
carrier must retain owner and height information through the exit.

## Honest verdict

No LRC(14) disproof was found. HYP-8815 is retained as a corrected heuristic
program: search primitive representatives of valid thirteen-speed sets using
AP distance, continued-fraction structure, autocorrelation profiles, and
resolved phase slack. Its rigorous contribution is an exact application of the
prior THM-1002 pair-sum lemma, giving fifteen finite safety certificates,
joined to the necessary THM-731 discrepancy obstruction above. LRC(14), AP
extraction, and the claimed joint-order autocorrelation extremality remain
open.

Links: HYP-8815, HYP-7310, THM-730, THM-731/732,
THM-1002-pair-sum-denominator-bound, THM-1017, THM-2043,
MISTAKE-221,
[[fibonacci-is-the-covering-mins-foil-not-its-lever-the-anti-golden-eisenstein-sibling-klein-S124]],
[[cf-descent-explains-the-far-element-lcm-13-14-and-stops-three-elementary-tools-converge-the-AP-core-is-open-boxeph-S103]].
