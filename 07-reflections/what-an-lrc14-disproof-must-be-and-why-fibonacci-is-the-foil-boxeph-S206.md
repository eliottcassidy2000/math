# Necessary anatomy of an LRC(14) disproof — and a Fibonacci foil heuristic

> **SCOPE CORRECTION — MISTAKE-221 (2026-07-21).** The rigorous necessary
> conditions recovered here are: after gcd normalization, a counterexample may be taken
> primitive; it must be Cover14 with `M<1/14`; its maximum-deletion core is
> non-AP by THM-1017, hence has a strict Schur-triple deficit by THM-730. The
> stronger “near-AP, anti-golden, higher-order-autocorrelation” description is a
> search heuristic, not a characterization. THM-731's `disc_v` is peel-specific
> and **smaller** discrepancy strengthens its sufficient safety certificate.
> The script's `q<=Qmax` scan returns a lower bound for `M`, not an exact maximum
> without a separate breakpoint-completeness proof. Its values above `1/14`
> safely exclude the displayed examples; they do not rank all Fibonacci or
> covering families. Read
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

## Finite candidate scan (denominator-truncated lower bounds)

A rational-time scan over the AP, near-AP perturbations, generalized/dilated
APs, and a few Fibonacci/Zeckendorf-structured covering sets finds for each
displayed family a witness above `1/14`. Except where canon separately proves
the exact value, the table entries are sampled lower bounds, not exact maxima:

| set | sampled lower bound | note |
|---|---|---|
| deep well `{1..12,182}` | `14/183 ≈ 0.07650` | exact AP-deletion benchmark, at `t*=14/183` |
| `2·AP = {2,4..24,182}` | `7/92 ≈ 0.07609` | **below `14/183` but > `1/14`** — see below |
| `AP12 + far=364/5460` | `≈ 0.0767` | sampled witness ≥ deep well |
| Fibonacci-12 `+` blocker | `5/29 ≈ 0.172` | displayed foil is loose |
| Fib-ish covering | `2/23 ≈ 0.087` | displayed packet is loose |
| random covering (×6) | `0.13 – 0.24` | finite sampled controls |

Two things sharpen the experiment. **(a)** The displayed Fibonacci-flavored
sets have large safety witnesses, supporting their use as hostile controls but
not a universal ranking. **(b)** Normalize primitively. The one sampled set
below `14/183` — the dilation
`2·AP = {2,…,24,182}` with lower bound `7/92` — is **non-primitive** (`gcd=2`); since covering-min is
dilation-invariant, `M(2·S)=M(S')` with `S'={1,…,12,91}` a set that omits a multiple of 14 and is
therefore lonely at `t=1/14`. It stays **above** `1/14`. This example explains
why covering is not dilation-invariant even though `M` is, and why normalization
must precede the Cover14 test.

So the search excludes its finite candidate list and motivates an anti-golden
comparison experiment. It does **not** prove that every counterexample is
near-AP or anti-golden, that Fibonacci families are uniformly loose, or that
Wall A is equivalent to joint-order autocorrelation extremality. Those are
questions exposed by the experiment.

Links: HYP-8815, HYP-7310 (Wall A), THM-730, THM-731/732, THM-1017, THM-724/726,
[[fibonacci-is-the-covering-mins-foil-not-its-lever-the-anti-golden-eisenstein-sibling-klein-S124]],
[[cf-descent-explains-the-far-element-lcm-13-14-and-stops-three-elementary-tools-converge-the-AP-core-is-open-boxeph-S103]].
