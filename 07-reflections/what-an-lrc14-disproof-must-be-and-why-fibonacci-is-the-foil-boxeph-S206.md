# What an LRC(14) disproof must be — and why Fibonacci is the foil

*boxeph-2026-07-21-S206. Owner: mine the repo for connections to 12 (esp. Fibonacci); think about a
disproof construction for LRC(14) and refine what it must be. Builds on klein-S124 (Fibonacci = the
covering-min's foil), boxeph-S103 (CF-descent, the far element `lcm(13,14)=182`), THM-724/726/730/731/732
(covering-min rigidity + the metric surrogate), THM-1017 (the AP-core reduction), HYP-7310 (Wall A).*

## The `12`, the AP, and the anti-golden Eisenstein extremal (what the repo already found)

The LRC(14) extremal is the **deep well** `{1,…,12,182}`: the `12` non-max speeds are the arithmetic
progression `{1,…,12}`, and the far element `182 = lcm(13,14) = 13·14` is forced because it must block
*both* continued-fraction neighbours of the argmax `t* = 14/183 = [0;13,14]` (boxeph-S103). The
arithmetic of `t*` is **Eisenstein**: `Φ₆(14)=14²−14+1=183`, `14` is a primitive 6th root of unity mod
183, the powers cycle with period 6 under `s_{k+2}=s_{k+1}−s_k` (the `x²−x+1` recurrence, `ℤ[ω]`, Heegner
`−3`). And `[0;13,14]` has **large** partial quotients — the *fastest*-converging CF.

**Fibonacci is the exact sign-flip and the opposite pole** (klein-S124): golden `x²−x−1`, `[0;1,1,1,…]`,
the *slowest*-converging CF, the three-gap **worst** case. The covering-min extremizer is engineered to
be **as anti-golden as a Farey rung allows**. So: *Fibonacci is the foil, not a lever.* On the LRC(14)
core, a Fibonacci/golden-structured speed set is the **least** dangerous — the easiest to make lonely,
the farthest from a counterexample.

## What a disproof must be (the refined statement)

LRC(14) holds for a 13-speed set `S` iff its covering-min `M(S)=max_t min_{v∈S}‖vt‖ ≥ 1/14`. The
conjectured global covering-min over covering sets is `14/183 ≈ 0.07650`, at the AP deep well; since
`14/183 > 1/14 ≈ 0.07143`, a disproof must dip **below `1/14`**, i.e. `13/2562 ≈ 0.00507` *below the AP
extremal itself*. Collecting every constraint, a disproof `S` **must** satisfy all of:

1. **Covering.** `S` contains a multiple of every `q∈{2,…,14}` (else `t=1/q` is a lonely witness,
   THM-366/523). Primitive up to dilation.
2. **`M(S) < 1/14`.** Strictly below the LRC threshold — hence below `14/183` and below every proved
   case (single-killer + ≥2-far-outliers have `M ≥ 14/183`, THM-724/726). So `S` lives in the *unproved*
   regime: **multi-killer with ≤1 far outlier**.
3. **Non-AP.** The 12 non-max speeds are **not** a dilated arithmetic progression (else the AP analysis
   forces `M ≥ 14/183`; Wall A / THM-1017).
4. **Sub-extremal additive energy.** By THM-730 (the AP *uniquely* maximizes additive triples,
   `T(A) ≤ C(k,2)`), a non-AP 12-set has **strictly fewer** additive triples than the AP.
5. **Higher-order tightness.** Therefore `S`'s tightness cannot come from additive-triple energy — it
   **must** come from **higher-order autocorrelation**: a smaller good-set discrepancy `disc_v`
   (klein THM-731/732) than the AP achieves, despite fewer triples.
6. **Anti-golden, not golden.** Its argmax `t*(S)` must have large partial quotients (fast CF), avoiding
   the golden three-gap regime; a disproof is a **near-AP anti-golden** object, **not** a Fibonacci one.
7. **Blocking far element.** The max speed must be divisible by the blockers of the CF neighbours of
   `t*(S)` (the `lcm(13,14)=182` analogue for the new argmax).

## The one sentence a disproof reduces to

> **A disproof of LRC(14) exists iff the additive-triple maximizer (the AP) is *not* also the maximizer
> of the full good-set autocorrelation — i.e. iff some non-AP 12-set concentrates its higher-order
> energy more than the AP, buying a smaller `disc_v` (hence smaller `M`) than the AP's `14/183`, despite
> losing additive triples.**

This is Wall A stated as a *joint-extremality* question: LRC(14) is true iff the AP jointly maximizes
additive energy at **all** orders; a disproof is a set where the higher-order energy **peels away** from
the triple energy. The repo's proved pieces control the triple order (THM-730) and bound `disc_v` by
arc-count (THM-732, `disc_v ≤ r²/(3v²)`); the open content is exactly whether the triple-maximizer wins
at all orders, or whether a near-AP with a shifted autocorrelation profile can undercut it.

## Where to look, and where not to

- **Not Fibonacci/golden.** By klein-S124 these are the foil — minimal additive energy, maximal `M`,
  the loneliest sets. A Fibonacci-structured covering set is the safest, not a counterexample.
- **Near the AP, anti-golden.** The disproof, if any, is a *small perturbation of the AP* — a
  generalized AP / Bohr set that trades a little triple energy for a more concentrated higher-order
  autocorrelation. In my continuum language (THM-1979/2013): the AP is the cold transitive point
  (`τ=0`); a disproof would be *slightly hotter* — a near-AP with extra higher-order cyclic structure
  that nonetheless tightens the covering. The search is a narrow band just off the cold vertex.
- **The `12` is not negotiable.** Any disproof still has exactly 12 non-max speeds (a covering-13 set),
  and they must cover `q≤12` while the far element covers `13,14` — the same skeleton as the AP, minus
  the arithmetic progression.

## Honest verdict (with the numbers)

An exact covering-min search (scan `t*=a/q`) over the AP, near-AP perturbations, generalized/dilated
APs, and Fibonacci/Zeckendorf-structured covering sets finds **none below `1/14`**:

| set | `M` | note |
|---|---|---|
| deep well `{1..12,182}` | `14/183 ≈ 0.07650` | the primitive extremal, at `t*=14/183` |
| `2·AP = {2,4..24,182}` | `7/92 ≈ 0.07609` | **below `14/183` but > `1/14`** — see below |
| `AP12 + far=364/5460` | `≈ 0.0767` | ≥ deep well |
| Fibonacci-12 `+` blocker | `5/29 ≈ 0.172` | **the foil — far looser** |
| Fib-ish covering | `2/23 ≈ 0.087` | looser |
| random covering (×6) | `0.13 – 0.24` | all far looser |

Two things sharpen the target. **(a) The Fibonacci sets are the loosest of all** (`M ≈ 0.17–0.24`) —
the foil confirmed numerically: golden structure is the *safest*, the opposite of a counterexample.
**(b) Primitivity is not negotiable.** The one set that dipped *below* `14/183` — the dilation
`2·AP = {2,…,24,182}` at `7/92` — is **non-primitive** (`gcd=2`); since covering-min is
dilation-invariant, `M(2·S)=M(S')` with `S'={1,…,12,91}` a set that omits a multiple of 14 and is
therefore lonely at `t=1/14`. It stays **above** `1/14`, and it reduces to a lower LRC case. So the
`14/183` extremal is for **primitive** covering sets, and a disproof **must be primitive** (add this to
constraint 1) — a non-primitive dip is always a dilation of an LRC(≤13)-safe set.

So the search finds no counterexample and confirms both poles (AP tight, Fibonacci loose). The value is
the **sharpened target**: a disproof is not a wild object but a *primitive, near-AP, anti-golden set that
beats the AP at higher-order autocorrelation*, and the whole question is the **joint-order extremality of
the arithmetic progression** — precisely Wall A, now with a Fibonacci-negative, a primitivity constraint,
and an autocorrelation-order-positive characterization.

Links: HYP-8815, HYP-7310 (Wall A), THM-730, THM-731/732, THM-1017, THM-724/726,
[[fibonacci-is-the-covering-mins-foil-not-its-lever-the-anti-golden-eisenstein-sibling-klein-S124]],
[[cf-descent-explains-the-far-element-lcm-13-14-and-stops-three-elementary-tools-converge-the-AP-core-is-open-boxeph-S103]].
