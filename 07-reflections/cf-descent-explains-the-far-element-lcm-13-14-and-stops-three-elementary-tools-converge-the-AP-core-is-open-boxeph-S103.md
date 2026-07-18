# Continued-fraction descent explains the far element `lcm(13,14)=182` and stops there: three elementary tools converge on one layer; the AP core is the open inverse theorem

*boxeph-2026-07-18-S103. Owner: prove the medium-modulus inverse theorem via continued-fraction descent.
Outcome: CF descent delivers a clean, genuine result — the far element is `lcm(13,14)=182` because it must
block **both** the CF convergent `1/13` and the adjacent mediant `1/14` — but it reaches only the far
element's **divisibility**, not the AP core. This is now the **third** elementary tool (after maximality
S101 and the sieve S102) to reach exactly the same 13/14-divisibility layer and provably stop. The AP core
is irreducibly additive — the open inverse theorem (= LRC(14), S94). Verified S103 computation.*

## What CF descent gives (a clean, genuine result)

For the primitive deep well `V={1,…,12,182}`, the maximizer is `t* = 14/183 = [0;13,14]` (verified). Its
Farey/CF structure near `t*`:

- **CF convergent `1/13`** (penultimate): for `t*` to be the maximum, `1/13` must not beat it, i.e.
  `min_v‖v/13‖ ≤ M < 1/13 ⟹ min_v‖v‖_{13}=0 ⟹ 13∣` some speed. Verified: killed by `182` (`13∣182`).
- **First mediant `1/14`** (`j=1` intermediate fraction `j/(13j+1)` between `1/13` and `14/183`): for
  `M<1/14`, `min_v‖v/14‖<1/14 ⟹ 14∣` some speed. Verified: killed by `182` (`14∣182`).

> **CF-descent result.** The far element is forced to be a multiple of **`lcm(13,14)=182`** precisely
> because it must **simultaneously block the CF convergent `1/13` and the adjacent mediant `1/14`** of the
> maximizer. `182 = 13·14` is the CF-descent signature of the deep well.

This is the cleanest explanation the project has of *why* `182` — it is the unique small blocker of the two
best rational approximations flanking `t*`. (For dilations the CF is scaled/messier, e.g.
`t*(3V)=14/549=[0;39,4,1,2]`; the statement is about primitive families, dilation-invariant.)

## Where CF descent stops (the wall, made precise)

CF descent, via the convergents/mediants of `t*`, reaches only the **far element's divisibility**
(`182∣v_max`). It does **not** deliver the AP core `V∖{v_max}=v_+·{1,…,12}`:

- The convergents of `t*=[0;13,14]` are just `1/13` and `14/183` — a single non-trivial small denominator
  (`13`). It yields one exact condition (`13∣v_max`) and, via mediants, a chain of increasingly weak
  *near*-divisibility conditions (`min_v‖jv‖_{13j+1}≤j`), none of which constrains the 12 core speeds.
- The medium-modulus witnesses that actually beat non-AP configs (S102, e.g. `q=24`) are at
  **non-convergent** denominators of `t*` — CF descent of `t*` does not even see them.

So CF descent is **insufficient** for the AP core, exactly as the sieve (S102) and maximality (S101) were.

## The synthesis (S101–S103): three tools, one layer

| tool | reaches | proved insufficient |
|---|---|---|
| **maximality / perturbation** (S101) | active-runner pinning at `t*` | interior small gaps are invisible to perturbation |
| **sieve completeness** (S102) | `q'∣` some speed for `q'≤13` | sieve-complete families can have `M>1/13` (beaten at `q'>13`) |
| **CF / Farey descent** (S103) | `lcm(13,14)∣v_max` (convergent+mediant) | reaches only far-element divisibility; core untouched |

**All three elementary tools converge on the same layer** — the far element's `13/14`-divisibility — and
**provably stop there**. The remaining content, the AP core (the 12 core speeds forming a dilated AP), is
**irreducibly additive**: it is the inverse theorem `M<1/13 ⟹ V∖{v_max}` is a dilated AP, equivalent to
Tao's `n=12` optimistic conjecture (S89/S94), and to the project's `≥6`-linear / additive-dimension-2
content (klein-S279, S92). It is **open**.

## Honest frontier statement

Across S96–S103 the proof of LRC(14) is now cleanly partitioned:

- **Everything elementary is done or discharged.** Non-covering ⟹ sieve witness (proved). The density
  route ⟹ closes for separated far elements via `|Error|≤κ'R_G/w`, deep well included (S100). The far
  element's structure ⟹ `lcm(13,14)∣v_max` by three independent routes (S101–S103). The full
  witness/descent machinery (`sieve_frac`, `fill1`, `descent`, `dilated_sieve`, THM-1017 elementary half)
  is kernel-pure Lean.
- **One thing remains, and it is the open inverse theorem.** `M<1/13 ⟹` the core is a dilated AP — a
  sharp additive-combinatorics statement that provably cannot be reached by maximality, sieving, or CF
  descent. Closing it is a research-level additive breakthrough (Tao `n=12`), not a further reduction.

I did not prove the crux and did not fabricate a proof. This session contributes the clean CF-descent
explanation of `182=lcm(13,14)` and the meta-result that the three natural elementary tools are exhausted,
pinning the open content to the additive inverse theorem.

Cross-links:
[[the-free-sieve-window-lever-is-refuted-the-crux-forces-at-medium-moduli-not-the-sieve-boxeph-S102]],
[[the-route-B-crux-is-the-open-inverse-theorem-what-covering-gives-and-why-maximality-cannot-finish-boxeph-S101]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
[[the-abandoned-attempts-decoded-the-crux-is-additive-dimension-two-not-any-scalar-boxeph-S92]].
