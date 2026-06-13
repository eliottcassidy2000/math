---
id: THM-417
title: The signed face IS the symmetric-closure / gauge-orbit reading of the pair-sum sieve — it refines but cannot move the gap
status: PROVED (elementary, parts 1-4) + the AP-global-energy-max is VERIFIED n<=12 (not yet proved)
source: monad-explorer-2026-06-06-S707
depends_on:
  - THM-401   # pair-sum sieve modulus C=2n-1, summand shells P_a={a,-a}
  - THM-403   # AP floor witness = primitive n-th roots ((Z/n)^x orbit)
  - THM-414   # signed additive face = multiplicative energy; E_+, r_+(0) dilation-invariant
  - HYP-2262  # signed LRC theory T1 (gauge invariance of M), T2 (sign=cut), T3 (zero-clock<=>shell-partner)
related:
  - HYP-2272  # is (r_+(0),E_+) a complete worry-set invariant (CORRECTED here)
  - THM-407   # twisted-involution shell reduction
---

# THM-417 — the signed face is the symmetric closure of the pair-sum sieve

**Angle (dispatched):** structural reduction — how the *signed*-LRC additive/shell
structure (mod `C = 2n−1`) transfers to the *regular* (unsigned) LRC gap
`M(S) = max_t min_i ‖v_i t‖` (which lives mod `n`, THM-403). The result is a clean
identification of three a-priori-different objects, plus a precise statement of what the
signed framework can and cannot do for the unsigned problem.

## Setup

`AP = {1,…,n−1}` (`k = n−1` runners, floor `δ = 1/n`). `C = 2n−1`. THM-401 fixes the
pair-sum sieve modulus at `C` with summand shells `P_a = {a, C−a} = {a, −a}`,
`a = 1,…,n−1`. The signed LRC (HYP-2262) gauges the speeds by signs `ε ∈ {±1}^k`; T1
says `M` is gauge-blind, so the content is the *pairwise* (additive) data exposed by the
sign group.

## Theorem

**(1) Modulus = size of the symmetric closure.** Let
`sym(AP) = {0, ±1, …, ±(n−1)}`. Then `|sym(AP)| = 2n−1 = C`, and reduced mod `C`,
`sym(AP) = ℤ/C` (the full group). *Proof:* the `2(n−1)+1 = 2n−1` integers
`−(n−1),…,n−1` are pairwise incongruent mod `2n−1` and there are exactly `2n−1` of them,
so they are a complete residue system. ∎ **The pair-sum sieve modulus `2n−1` is exactly
the cardinality of the symmetric closure of the AP** — the `(2n−1)`-th-roots witness orbit
of THM-403's complement is `sym(AP)` filling `ℤ/C`. (Verified n=3..9: `sym(AP)=ℤ/C`,
`#shells = n−1`.)

**(2) Shells = sign-gauge orbits.** The summand shells `P_a = {a, −a}` are *exactly* the
orbits of the sign-gauge group `{±1}` acting on `(ℤ/C)∖{0}` by `x ↦ −x`. *Proof:* the
orbit of `a≠0` under `{±1}` is `{a,−a}=P_a`; `C` odd ⟹ `a≢−a` ⟹ every orbit has size 2;
there are `(C−1)/2 = n−1` of them. ∎ So **the gauge group of the *signed* LRC (T1) and the
summand-shell structure of the *additive* sieve (THM-401) are one and the same `{±1}`
action**, and `r_+(0)` (the signed zero-clock count, T3) `=` the number of shell
*collisions* of the config (pairs `{i,j}` with `v_i ≡ −v_j`, i.e. both elements in one
shell).

**(3) Symmetrization bridge (sum-face of `S` = difference-face of `S∪(−S)`).** For any set
`S`,
```
   Sum(S) := {v_i + v_j}  =  {v_i − (−v_j)}  ⊆  Diff(S ∪ (−S)),
```
and `{i,j}` is a **shell-partner** (`v_i+v_j ≡ 0`) `⟺` `v_i = −v_j` in `ℤ/C`
`⟺` the `+`-copy of `i` and the `−`-copy of `j` **collide** in `S∪(−S) ⊂ ℤ/C`. *Proof:*
direct. ∎ Hence the *signed* pair-sum structure at parameter `n` is the *unsigned*
difference (relative-speed) structure of the doubled symmetric config of size `2n−1`, and
the shell-partners are precisely its `±`-collisions. **AP is shell-perfect** (its symmetric
closure tiles `ℤ/C` with no collision, by (1)); a shell-partner is a defect of that
tiling (one shell doubled, one empty — e.g. `V*` at `n=14`: `3+24≡0 mod 27`, shell `{3,24}`
doubled, shell `{12,15}` empty).

**(4) The two-clock reduction (what transfers, and the limitation).** There are two
coprime clocks (`gcd(n, 2n−1)=1`):
- the **witness clock mod `n`** governs the gap `M` (THM-403: a witness `t=j/n` requires
  only `v_i ≢ 0 mod n`); it is what the *unsigned* LRC sees.
- the **shell clock mod `2n−1`** carries the signed additive face `(E_+, r_+(0))`.

Every signed additive invariant is **sign-gauge invariant** (T1) and **dilation-invariant**
(THM-414), hence descends to a function on unsigned speed-sets mod `(ℤ/C)^× ⋊ {±1}`.
Therefore:

> **The signed face is a genuine invariant of the *unsigned* config, but it lies in the
> kernel of "determining `M`": by T1 it cannot change the gap. The signed framework can
> only REFINE / classify gap-equivalent configs — it can never, by itself, bound `M`.**

This is confirmed sharply at the `n=14` floor: `AP`, `2·AP`, `V*` are **all tight**
(`M = 1/14`), yet the shell clock separates them (`r_+(0) = 0,0,1`; `E_+ = 328,328,290`)
while the witness clock mod 14 is blind to the `V*` shell collision (it only forbids
residue 0). The signed face is *strictly finer than* `M`, and *strictly weaker than* `M`
for the purpose of proving loneliness — a refinement, not a reduction of the gap.

## Corrected corollary (sharpening HYP-2272)

**AP is the GLOBAL additive-energy maximizer**, not merely the worry-set maximizer:
`E_+(S) ≤ E_+(AP)` for **every** `(n−1)`-subset `S` of `(ℤ/C)∖{0}` (VERIFIED exhaustively
`n = 4..12`, `C = 7..23`, **0 exceptions**; an instance of the classical
"interval maximizes additive energy"). Two consequences correct the reading of HYP-2272:
- the inequality "`tight ⟹ E_+ ≤ E_+(AP)`" is a **general fact**, true for all configs, so
  it is **not** a worry-set separator;
- the equality case is **not** "dilate of AP": the global-energy maximizer set is far
  larger than the AP dilation-orbit (e.g. `n=12`: 132 maximizers vs 22 in the AP-orbit).

The genuine worry-set content is therefore the **conjunction**: among *tight* configs, the
energy-maximal ones appear to be exactly the AP-dilates (verified small `n`: at `n=5,6` the
unique tight energy-maximizer is the AP; sporadic tight configs sit strictly below). This
is restated as the sharpened HYP-2272 (the inequality is general; tightness is what makes
equality characterize AP).

## Honest scope

- Parts (1)–(4) are **proved** (elementary) and verified `n ≤ 9` (Part A/D of the script).
- "AP is the global `E_+`-maximizer" is **VERIFIED `n ≤ 12`, not proved.** A general proof
  (interval-maximizes-additive-energy in `ℤ/(2k+1)` at the half-size `k=(C−1)/2`) is left
  open. Per repo discipline (patterns break at `n≥7/8`) this is recorded as VERIFIED.
- This is **NOT** a proof of `LRC(14)` or any new bound on `M`. By part (4) the signed face
  cannot, in principle, supply such a bound; it organizes the worry-set.

**Artifacts:** `04-computation/signed_lrc_shell_reduction_s707.py` (+`.out`),
`04-computation/signed_lrc_energy_globalmax_s707b.py` (+`.out`),
`07-reflections/the-signed-face-is-the-symmetric-closure-of-the-sieve-s707.md`.
Builds on THM-401/403/414, HYP-2262 (T1–T3), THM-407.
