---
source: monad-explorer-2026-06-06-S709 (deep-research, signed→unsigned structural-reduction lane)
status: RIGOROUS structural reduction. Recalls canon (THM-369 n-grid witness; THM-398 LRC⟺C′;
  T1 gauge). NEW synthesis (HYP-2286): the signed/2n−1 apparatus is provably BLIND to C′ — the
  open core of LRC — for three independent, compounding reasons (modulus transversality; domain
  confinement to the worry-set; gauge-blindness of C′). Sharpens HYP-2281 (S708b: "sign group
  adds no M-content") by aiming it at the open conjecture, and orients where a measure-side tool
  must live.
tags: [signed-lrc, lrc, C-prime, THM-398, THM-369, gauge, T1, transversality, 2n-1, mod-n,
  worry-set, construction-vs-measure, open-core, structural-reduction, credits-HYP-2281, HYP-2286]
---

# The signed LRC is blind to the open core C′

**Dispatched angle.** *Structural reduction — how do signed-LRC witnesses/tight configs transfer
to the UNSIGNED LRC; what does a signed counterexample or tight configuration imply for the
unsigned problem?*

**One-line answer.** For the part of the LRC that is **already settled**, everything transfers
trivially; for the part that is **open**, the signed apparatus carries **no direct information at
all** — and this is a theorem, not an impression. The open core is gauge-invariant and lives on a
grid (mod `n`) provably transverse to the signed shell grid (mod `2n−1`).

## The two halves of the LRC, recalled (canon)

For a primitive speed set `S = {v_1,…,v_{n-1}}`, `M(S)=max_t min_i‖v_i t‖`; LRC(n) is `M(S)≥1/n`.

- **THM-369 (PROVED, Lean).** If `n ∤ v` for every `v∈S`, the time `t=k/n` already witnesses
  loneliness: `‖v·(k/n)‖ = min(r,n−r)/n ≥ 1/n` where `r=v mod n ∈{1,…,n−1}`. So the entire
  **construction side** (`n∤v`) is settled.
- **THM-398 (PROVED).** LRC(n) reduces *entirely* to

  > **C′ (CONJECTURE, the open core):** if `n | v` for some `v∈S`, then `M(S) > 1/n` (`S` is loose).

  i.e. a multiple of `n` forces *strict* looseness. (S564: every **tight** config has `n∤v`.)
- **T1 (S699 gauge).** `M({ε_i v_i}) = M(S)` for all signs `ε∈{±1}^{n-1}`: `M` is **sign-blind**.

So LRC's difficulty is concentrated **exactly** on the **measure side** `n|v`, where the `n`-grid
fails (the runner `v≡0 (mod n)` sits on the observer at every `t=k/n`) and loneliness must be
certified **off-grid**.

## The result: three reasons the signed/2n−1 theory cannot touch C′

All three are verified in `signed_lrc_blind_to_Cprime_s709d.py`.

> **(1) Modulus transversality.** The floor / C′ lives on the grid `G_n={k/n}` (THM-369). The
> signed-LRC shells live mod `C=2n−1` (T3/THM-401). Since `gcd(n,2n−1)=1` (verified `n≤1999`),
> `G_n ∩ G_{2n−1} = {0}` in `ℝ/ℤ` (verified `n≤59`). The two grids meet only at `0`. Moreover a
> **signed shell-partner** (`v_a+v_b ≡ 0 mod 2n−1`) is **never** the **floor-active pair**
> (`v_a+v_b ≡ 0 mod n`, canonically `(1,n−1)`): being both needs `v_a+v_b ≡ 0 mod n(2n−1)`,
> impossible for worry-set speeds (`0` such pairs over all tight configs `n≤8`).

> **(2) Domain confinement.** The signed theory's object of study is the **worry-set** (tight
> configs, e.g. `AP_n` and `V*`). Every tight config satisfies `n∤v` (S564; reverified: `0/9` tight
> configs `n≤8` carry a multiple of `n`). So the signed apparatus lives **entirely on the
> construction side**, the half *already* closed by THM-369. It never visits the `n|v` measure side
> that C′ is about.

> **(3) Gauge-blindness of C′.** C′ is a statement about `M`. By T1 the sign group `(ℤ/2)^{n-1}`
> fixes `M`, hence fixes C′'s truth value: every signed observable (sign-orbit size, shell-partner
> count, homometry class) varies across the sign group while the quantity C′ asks about does not
> move. Verified: on explicit C′ instances `(8,1,2,3,5,6,7)`, `(6,1,2,4,5)`, `(4,1,2)` the **full**
> `2^{n-1}` sign-orbit gives a single `M` value (`2/13, 1/5, 1/3` resp., all `>1/n` = loose), and
> `1418/1418` random sign-flips preserve `M`.

**Therefore** the signed/`2n−1` apparatus — gauge invariance, the cut/shell-partner structure,
the homometry/Patterson reframing (THM-413/415/417, HYP-2273), the silent-move chain calculus
(HYP-2280/2281) — carries **no direct information about C′**, the open core of LRC. What it *does*
refine is the **second Farey gap**: `2/(2n−1)` is the immediate Farey successor of `1/n`
(THM-401(i)), so the `2n−1` shells govern the loneliness landscape **one pinch above** the floor,
on already-settled construction-side configs. (This sharpens HYP-2281's "the sign group adds no
M-content" by pointing it at the open conjecture: not only does it add no `M`-content on the
worry-set, it provably cannot reach C′.)

## What this says about "signed counterexample ⟹ unsigned counterexample"

By T1 a *signed* counterexample (`M<1/n` for some `{ε_i v_i}`) is *literally* an unsigned one with
identical `M`. And **every** counterexample (signed or not) must, by THM-398, carry a multiple of
`n` — i.e. live on the measure side, exactly where the signed refinement is blind. So the transfer
in the counterexample direction is real but degenerate: signs are a gauge, and a gauge cannot
manufacture or detect the obstruction; the obstruction is the `n|v`, off-grid, sign-invariant
content.

## The forward handoff: a *measure-side* tool must live mod n and be gauge-blind

The honest consequence is a redirection, not a dead end. To bite **C′** one needs an instrument
with the right three properties, opposite to the signed one:

1. **mod n, not mod 2n−1** — C′ is an `n`-grid obstruction.
2. **on the `n|v` side** — the worry-set/tight analysis (where the signed theory lives) is the
   wrong domain; the relevant configs are the off-grid ones with a multiple of `n`.
3. **gauge-blind by design** — since `M` is sign-invariant, any sign-orbit-varying quantity is
   measuring noise relative to C′.

A concrete candidate worth probing (open): the off-grid loneliness of a C′ config is certified at
**pair-sum pinch times** `t=m/(v_a+v_b)` (THM-401's pinch lemma). The sum-clocks `v_a+v_b` are the
*only* signed datum that is **not** mod-`2n−1`-shell noise — they are the off-grid denominators
themselves. So the productive reading of "signed" for C′ is **not** the homometry mod `2n−1`, but
the bare **sum-set `{v_a+v_b}` reduced mod `n`** of a multiple-of-`n` config: does some pair-sum
avoid the obstruction created by the `n|v` runner? That is the measure-side, mod-`n`, gauge-blind
question THM-398 actually asks, and it is untouched by the present signed program.

## Honest status

- **Recalled canon (PROVED):** THM-369 (`n∤v ⇒ M≥1/n` via `t=k/n`); THM-398 (LRC ⟺ C′); T1.
- **PROVED here:** `gcd(n,2n−1)=1` ⇒ `G_n ⊥ G_{2n−1}`; floor-active and shell-partner pairs are
  disjoint on the worry-set; C′ is gauge-invariant (immediate from T1).
- **VERIFIED:** transversality (`n≤59`), domain confinement (`0/9` tight configs `n≤8` have `n|v`),
  gauge-constancy of `M` on C′ instances (full `2^{n-1}` orbits + `1418` random checks).
- **NOT proved (and not claimed):** C′ itself (the open core); that *no conceivable* sign-based
  device helps (only that the **present** signed observables, being `M`-functions across a gauge,
  cannot — a precise and, I believe, decisive obstruction for the program as formulated).
- **CONJECTURE/handoff (HYP-2286):** the measure-side instrument for C′ is the pair-sum-set mod `n`
  of `n|v` configs, gauge-blind; the `2n−1` homometry program is provably orthogonal to it.

**Relation to neighbours.** Credits and sharpens **HYP-2281** (S708b: shell-partner is unsigned,
sign-orbit-constant). Uses **THM-369/THM-398/S564** (the `n|v` Vitali handoff; see
`lrc-the-vitali-handoff-equation-is-multiple-of-n-s573.md`) and **THM-401** (the `2n−1` second-gap
modulus). Complements the homometry view (`signed-lrc-as-cyclic-homometry-s705.md`), which already
flagged "the observer gap `M` is untouched by design" — this note proves *why* it must be, and
that the untouchedness extends to the entire open core.

**Artifacts:** `04-computation/signed_lrc_unsigned_transfer_s709.py`,
`lrc_ngrid_floor_completeness_s709b.py`, `lrc_offgrid_tight_census_s709c.py`,
`signed_lrc_blind_to_Cprime_s709d.py` (+ `05-knowledge/results/*.out`). HYP-2286.
