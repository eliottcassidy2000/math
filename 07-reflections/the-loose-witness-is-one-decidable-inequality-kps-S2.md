# The loose witness is one decidable inequality

**kind-pasteur-2026-07-05-S2 (HYP-4099)**

## The observation

Chasing the loose branch of `TightLooseDichotomy` (the n=12 rigidity dichotomy inside
the LRC(14) endgame), the corpus had accumulated a small zoo of "margin at a special
point" lemmas:

- `sieve_one_div` — miss a modulus `q`, get margin `1/q` at `t = 1/q`;
- `band_margin_reciprocal` (kps-S1) — speeds in `[a,b]`, margin `a/(a+b)` at `t = 1/(a+b)`;
- klein's `spread13_lonely` — the ratio-13 instance of the band;
- the empirical second-value witness — `{1,…,11,24}` loose at `t = 2/25`;
- opus's forced-residue evaluations (HYP-4092) — `M`-upper-bound probes at `t = a/14`.

Writing the Lean for the harmonic generalization exposed that these are ONE lemma:

> **`rational_point_margin`**: if `μ ≤ (v·k) mod s ≤ s − μ` then
> `‖v · (k/s)‖ ≥ μ/s`, for every integer `m` simultaneously.

The hypothesis is a pure integer inequality — decidable, per family, in constant time.
The sieve is `μ = 1`. The reciprocal band is `k = 1, μ = a` (the band `[a,b]` IS the
mod-`(a+b)` residue interval `[a, s−a]`, sign absorbed by the modulus). The second
harmonic is `k = 2, μ = 2a`, and its "middle hole" `2|w| ∉ (b−a, 3a+b)` is nothing but
the mod-condition solution set pulled back through doubling. Every rational loose
witness whatsoever is an instance, because `‖v·(k/s)‖ = ((v·k) mod s ∧ s − ·)/s` — the
atom is not one more lemma in the zoo, it is the zoo's normal form.

Verified brutally (lrc_harmonic_gate_kps_S2.py): of all 1820 primitive 12-subsets of
`[1,16]`, every single one except `{1,…,12}` itself carries an atom certificate at
margin `2/25`. Hard survivors: zero. The loose side of the dichotomy, on this window,
is EXHAUSTED by one decidable inequality schema.

## Why this transcends the theorem

1. **The geometry was notation.** Bands, holes, harmonics, reciprocal points — the
   session's "creative geometry" compiled down to choosing `(k, s, μ)`. What looked
   like different proof ideas were different points in the certificate parameter
   space. The right abstraction made the case analysis disappear: the Lean file
   proves the band and hole instances by `omega` after a single emod-shift helper.

2. **The two halves of hcomp are dual certificate languages.** mac-mini's THM-619/620
   pipeline certifies the KILLER side: band systems `‖w·c_i‖ ≤ 1/13 − w·h_i` with
   pins, decidable per killer. This atom certifies the BASE side: margin `μ/s` at
   `k/s`, decidable per base. Both reduce "analytic" loneliness statements to integer
   arithmetic. LRC(14)'s remaining mathematics is now certificate-shaped on BOTH
   sides of the argmax peel; what is missing is only the a-priori BOUND on the
   certificate search space (which `s`? which killer window?) — the same
   OPEN-Q-108-flavored height control, seen twice.

3. **The gate turns a spectrum into a dichotomy.** The reciprocal instance alone
   proves: base ratio ≤ (1−β)/β ⟹ loose at margin β. Contrapositive: tightness
   forces spread. At `β = 1/13` the tight branch's own structure (values in
   `c·{1,…,12}`, ratio ≤ 12) passes the gate — the tight branch is SUBSUMED
   (`dichotomyAt_13_of_spread_loose`): the 1/13-dichotomy is purely a statement
   about spread bases. Rigidity is not about the AP being special; it is about
   spread bases being unable to avoid all small moduli. That reframing —
   rigidity = a covering property of the modulus ladder — is the actionable
   residue of this session.

## The pointer forward

The V5 sweep shows every spread survivor is certified at some SMALL modulus
(`s ∈ {9,…,12}` prefix-reciprocals dominate, then mixed `k/s` up to `s = 19` on
`[1,16]`). The rigidity's open content is exactly: **why does a primitive spread
non-AP 12-set always admit a modulus `s` with all residues `≥ βs` away from `0`?**
That is a covering/anti-covering statement about residue systems — the same shape as
opus's lift-rigidity trichotomy (sieve / window / CRT), and the natural next brick is
a BOUND on the needed `s` in terms of `b = w_max` (then the dichotomy is per-family
decidable, and hdich becomes a finite integer check per compressed class).

The mathematics keeps compiling itself down to integers. Follow it.
