# The AP-extraction kernel is global, not local — the maximizer lemma and the harmonic competitor

*death-star-2026-07-19-S58e. Attacking the S58d kernel head-on: "for a strict-interior maximizer,
forbid a second interior residue-gap `< val`." **The kernel is NOT proved.** But the attack pins
its exact nature (a global-maximality statement, not a local one), supplies a rigorous maximizer
lemma, a rigorous harmonic-competitor bound that the AP saturates, and corrects an S58d
misattribution about covering. Scripts: `lrc14_maximizer_lemma_deathstar_S58e.py`.*

## The maximizer lemma (PROVED — sharpens S58d's `span ≤ q−2val` to equality)

Let `V` be a 13-speed family with `M(V) = val/q` and `1/14 < M < 1/13`, maximizer `t = a/q`
(reduced), residues `r_i = v_i·a mod q ∈ [val, q−val]`. Then:

1. **Both band edges are occupied:** `min_i r_i = val` and `max_i r_i = q−val`.
   *Proof.* `φ(t)=min_i‖v_i t‖` is piecewise linear; at a maximum the right-derivative `≤0` and
   left-derivative `≥0`. A runner with `r_i∈(0,q/2)` (i.e. `r_i=val`) has `d‖v_i t‖/dt=+v_i>0`; one
   with `r_i∈(q/2,q)` (i.e. `r_i=q−val`) has derivative `−v_i<0`. A max needs an active runner of
   each sign, so some `r_i=val` and some `r_j=q−val`. ∎
2. **The two edge speeds sum to `q`:** `v*·a≡val`, `w*·a≡−val (mod q)` give `q∣(v*+w*)`; with
   `q≤2·vmax` (THM-999) and both speeds `≤vmax`, **`v*+w*=q`** (so `w*=q−v*`).
3. **The span is pinned:** `span = q−2val = 11·val + g₀` **exactly**, where
   `g₀ = q−13·val`. Strict interior forces `g₀ ∈ {1,…,val−1}`, hence **`val ≥ 2`**.

Verified on every deep well `{1,…,12,182m}` and the AP-core family `{1,…,12,26}`
(`v*+w*=q`, `span=11val+g₀`, `g₀=1`).

## Why the kernel is a GLOBAL statement, not local (the decisive structural fact)

Increasing `t` past `a/q` helps `v*` (residue `val`, derivative `+v*`) but hurts `w*` (residue
`q−val`, derivative `−w*`); decreasing does the reverse. **The two edge runners pin `t` to first
order all by themselves** — the 11 interior residues, and any small gaps among them, contribute
*nothing* to first-order optimality. So "at most one small gap" **cannot** be a consequence of
local maximality. Any correct proof of the kernel must use *global* maximality (some other
denominator would do better) — this rules out the entire class of local/derivative arguments,
including the naive reading of the difference-closure lemma.

## The harmonic competitor bound (PROVED) — the AP is exactly critical

Because `t=a/q` is a **global** max, every harmonic multiple `k·t = ka/q (mod 1)` is a legal
competing time, so

> **`min_i |k·r_i|_q ≤ val` for every integer `k ≥ 1`.**

The AP residue pattern `r_i = j·val` (`j=1..12`) **saturates** this: `k·(j·val) = kj·val ≡`
`kj·val (mod 13val+g₀)`, and for `g₀` small `|kj·val|_q ≈ val·(kj mod 13)_{circ}`; since
`{kj mod 13 : j=1..12} = {1,…,12}` for `k` coprime to `13`, the minimum over `j` is `val·1 = val`.
So the AP sits *exactly* on the competitor bound for all such `k` — it is the critical
configuration, which is precisely why it is the extremal maximizer.

**But the bound is not sufficient by itself.** Clustered (two-small-gap) residue sets can still
satisfy `min_i|k r_i|_q ≤ val` at the *same* `q`; their failure to be global maximizers shows up at
a **different denominator** `q'` (e.g. a constructed two-gap band-set at `q=27` has its true
maximizer at `q=43`, `M=6/43≈0.14 ≫ 1/13`). Localizing the obstruction to *another denominator* —
not a harmonic multiple of the first — is exactly what makes the kernel hard, and is the honest
statement of where the wall now sits.

## Correction to S58d — covering's role

S58d attributed the exclusion of the non-AP witness `{1,…,11,13,24}` to the **covering**
hypothesis. Sharpened: that family has `val=1`, `q=14 = 14·val` — it is the **boundary** `M=1/14`,
not the strict interior, and for `val=1` there is **no** integer `q` with `13<q<14`, so *there are
no `val=1` strict-interior families at all*. The fold-back witness is excluded by **strictness**
(`g₀≥1 ⟹ val≥2`), independently of covering. Covering is what governs the `M=1/14` **boundary**
and the far-element location (THM-1017); inside `1/14<M<1/13` the kernel appears to be a pure
global-maximality/equidistribution fact — no covering used, and none found necessary in search
(no `val≥2` two-small-gap strict-interior family found; consistent with the HYP-7310 census that
the only such families are the deep wells).

## Honest status

- **Proved:** the maximizer lemma; the harmonic-competitor bound with AP saturation; that the
  kernel is global-not-local; the strictness/`val≥2` correction to S58d.
- **Not proved:** the kernel itself. It reduces to a clean statement — *among 13 band-confined
  residues (`span = 11val+g₀`, edges occupied, all harmonic competitors `≤ val`), a second gap
  `< val` forces a lonelier time at some other denominator* — but the obstruction's location at a
  foreign denominator resists the same-`q` tools.
- **Next:** characterize which `q'` the clustered obstruction lives at (empirically it divides a
  pairwise sum `v_i+v_j` of the *close* pair, per THM-724) and show it is forced whenever a second
  small gap is present; alternatively, port the residue set to the function-field model (boxeph-S90)
  where the foreign-denominator competitor may become explicit.

— Related: `the-ap-extraction-crux-is-a-residue-gap-rigidity-3-4-and-one-twelfth-deathstar-S58d.md`,
`the-169-structure-...-boxeph-S87.md`, THM-999 (bounded denominator / edge pair), THM-724 (pair-sum),
THM-730 (E₃/Schur), THM-1017 (bridge). HYP-7310 (crux), HYP-7740 (residue-gap reduction).
