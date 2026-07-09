# The 7-structured μ-floor is the Z/7-coloring problem, resurfacing

*mac-mini-2026-07-09-S64. Owner: find creative angles on the last residual (the dissociated
7-structured good-period μ-floor) and mine past work for connections. The residual turned out to
be an old friend — the Z/7 vertex-coloring the project already mapped (mac-mini-S6, HYP-2703) —
seen from the finite-`Vmax` side.*

## The residual, and why moments died

LRC(14)'s covering case is one lemma from done: for **dissociated** co-offset clusters
(`longest-AP ≤ k−6`), a good period exists (`#{good grid j} > 0`). The hard sub-case is the
**7-structured** sets — many co-offsets `≡ 0 mod 7`. There (MISTAKE-128, klein-S200) the arc-count
`c = #arcs/spread` spikes to `0.88`, the moment lower bounds `B_3..B_6 = 0.63..0.77` all stay below
it (S63: moments converge to `μ` too slowly on the resonance), so **no moment certificate** `c < B_d`
survives. klein-S200 named the right invariant — the *maxgap margin*, not the arc count — but the
a-priori `μ`-floor stayed open. The obstruction is arithmetic, and it is `7`.

## The connection: this is the coloring, from the other side

mac-mini-S6 / HYP-2703 already established that the density-floor object is a **Z/7 vertex-coloring**
`c(e,x) = ⌊7·frac(ex)⌋`, and that "does the cluster leave a `>1/7` gap" is *the arithmetic of
residues mod 7, not the geometry of gaps on the circle* — with a proven **slope-band decomposition**
`measS7(E) = (1/7) Σ_{s=0}^{6} bandcover_s(E)`, each band the same Sturmian cover twisted by the
multiplier `e ↦ s·e mod 7`. The seven matters *because it is the modulus*. The 7-structured clusters
are exactly the ones that are **arithmetically special for that coloring**: their `≡0 mod 7` elements
pin colour `0`, so they resonate with the `1/7` grid the coloring is built on. The finite-`Vmax`
μ-floor is the same Z/7 object, now asked on the ruler grid `{j/Vmax}` instead of the continuum.

## The mechanism at the resonances `x = m/7`

At `x = m/7` every phase lands on the 7-point grid: `frac(e_i·m/7) = (e_i·m mod 7)/7`. So the maxgap
there is decided by mod-7 arithmetic alone:

- **Missed residue.** If `{e_i mod 7} ≠ ℤ/7` — some residue `r` is absent — the slot `[r/7,(r+1)/7)` is
  empty, forcing a gap `≥ 2/7 > 1/7`. (Original MISTAKE-128 set: residues `{0,1,2,4,5}`, misses
  `{3,6}` ⟹ `maxgap = 2/7` at `x = 1/7`, wide-open good.)
- **÷7 collapse ⟹ the `k ≤ 7` pigeonhole (THM-530).** The `|S₇| = #{e_i ≡ 0 mod 7}` elements *all
  collapse to phase `0`* at `x = m/7`. If `|S₇| ≥ k−6`, the effective phase count drops to `k−|S₇|+1 ≤ 7`,
  and THM-530's `k ≤ 7` pigeonhole gives `maxgap ≥ 1/(k−|S₇|+1) > 1/7`. (klein's worst set:
  `|S₇| = 10 ≥ 7` ⟹ ≤ 5 effective phases ⟹ `maxgap ≥ 1/5`.)

Every hard 7-structured set satisfies **one** of these (verified): the missed-residue one when the
resonance is sparse, the collapse/pigeonhole one when it is dense. So at the *continuum* resonance
`x = m/7` the maxgap is not just `> 1/7` — it is `≥ 1/5` to `2/7`, a **huge** margin. This is klein's
"7-structure fragments Good_E without lowering the best maxgap," now with the arithmetic reason:
the fragments are the tiling generic `x`; the *resonances themselves are the widest good arcs*.

## The `gcd(7, Vmax)` split — where the grid meets the coloring

The good *period* needs a *grid* point in a good arc, so the arithmetic of `gcd(7, Vmax)` is the hinge:

- **`7 | Vmax`:** the resonance `x = m/7 = (m·Vmax/7)/Vmax` **is a grid point**, and the mechanism
  above makes it good (margin `≥ 1/5`). The 7-structure that broke the arc-count now *hands us the
  good period for free* — at the very point it resonates. (Verified: `m·Vmax/7` good for the vast
  majority; the exceptions cover all residues *and* have `|S₇| < k−6`, i.e. are not the hard
  high-`c` sets.)
- **`gcd(7, Vmax) = 1`:** no grid point sits on a resonance, so the grid *decorrelates* from the
  `1/7`-structure and samples the good bulk (`μ ≈ 0.94`); existence held on 100% of the sample.

So the 7-structure is not an obstruction to *existence* — it is the reason existence is *easy* when
`7 | Vmax` (the resonance is a wide good arc on the grid), and irrelevant when `gcd(7,Vmax)=1` (the
grid never sees it). The arc-count read the fragmentation and got scared; the arithmetic says the
best period is sitting on the resonance.

## What this buys, honestly

Not yet a closed proof — the `gcd(7,Vmax)=1` branch still leans on the bulk `μ` (the general
decorrelation), and the `7|Vmax` "covers-all-residues, `|S₇|<k−6`" corner needs the general argument
too. But it **reframes the last residual into the project's oldest and best-understood object** (the
Z/7 coloring, the slope-band decomposition, THM-530's `k≤7` pigeonhole, THM-536 Sturmian), splits it
by the one arithmetic quantity that matters (`gcd(7,Vmax)`), and turns the `7|Vmax` hard case — the
one that broke every generic tool — into a *one-line* pigeonhole at the resonance. The moment route
died because it ignored the `7`; the way through is to *use* it.

Follow the modulus.
