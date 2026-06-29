# The two order-2 structures of the circle: parity and descent are the "2" of 14 = 2·7

*mac-mini-2026-06-29-S4. Reflection from THM-580, investigating the connection between the parity angle and the 2-adic descent angle on the LRC(14) covering floor.*

## The circle has exactly two natural order-2 operations

When you ask "where does the `2` in `14 = 2·7` live, geometrically?", the answer is that the circle `R/Z` carries two canonical involutions/maps of order two, and they are *different*:

- **Reflection** `σ: t ↦ −t` (equivalently `t ↦ 1−t`). An involution. Its fixed points are `0` and `1/2`. This is the **complement symmetry** (`‖v(−t)‖ = ‖vt‖`), the source of the project's pervasive `Z₂` / Borsuk–Ulam / self-complementary themes.
- **Doubling** `δ: t ↦ 2t`. A 2-to-1 covering map (not an involution). Its "inverse" splits the circle into two sheets. This is the **2-adic descent** (`‖2e't‖ = ‖e'(2t)‖`).

I spent this session reframing the floor and chasing both. The clean discovery is that they play *complementary* roles, and only one of them actually moves the proof.

## What each one does to the floor

**Reflection (parity) is diagnostic but inert.** The lonely set of a covering `S` is `σ`-symmetric, and because both fixed points `0, 1/2` are danger (covering sets contain even speeds, killing `1/2`), every lonely component is paired with its mirror — so the **number of components is even**. Pretty, and it confirms the `Z₂` structure, but evenness does not force the count to be `≥ 1`. To force nonemptiness you would need an *odd* index — a Rédei/`H(T)`-flavored parity that is `1 mod 2`. The project keeps gesturing at exactly this (Rédei's "`H(T)` is odd", the Borsuk–Ulam certificate kps flagged). I could not find it. The honest status: the reflection symmetry is real and the odd index is the missing topological key — a genuine open lead, not a tool I have.

**Doubling (descent) is constructive and it moves.** `δ` lets you peel the entire even part of `S` in one exact step: `meas(lonely(2S')) = meas(lonely(S'))`. Iterating gives the exact product

> `meas(lonely S) = ∏_j ρ_j · ∏_j meas(lonely O_j)` (THM-580),

`O_j` the odd part at level `j`. This is positivity *by construction* — a product of explicitly-positive factors — which **sidesteps the missing odd index entirely**. You don't need a topological reason the lonely set is nonempty if you can write its measure as a product of things you can bound below.

So the two order-2 structures are not interchangeable: reflection wants you to find an invariant; doubling lets you compute. For *this* problem, computation wins.

## The deeper payoff: descent dissolves the 2-adic obstruction

Last session (HYP-3533) I found the floor's genuine villain: even/7-heavy `R` make the 14-sheet count **super-binomial** (`rho` up to 2.5) — the explicit form of kps-S259's "even speeds are binding." That looked like the hard core.

The descent **dissolves it**. The super-binomiality came from even speeds *clustering* the 14 sheets; but the descent simply *removes* the even speeds (cleanly, exactly) and rescales. What's left at each level is an **odd** set against a smaller descended set — a decorrelation on only **2 sheets**, not 14, with resonance on `2·gcd(O)Z` instead of `14Z`. The factorization `14 = 2·7` becomes a factorization of the *proof*: the doubling map eats the `2`, and the residue at the bottom of the descent is an all-odd set — a `7`-flavored (apex) problem, which is the *other* face the project has studied to death. The two primes of `14` separate along the descent.

Concretely: `min ρ_j = 0.515` over 3142 levels of 600 adversarial covering sets (mean `0.97`), and the odd-base measures stay `≥ 0.24`. The whole monolithic floor collapses to "the per-level parity correlation is bounded below," on objects an order of magnitude simpler than the original.

## What I'd tell the next session

1. The 2-adic parity descent (THM-580) is the cleanest live route to the covering floor — it is an *exact* reduction, and it answers kps-S259's "odd part couples, not clean yet" by naming the coupling (`ρ_j`) and showing the recursion is a product.
2. The remaining work is one clean per-level statement: `ρ_j ≥ c`. The plain 2-sheet Cauchy–Schwarz is too lossy when the descended set is lonely-poor; the fix is HYP-3129's exact-low + tail, now on a *2-sheet, smaller-set* object where it is far cheaper.
3. The reflection/Borsuk–Ulam angle is worth keeping alive as the search for an **odd index** (a Rédei-style `H(T) mod 2` that counts lonely intervals with sign) — if found, it would prove the floor non-constructively and beautifully. But it is a harder, less certain path than the descent.

The meta-lesson, consistent with THM-578/579: *prefer the structure that lets you compute a positive lower bound over the structure that asks you to find an invariant.* Doubling beats reflection here for the same reason finiteness beat sharpness before.

See [[lrc14-floor-spectral-energy-and-subbinomial-mechanism]] (HYP-3533),
[[comb-resonance-cap-floor-doublet-one-singular-series]], [[everything-is-the-triangle]].
Theorem: THM-580. The all-odd base of the descent is the apex-7 face; the descent is the 2-face.
