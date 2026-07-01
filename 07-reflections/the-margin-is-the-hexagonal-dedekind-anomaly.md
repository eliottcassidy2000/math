# The LRC14 margin is the hexagonal Dedekind anomaly — and B2 is the road not taken

*mac-mini-2026-06-30-S64. Prompted by "consider B2/E2/Dedekind objects and work on the hard residual."*

Three classical objects were named — the B2 root system, the weight-2 Eisenstein series E2, the Dedekind
eta/sums — and asked to touch the hard residual. They turned out to be one object, and it is the **margin**.

The covering-min of LRC14 is `14/183 = n/Phi_6`, and it beats the conjectured floor `1/14` by a margin
`13/2562`. That margin has always looked like a small arithmetic accident. It is not. It is *exactly* a
Dedekind sum:

    margin(n) = n/Phi_6 - 1/n = -12 · s(n, Phi_6) / n²,   Phi_6 = n² - n + 1.

The Dedekind sum `s(n, Phi_6)` has a clean closed form, `12 s = -(Phi_6 - 1)/Phi_6`, and the *reason* is that
`n^3 ≡ -1 (mod Phi_6)`: the observer speed `n` is the **order-6 element**, the sixth root of unity `omega`,
the 60-degree rotation of the Eisenstein/hexagonal lattice that the whole "everything is the triangle" picture
has been circling. The margin is a special value of the Dedekind sum at `omega`.

Now B2 earns its place by **contrast**. Change the observer's order from 6 to 4 — from the hexagonal `omega`
(60°) to the square `i` (90°), from A2 to B2, from `h^3 = -1` to `h^2 = -1` — and the Dedekind sum
**vanishes**: `s(i, k) = 0` for every order-4 element. A square-lattice world has *zero margin*: its
covering-min would sit exactly on the floor, with no room, degenerate. The reason LRC14 holds with room to
spare is that its geometry is hexagonal, not square — and "hexagonal, not square" is precisely "`s ≠ 0`",
the non-vanishing of the Dedekind anomaly. B2 is the road not taken; on it, the conjecture would be tight and
airless.

And E2 is where the anomaly lives as a modular object. `E2(τ) - 14 E2(14τ)` is a weight-2 Eisenstein form on
`Gamma_0(14)`; its quasimodular anomaly *is* the Dedekind sum, and its cusp data carries the floor. Three of
its cusps are Eisenstein (the computable bulk); the fourth direction — the genus-1 cusp form `f_14` at the
apex cusp `d=7` — is the obstruction the program has named for a year. So the chain closes on itself:

    the margin  =  s(n, Phi_6)  =  the E2 anomaly on Gamma_0(14)  =  positive because A2, zero on B2.

This reframes the hard residual in one sentence. The residual is: a covering set missing a core speed, patched
by a large `±k`, must not erase the margin. In the new language, **erasing the margin means driving the
hexagonal Dedekind anomaly to zero — restoring the square/B2 world.** A patch is a single integer; its
residues across moduli are CRT-linked; and Dedekind–Rademacher reciprocity (the repo's own far-coherence tool)
forbids its Dedekind sums from going anomaly-free at every modulus at once. The anomaly is *rigid*. That
rigidity is "the hole moves but never vanishes," and it is why B2 cannot be reached from A2 by a finite patch.

The lesson that keeps recurring: every time we reduce LRC14 to a smaller gadget, the gadget turns out to be a
known modular/arithmetic invariant we had already been staring at from another side. The margin was never an
accident. It was `s(omega)` — the anomaly that separates the hexagon from the square.

*See [[HYP-3768]], [[HYP-3715]] (the zeta_6 hexagonal line), [[HYP-3586]] (the apex cusp / f_14), [[HYP-2808]]
(Dedekind–Rademacher reciprocity), and [[everything-is-the-triangle]].*
