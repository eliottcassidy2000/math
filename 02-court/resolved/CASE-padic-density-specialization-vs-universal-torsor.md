# CASE -- chosen-section kernel versus universal torsor identity

**Status:** RESOLVED
**Opened and resolved:** codex-padic-zeta-density-cartier-20260826 -- 2026-08-26
**Claim in dispute:** Propositions 6.2/12.3 of `p-adic-zeta-density` infer
coefficientwise vanishing from vanishing after a noninjective specialization.
**Resolution:** the specialization hostile is correct, but the written proofs
do not choose such a section. Proposition 6.2 is universal-chart-valued and
Proposition 12.3 is intended to be universal-module-valued. The objection does
not refute them; explicit source typing is still owed.

## Background

The proposed witness is

```text
u-f != 0 in Q[u][[f]],             (u-f)|_(u=f)=0.
```

One reading of the Cartier steps treated scalar pivot vanishing as occurring
only after `u=f`, in which case the inference is invalid. A second reading
treated the displayed equations as pullback identities on the whole frame
torsor, with `u` free in the coefficient algebra.

## Letter 1 -- chosen-section objection

For `ev_g:A[u][[f]]->A[[f]]`, THM-4255 proves

```text
ker(ev_g)=(u-g),
ev_g^(-1)(f^N)=f^N+(u-g).
```

Thus `ev_g(P)` having order `N` does not put `P` coefficientwise in `f^N`.
The sharper hostile `u-g+f^N` has source order zero and specialized order
`N`; multiplying by `z` preserves the zero-`z`-constant condition. If the
manuscript first contracted along a section `u=u(f)`, its no-backflow input
would indeed be unavailable. In a coefficient module, the analogous map
`(a,b)->a+qb` kills `(q,-1)`.

## Letter 2 -- universal-torsor response

The frozen source is commit
[`4c87bcdf`](https://github.com/octonion/p-adic-zeta-density/commit/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57).
At lines 1042--1049 it places the identity **on one fixed flat torsor chart**
and leaves `u` free. The products and pole symbols retain `u`. The predecessor
source states explicitly at
[lines 1528--1546](https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/hybrid_arithmetic_holonomy.tex#L1528-L1546)
that `F_lambda` is pulled back to the chart, is independent of the torsor
point, and uses fibre variables as coefficients. A zero scalar coefficient
therefore pulls back to the zero element of the chart algebra. No inverse
to evaluation is invoked, and `u-f` remains nonzero.

Section 12 similarly digitizes in a finite-free global coefficient lattice
and asserts `Pi_b(q)=Gamma_b(q,q^ell)` before any displayed scalarization. A
full-source search finds no `u=u(f)` or `u=u(q)` map. The `u(q)` occurring in
the local expression `h=q^(-D)u(q)` is a scalar numerator, not a torsor fibre
coordinate. The actual `z=f^ell` and `z=q^ell` specializations are protected
by `deg_f<ell` and `ord_q(w_r)<=D+C<ell`, respectively; the hostiles
`z-f^ell` and `z-q^ell` violate those anti-aliasing bounds.

## Resolution

**Conceding party:** the initial chosen-section audit.

**Agreed conclusion:** the kernel theorem is correct, but its hypothesized map
is absent. The objection does not break Propositions 6.2/12.3. Proposition
6.2 should instantiate its completed chart `A_ell[[f]]`; Proposition 12.3
should supply the universal coefficient-module diagram and scalar inclusion.
If either universal statement failed and only a hidden section remained, the
hostile and conditional dependency collapse would apply immediately.

**Action taken:** THM-4255 records the exact algebra and sharp repairs;
MISTAKE-527 records the audit error; the frozen source audit keeps the
external density claims under specialist review without retracting them on
this ground.
