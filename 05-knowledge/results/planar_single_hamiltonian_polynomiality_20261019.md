# Polynomiality of a single source Hamiltonian flow

2026-10-19. **PROVED FROM CITED INPUT.** This is an application of
Gwoździewicz's fixed-line theorem. No historical novelty claim is made.

## Statement

Let K be a characteristic-zero field, let

    H=c+t²R(x,t),   c∈K,   R∈K[x,t],

and use the bracket `{A,B}=A_xB_t−A_tB_x`. Set `δ_H={−,H}`. For any
nonzero `λ∈K`, let

    Φ_λ=(exp(λδ_H)x, exp(λδ_H)t)

be the original source t-adic flow, with coordinates in `K[x][[t]]`.
Then the following are equivalent:

1. Both coordinates of Φ_λ belong to `K[x,t]`.
2. The second coordinate of Φ_λ equals t.
3. `H_x=0`, equivalently `H∈K[t]`.

In these cases the flow is the polynomial shear

    Φ_λ=(x+λH′(t), t).

In particular, if `H_x≠0`, no nonzero scalar-time source flow is a
polynomial map. This conclusion has no Newton-support or second-boundary
assumption.

## Proof

First,

    δ_H x=H_t∈tK[x,t],
    δ_H t=−H_x∈t²K[x,t].

Consequently δ_H raises t-adic order by at least one. Its exponential is
a continuous automorphism of `K[x][[t]]`, with inverse `exp(−λδ_H)`. Write
its coordinate images as `(X,T)`. They satisfy

    X≡x mod t,   T≡t mod t²,   {X,T}=1.                 (1)

The last identity follows because a Hamiltonian derivation is a derivation
of the Poisson bracket, so its exponential preserves that bracket. All
sums here are coefficientwise finite in the t-adic topology.

Suppose X and T are polynomials. Extend scalars to an algebraic closure
`K̄`. Equation (1) makes `(X,T)` a polynomial Keller map over `K̄` fixing
the line `t=0` pointwise. Gwoździewicz's fixed-line lemma therefore makes
it a polynomial automorphism. This is precisely Lemma 2.3 of
[Janusz Gwoździewicz, *Injectivity on one line*](https://arxiv.org/abs/alg-geom/9305008);
the more general injective-on-one-line statement is Theorem 1.1. Both are
stated over algebraically closed characteristic-zero fields.

As the image of the irreducible coordinate t under a polynomial ring
automorphism, T is irreducible in `K̄[x,t]`. But (1) says that t divides T.
Hence `T=αt` for some `α∈K̄*`; its first source jet in (1) gives `α=1`.
Thus polynomiality implies `T=t`. One may also recover the shear form
directly: `{X,t}=X_x=1`, whence `X=x+f(t)` with `f∈K[t]` and `f(0)=0`.
No descent theorem for automorphisms is needed.

Next suppose `H_x≠0`, and let `m=ord_t(H_x)≥2`. Since δ_H raises t-order,

    ord_t(δ_H^n t)≥m+n−1   for every n≥1.

Therefore

    T−t=−λH_x mod t^(m+1),                             (2)

which is nonzero: its degree-m coefficient is −λ times the nonzero
degree-m coefficient of H_x. Thus `T=t` is impossible. This proves
`(1)⇒(2)⇒(3)` without assuming that an arbitrary formal exponential has
no periods; the required absence of a period is supplied by its first
nonzero source jet.

Finally, if `H_x=0`, characteristic zero gives `H∈K[t]`, and
`δ_H=H′(t)∂_x`. Hence

    exp(λδ_H)x=x+λH′(t),   exp(λδ_H)t=t,

which proves `(3)⇒(1)` and the formula. ∎

## Finite commuting words with one common Hamiltonian

Assume additionally `H_x≠0`. Let `f_1,…,f_r∈K[s]` and
`λ_1,…,λ_r∈K`, and put

    F(s)=Σ_i λ_i f_i(s),
    Ψ=∏_i exp(λ_i δ_(f_i(H))).

Each `f_i(H)` has the form `f_i(c)+t²S_i(x,t)`, so all factors exist in
the same source t-adic completion. Since

    δ_(f_i(H))=f_i′(H)δ_H,   δ_H(H)=0,

these derivations commute. Thus the finite word is exactly

    Ψ=exp(δ_(F(H))).

The preceding theorem applies to `F(H)=F(c)+t²S(x,t)`. Moreover,

    ∂_x F(H)=F′(H)H_x.

Because `H_x≠0`, H is nonconstant; substitution `K[s]→K[x,t]`,
`s↦H`, is injective. Characteristic zero therefore gives

    Ψ is polynomial
      ⇔ F′(H)H_x=0
      ⇔ F is constant
      ⇔ Ψ is the identity.

This permits cancellation between factors and includes zero scalar times.
If instead `H∈K[t]`, every such word is a polynomial shear, explicitly
`(x+(F(H(t)))′,t)`.

## Scope and verification

- The result concerns both coordinates being polynomial in the original
  affine source coordinates. It does not classify rational, algebraic, or
  transcendental coordinate images.
- It does not assert that each coordinate is individually nonpolynomial
  when `H_x≠0`. For example, `H=xt²` has
  `Φ_λ=(x(1+λt)²,t/(1+λt))`: its first coordinate is polynomial.
- The commuting-word corollary requires a single common H. Arbitrary
  different-invariant compositions are not classified here.
- This is not a proof of the planar Jacobian conjecture; the fixed-line
  hypothesis is essential in the cited input.

The proof was checked symbolically and logically, including arbitrary
characteristic-zero base fields and the scalar-time completion issue.
The exact primary-source statements were read on 2026-10-19 in the
[author's two-page preprint](https://arxiv.org/pdf/alg-geom/9305008).
No numerical census or finite-range computation is used.
