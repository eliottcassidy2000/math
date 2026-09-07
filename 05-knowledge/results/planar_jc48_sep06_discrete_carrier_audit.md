# Independent audit: one-way polynomial carrier descent

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.** September 6, 2026.
Auditor: `three_ray_geometry`; producer: `orthogonal_returns`, with root's
exceptional-fibre strengthening. No producer file was edited.

Audited artifacts:

- [Proof](planar_jc48_sep06_discrete_carrier.md).
- [Source](../../04-computation/planar_jc48_sep06_discrete_carrier.py).
- [Frozen output](planar_jc48_sep06_discrete_carrier.out).

## 1. Main classification and quantifiers

I independently reconstructed the carrier map for
`p=t(1+x^2t)`, `y=xtp`, `Delta=p^3-y^2`. The identities
`Delta=t p^2`, `t=Delta/p^2`, and `x=yp/Delta` give exactly the asserted
image `{(0,0)} union {p Delta!=0}`. The inverse reconstructs every point
of the stated open set in both directions. The only positive-dimensional
fibre is the disjoint union `L={t=0}` and `M={1+x^2t=0}`. Its two
components are respectively `A1` and `G_m`; all other fibres are empty
or singletons.

One-way inclusion `Phi^*(k[p,y]) subset k[p,y]` indeed gives a polynomial
map `f` with `pi Phi=f pi`. The proof never assumes that this induced
map is invertible or surjective. Since `Phi(L)` is a closed curve inside
the fibre of `f(0)`, that image value must be zero. This forces
`Phi(L union M) subset L union M`.

Each component image is a whole component, since it is closed,
irreducible, and one-dimensional. The two images are different because
`Phi` is injective. The unit groups of `k[z]` and `k[z,z^-1]` forbid
swapping them. Thus each component is preserved, and the two principal
prime ideals give

```text
Phi^*(t)=ct,
Phi^*(1+x^2t)=d(1+x^2t),          c,d in k*.
```

Restriction to `t=0` gives `d=1`, hence `c Phi^*(x)^2=x^2` in the
polynomial domain. Factorization yields `Phi^*(x)=b x`, `c=b^-2`.
The resulting unique `lambda=b^-1` is exactly the claimed scaling,
with source Jacobian `lambda` and full carrier equality. The converse
is immediate. No hidden finiteness, dimension, or component-permutation
case is omitted.

The hypotheses remain polynomial **automorphism** and the **entire**
carrier. This is not an argument that an arbitrary Keller endomorphism
is invertible. Requiring Jacobian one removes every nonidentity scaling.

## 2. Source-form, rational, and formal checks

For preservation of every `G_H=-u/2+H(p,y)`, the three tests `H=0,p,y`
give carrier inclusion by subtraction. Applying the classification then
fixes `u` and supplies the asserted iff. Fixing just one `G_H` is
explicitly weaker, as required by the rational hostile.

The two-way Poisson consistency route is sound. With constant carrier
Jacobian `j` and source Jacobian `kappa`, covariance gives
`j Delta Phi^*(p)=kappa p Phi^*(Delta)`. Unique zero/pole divisors force
the irreducibles `p` and `Delta` to be preserved individually up to
units. The polynomial carrier automorphism therefore has the form
`p->alpha p`, `y->beta y+h(p)`. Preserving `Delta` yields `h=0` and
`beta^2=alpha^3`; covariance yields `kappa=alpha^2/beta`, agreeing with
the source scaling. This route is correctly restricted to the two-way
subcase.

I checked the Hamiltonian formulas for `S=p^2 Delta`. They show
`L(A) subset A` and `L(-u/2) in A`. These inclusions, rather than the
finite iterate controls, prove that `exp(epsilon L)` is an automorphism
of the completed source and takes every affine source form into
`-u/2+A[[epsilon]]`. Characteristic zero permits the exponential
coefficients; the factor `epsilon` makes the series well defined in
the completion. It has inverse `exp(-epsilon L)` and is symplectic.
No polynomial scalar-time specialization is inferred. The stated
nonconstant-LND stopping point is also consistent with factorial
closure: a nonconstant `S` in `k+Delta A` has its invariant nonconstant
part divisible by both `t` and `1+x^2t`, which would force `t` and then
`x` into an LND kernel, making the derivation zero.

For the separate rational action, direct substitution gives its inverse,
additive law, Jacobian one, and exact invariance of `u`. The carrier
expressions have the declared denominators `Delta^2` and `Delta^3`.
Reduction of the numerators modulo `Delta` gives `w^2p^6` and
`-w^3p^9`; at every nonzero scalar `w` these do not vanish generically
on the irreducible cusp divisor. Thus the claimed pole orders are
exact, not merely denominator upper bounds. The original rational
`x` coordinate also has an uncancelled pole. Fixing `G_0` alone does
not imply whole-carrier descent.

Finally, the polynomial symplectic translation sends the collapsed
`M`-fibre to carrier values `(2z+1)/z^4`, with the two displayed distinct
values `3` and `5/16`. This is a decisive exact failure of descent on an
actual whole fibre, not only a formal or fraction-field obstruction.

## 3. Source audit and independent replay

The standalone program uses always-active `check` calls, and its declared
universe is a fixed symbolic identity set. I inspected simultaneous
substitution, both rational inverse directions, the indeterminate
scaling, rational group law and residues, and ordinary-source versus
carrier-ring Hamiltonian iteration through order three. The literal
Poisson correction and second symplectic exponential coefficient are
checked separately. No finite-order control is used to infer the
all-order formal statement or the complete automorphism classification.

Independent runs were performed with

```sh
python3 -B 04-computation/planar_jc48_sep06_discrete_carrier.py
python3 -B -O 04-computation/planar_jc48_sep06_discrete_carrier.py
```

Both outputs equal the saved output byte for byte. All **67 gates** pass.
The independently recomputed raw pins agree:

```text
source SHA256 f748d1b7f1ce09b2a1ab7f72bf9fb397b2bb4bc3ae1a24951041ad1b90ca8dcd
output SHA256 528bb8a3a7ca013d17e80b32cea8cfd0b7c31a4ebc3ddd692a151f881808884e
semantic SHA256 250f154a1a40184013f186e918a3f2026fdf4f4a35b1f0a7e9b0c5b8ff88c7ec
```

No mathematical or source correction is requested. The producer may
promote its independent-audit status while keeping the source and output
frozen. The precise carrier and automorphism hypotheses remain essential.
