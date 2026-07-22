---
id: THM-2101
title: "Two additive orbit-residue proofs of one-variable DvdK"
status: >
  PROVED ON PAPER; FORMALIZATION PARTIAL. For a genuinely two-sided complex
  Laurent polynomial, all positive-power constant terms cannot vanish. Two
  independent additive proofs use the same barycentric residue observable.
  One analytically continues the small-root identity through transitive
  monodromy; the other specializes once at a small parameter transcendental
  over the coefficient field and applies a Galois root-packet lemma. Both
  contradict the zero full-root Lagrange sum and bypass the small-root product,
  Hensel factorization, logarithms, and Wiener--Hopf. Check A, irreducibility,
  additive incidence, and the full-root identity are kernel-checked. The
  analytic monodromy route, the transcendental contour/lift wrapper, and the
  final DvdK1 interface are not yet formalized.
source: codex-2026-07-22-GMC2-additive-orbit-residue
related:
  - THM-1550
  - THM-1765
  - THM-2022
  - THM-2067
  - HYP-8946
  - HYP-8960
  - MISTAKE-243
formalization:
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2PhiIrreducible.lean
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2LaurentShiftCheckA.lean
script: 04-computation/dvdk_additive_orbit_residue_codex_20260722.py
output: 05-knowledge/results/dvdk_additive_orbit_residue_codex_20260722.out
script_sha256: d2f67a19df74b6095b0ab5f8a08eed3ddb874e3422609f89c8b60e4c05222f29
output_sha256: d2f4ebf4607ceeec3936c5ff124be0326b916e6adae6069ed31dbe99e6e2ce30
hash_basis: repository blobs with LF line endings
---

# THM-2101 -- two additive orbit-residue proofs of the strict DvdK theorem

## 1. Statement

Let

```text
f(u)=u^(-M)R(u),
R(u)=r_0+r_1u+...+r_d u^d,
d=M+N,                    M,N>=1,          r_0 r_d!=0.   (1)
```

Then

```text
CT(f^m)!=0 for some integer m>=1.                       (2)
```

This is the one-variable DvdK statement needed as the seed in THM-2022. It is
bare existence; THM-2111 independently gives the effective bound

```text
m<=binom(M+N,min(M,N)).                                 (3)
```

If `CT(f)=r_M` is nonzero, (2) is immediate. We therefore assume, toward a
contradiction,

```text
CT(f^m)=0                         for every m>=1.        (4)
```

## 2. The additive small-cluster observable

Put

```text
Phi(u,t)=u^M-tR(u),
w_t(alpha)=alpha^(M-1)/Phi_u(alpha,t)                   (5)
```

at a simple root `alpha` of `Phi`. Choose a small circle `|u|=delta` and then
`|t|` small enough that

```text
|tR(u)|<|u|^M                       on |u|=delta.        (6)
```

Rouche's theorem gives exactly `M` roots of `Phi` inside the circle. Call
this set `S_0(t)`. Away from the discriminant, the contour integral and the
geometric series are both valid and give

```text
G(t)=sum_(m>=0) CT(f^m)t^m
    =(2*pi*i)^(-1) integral_(|u|=delta) du/[u(1-tf(u))]
    =(2*pi*i)^(-1) integral_(|u|=delta) u^(M-1)du/Phi(u,t)
    =sum_(alpha in S_0(t)) w_t(alpha).                  (7)
```

The sign and normalization in (7) are literal: the residue of
`u^(M-1)/Phi(u,t)` at a simple root is exactly the quotient in (5). Assumption
(4) turns (7) into the analytic germ identity

```text
sum_(alpha in S_0(t)) w_t(alpha)=1.                    (8)
```

This is the key change of viewpoint. The observable is the **sum** of the
small-root residues, before taking a product or applying a logarithm.

## 3. Monodromy transports the identity

Consider the rational map

```text
h:P^1 -> P^1,                    h(u)=u^M/R(u).          (9)
```

Because `R(0)!=0`, the numerator and denominator are coprime. Since
`deg R=M+N>M`, the map has degree `d=M+N`. Remove its finitely many branch
values together with `0,infinity`, and choose a sufficiently small finite
regular basepoint `t_b` in the remaining punctured disk. The restricted source
is a sphere minus finitely many fibers and remains connected, so the resulting
degree-`d` cover has a finite monodromy group `G_mon` acting transitively on
the fiber

```text
Omega={alpha:Phi(alpha,t_b)=0}.                         (10)
```

Equivalently, `Phi` is irreducible over `C(t)`. As a polynomial in `t` over
`C[u]`, it is linear with coprime coefficients `u^M` and `-R(u)`, hence is
irreducible in `C[u,t]`. As a polynomial in `u` over `C[t]`, it is primitive:
the coefficients `-t r_0` and `1-t r_M` are coprime. Gauss then gives
irreducibility in `C(t)[u]`.

Analytically continue (8) around a loop representing `sigma in G_mon`.
Permanence of a zero analytic germ and equivariance of (5) give

```text
sum_(alpha in S_0(t_b)) w_(t_b)(sigma alpha)=1
                                      for every sigma in G_mon.             (11)
```

This is the local-to-global bridge: one does not identify an unordered
analytic subset with a canonical subset of an abstract splitting field.
Instead, continuation produces exactly the orbit of that subset, and the
identity travels with it. A translated subset need not remain inside the
original small circle; only the seed subset is selected by that local place.

## 4. The full-root sum vanishes

At `t_b`, sum the finite residues of the rational differential

```text
u^(M-1) du/Phi(u,t_b).                                  (12)
```

Its coefficient at infinity is

```text
u^(M-1)/Phi(u,t_b)=O(u^(-N-1)).                         (13)
```

Since `N>=1`, the residue at infinity is zero. Therefore

```text
sum_(beta in Omega) w_(t_b)(beta)=0.                   (14)
```

This is also the elementary Lagrange sum for an exponent `M-1<=d-2`; no
monic normalization is needed because (12) uses the actual derivative of the
actual polynomial.

Now sum (11) over the finite group. Transitivity says that for fixed
`alpha in Omega`, the map `sigma |-> sigma alpha` reaches every root with
fiber size `|Stab(alpha)|`. Hence

```text
|G_mon|
 =sum_(alpha in S_0) sum_(sigma in G_mon) w(sigma alpha)
 =sum_(alpha in S_0) |Stab_G_mon(alpha)|
                       sum_(beta in Omega) w(beta)
 =0.                                                    (15)
```

This is impossible over `C`. Assumption (4) is false, proving (2). QED.

## 5. Exact scope and hostile boundary

Both signs in the Laurent support are load-bearing. If `N=0`, then
`M-1=d-1`, (13) has an `u^(-1)` term, and the full-root sum need not vanish.
For example `f=u^(-1)` has every positive constant term zero, while its small
cluster is the full root set and its residue sum is one.

No simplicity problem is hidden: the basepoint is chosen outside the finite
branch-value set. No choice of a global small-root subset is claimed. The
proof uses the local cluster only to seed its full monodromy orbit.

The wrapper to the repository's `DvdK1` interface is elementary. For an
injectively charged finite support, a charge zero gives the conclusion at
`m=1`. Otherwise `0` in the convex hull supplies exact extreme charges
`q_min<0<q_max`; set `M=-q_min`, `N=q_max`, and
`R(u)=sum_i c_i u^(q_i+M)`. Injectivity and nonzero support coefficients give
`r_0 r_(M+N)!=0`, while `f=u^(-M)R`. The Check A identity already formalized
in `GMC2LaurentShiftCheckA` identifies these Laurent constant terms with the
universal relation used by `DvdK1`. Only this concrete wrapper and the
analytic monodromy package remain to be assembled in Lean.

The result removes the small-root **product** `Pi`, THM-1550, Hensel lifting,
Wiener--Hopf, and the logarithmic local/global bridge from the paper proof of
the bare seed. The Lean formalization remains partial: Check A, additive orbit
incidence, and the full-root identity are kernel-checked after MISTAKE-243;
analytic continuation and the final `DvdK1` wrapper are still absent.

## 6. Abel integration explains the old bottleneck

For comparison, let `Pi(t)` be the product of the `M` small roots and put
`c=(-1)^(M+1)r_0`, so `Pi(t)/(ct)` tends to one. Implicit differentiation at
`Phi(alpha,t)=0` gives

```text
t alpha'(t)/alpha(t)=alpha^(M-1)/Phi_u(alpha,t).        (16)
```

Summing (16) over the small cluster and using (7) yields

```text
G(t)=t Pi'(t)/Pi(t),
log(Pi(t)/(c t))=sum_(m>=1) CT(f^m)t^m/m.              (17)
```

Thus the product route applies the local Abel operator

```text
A(G-1)=integral_0^t (G(s)-1) ds/s.                     (18)
```

The `1/m` and the product/Hensel bottleneck are created by integrating the
observable. Equation (15) stays before that integration and closes by additive
orbit incidence. This is a convergent-germ, hence coefficientwise formal,
identity; it is unrelated to Abel boundary summation or Dini--Bertrand tests.

## 7. Assumption challenge and Tournament Analysis

The challenged assumption is that Galois symmetry should act on roots or
subset products as the primary vertices. Here the faithful vertices are the
monodromy translates of one **root subset**, carrying additive residue weight.
The pairwise observable is subset intersection size, and transitivity makes
its incidence regular. There is no canonical tournament orientation on these
translates: the faithful finite object is a block-incidence hypergraph, in
fact a uniform `1`-design under the transitive action. Score histograms,
cycles, SCCs, edge flips, and Hamiltonian-path counts would discard the only
load-bearing data, namely uniform root incidence and the conserved subset sum.

The carrier

```text
(Omega, G_mon orbit of S_0, weights w, full sum zero)   (19)
```

preserves the constant-term predicate through (7)--(11). It forgets which
roots were locally small; the seed subset and analytic-continuation sidecar
restore exactly that information and no product data.

## 8. Independent transcendental-specialization proof

For a second proof write `Lambda=f` and again assume total vanishing. Equation
labels restart inside this proof.

```text
CT(Lambda^m)=0 for every m>=1                           (2)
```

is impossible. Equivalently, some positive power has nonzero constant term.
This is the strict two-sided one-variable DvdK statement needed after the
lowest-face reduction in THM-2022.

The proof has two independent parts: an algebraic root-packet lemma and one
fixed complex contour computation. No analytic continuation of a root subset
in the parameter is used.

### 8.1. Algebraic barycentric root-packet lemma

Let `F` be a characteristic-zero field, let `p in F[X]` be irreducible of
degree `d`, and let `L/F` be a splitting field. For `0<=k<=d-2` and a subset
`S` of the roots of `p` in `L`, put

```text
b_k(S)=sum_(alpha in S) alpha^k/p'(alpha).              (3)
```

> **Root-packet lemma.** If `b_k(S)` belongs to `F`, then `b_k(S)=0`.

Indeed, characteristic zero makes `p` separable and its Galois group `G`
acts transitively on the roots. The derivative weight is equivariant:

```text
sigma(alpha^k/p'(alpha))
  =(sigma alpha)^k/p'(sigma alpha).                    (4)
```

If (3) is in `F`, every translate `sigma S` has the same sum `b_k(S)`. Sum
over `sigma in G`. Transitivity makes the incidence multiplicity of each root
independent of that root, whereas Lagrange interpolation gives

```text
sum_(p(alpha)=0) alpha^k/p'(alpha)=0,     0<=k<=d-2.   (5)
```

Thus `|G| b_k(S)=0`, and characteristic zero proves the lemma. Formula (5)
does not require `p` to be monic: writing `p=a n` with `n` the monic nodal
polynomial merely multiplies every derivative denominator by `a`.

The same proof works after embedding an abstract splitting field into any
algebraically closed extension: pull the concrete subset back along the root
equivalence first. One must not pretend that an analytically selected subset
is itself Galois-stable.

### 8.2. One small transcendental specialization

Shift (1) exactly by its least exponent:

```text
R(z)=z^M Lambda(z) in C[z],
d=deg R=M+N>M,
Phi_t(z)=z^M-tR(z).                                    (6)
```

Let `k` be the subfield of `C` generated over `Q` by the finitely many
coefficients of `R`. It is countable. Put

```text
B=sum_(j=0)^d |r_j|>0.                                 (7)
```

The elements algebraic over `k` form a countable subset of `C`, so every
punctured disk contains a number transcendental over `k`. Choose `tau` with

```text
0<|tau|<1/B,             tau transcendental over k.    (8)
```

Evaluation `t -> tau` is injective on `K=k(t)`. By the irreducibility theorem
kernel-checked in `GMC2PhiIrreducible`, `Phi_t` is irreducible over `K` because
`M>=1` and `R(0)!=0`. Let `L/K` be a splitting field. Since `C` is
algebraically closed, the injective evaluation of `K` extends to an embedding

```text
iota:L -> C.                                           (9)
```

It identifies the abstract roots with the roots of

```text
Phi_tau(z)=z^M-tau R(z).                               (10)
```

No specialized collision occurs: irreducibility in characteristic zero makes
`Phi_t` separable, and injectivity of evaluation preserves its nonzero
discriminant. Thus all derivative weights below are defined.

### 8.3. The inside-root sum is one

On the unit circle, (7)--(8) give

```text
|tau Lambda(z)|=|tau z^(-M)R(z)|<=|tau|B<1.           (11)
```

Hence (10) has no root on the circle. Let `S_tau` be the finite subset of its
roots inside the unit disk. The rational function

```text
z^(M-1)/Phi_tau(z)=1/[z(1-tau Lambda(z))]              (12)
```

is regular at zero because `Phi_tau(0)=-tau R(0)!=0`. Partial fractions, or
equivalently the simple-pole contour formula, give

```text
sum_(alpha in S_tau) alpha^(M-1)/Phi_tau'(alpha)
 = (2*pi*i)^(-1) integral_(|z|=1) dz/[z(1-tau Lambda(z))]. (13)
```

The strict bound (11) makes the geometric series uniformly convergent on the
circle, so termwise integration is legitimate:

```text
right side of (13)
 =sum_(m>=0) tau^m CT(Lambda^m)=1.                     (14)
```

The last equality uses (2), with the `m=0` term equal to one. Rouché would
also show `|S_tau|=M`, but neither the proof nor its formalization needs that
cardinality.

Pull `S_tau` back through (9) to a subset `S` of the roots in `L`. Derivatives
commute with the embedding, so injectivity of (9) turns (14) into

```text
sum_(alpha in S) alpha^(M-1)/Phi_t'(alpha)=1 in K.     (15)
```

But `M-1<=d-2`, since `N>=1`. Applying the root-packet lemma to (15) says its
left side must be zero. This contradiction proves the theorem. QED.

### 8.4. Exact scope and hostile controls

- A neutral charge is dispatched first: if the constant coefficient of
  `Lambda` is nonzero, (2) already fails at `m=1`.
- The exact shift `M=-min(support)` is load-bearing. A larger shift inserts a
  zero root and destroys the irreducible/simple-derivative setup.
- Both endpoint coefficients must be nonzero and `N>=1`. For the one-sided
  hostile example `Lambda=z^(-1)`, all positive constant terms vanish, but
  the only root-subset sum and the full-root sum both equal one; here
  `deg Phi=M`, so (5) is inapplicable at exponent `M-1`.
- The proof handles arbitrary complex coefficients. Transcendence is required
  only of the chosen parameter `tau` over their finitely generated coefficient
  field; it imposes no arithmetic restriction on those coefficients.

### 8.5. Formalization boundary

The repaired `GMC2LaurentShiftCheckA` kernel-checks the exact Laurent shift,
constant-term coefficient identity, additive incidence contradiction, and
full-root Lagrange sum without `sorryAx`. `GMC2PhiIrreducible` kernel-checks
the full variable-swap/Gauss irreducibility over `F(t)`. Both are root-imported.

What remains for a Lean theorem is assembly, not a missing mathematical
bridge: construct the countable coefficient field and transcendental
specialization; lift the splitting field; prove the partial-fraction contour
identity and uniform geometric-series interchange; normalize the nonmonic
derivative; and connect the exact-shift theorem to `GMC2DvdKInterface.DvdK1`.
Mathlib already supplies the relevant Lagrange interpolation, circle-integral,
dominated-convergence, `RatFunc.liftAlgHom`, and splitting-field lift pieces.

THM-1550/HYP-8960's Henselian small-root-product route remains independently
interesting, and its power-series local-ring instance is now kernel-checked,
but it is not a dependency of this additive proof.
