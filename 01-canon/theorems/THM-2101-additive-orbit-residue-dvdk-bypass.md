---
id: THM-2101
title: "An additive orbit-residue proof bypasses the small-root product in one-variable DvdK"
status: >
  PROVED on paper. For a genuinely two-sided one-variable Laurent polynomial,
  total constant-term vanishing makes the sum of the barycentric residue
  weights on the small-root cluster equal to one. Analytic continuation sends
  this identity to every monodromy translate of that cluster. The rational map
  X^M/R(X) has transitive degree-(M+N) monodromy, while the sum of the weights
  over every root is zero by the residue at infinity. Uniform orbit incidence
  then says a positive group order is zero, a contradiction. This proves the
  bare one-variable DvdK seed without the small-root product, THM-1550,
  Hensel, logarithms, or Wiener--Hopf. The repaired algebraic incidence core is
  kernel-checked; the analytic monodromy bridge and final DvdK1 wrapper are not
  yet formalized in Lean.
source: codex-2026-07-22-GMC2-additive-orbit-residue
related:
  - THM-1550
  - THM-1765
  - THM-2022
  - THM-2067
  - HYP-8946
  - MISTAKE-243
formalization:
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2LaurentShiftCheckA.lean
---

# THM-2101 -- the additive orbit-residue bypass

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

Equivalently, `Phi` is irreducible over `C(t)`: viewed as a polynomial linear
in `t`, its two coefficients `u^M` and `-R(u)` are coprime, and Gauss applies.

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
 =M |Stab(alpha)| sum_(beta in Omega) w(beta)
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
its incidence regular. Orienting two translates by the first moved root gives
a tie Hamiltonian path, but score histograms, cycles, SCCs, edge flips, and
path counts are irrelevant: only the uniform root incidence and the conserved
subset sum survive the quotient.

The carrier

```text
(Omega, G_mon orbit of S_0, weights w, full sum zero)   (19)
```

preserves the constant-term predicate through (7)--(11). It forgets which
roots were locally small; the seed subset and analytic-continuation sidecar
restore exactly that information and no product data.
