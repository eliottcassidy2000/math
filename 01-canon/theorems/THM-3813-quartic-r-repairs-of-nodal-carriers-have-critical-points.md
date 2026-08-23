---
id: THM-3813
title: "Quartic R-repairs of nodal carriers have critical points"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  On the c=1 cubic pseudo-plane, every canonical
  nodal carrier A=e^2-z/3+r sum_(i=0)^4 b_i e^i has a critical point.  In the
  genuine quartic case, the degree-twenty-two residual resultant cannot be
  supported on e*g=0: after passing to the intrinsic weight-zero parameters,
  the first seven logarithmic-remainder coefficients have two incompatible
  Groebner normal forms.  This closes every pure-r profile of degree at most
  four, but not degree at least five or mixed corrections.  The candidate is
  not a proved dependency until independent audit promotion.
source: jc_sparse_direct_search / quartic logarithmic-divisibility lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The deterministic companion derives the
  universal Hamiltonian reduction and resultant, degree and leading
  coefficient, the localized Euclidean quotient, all seven invariant
  remainder equations, the 28-element exact Groebner basis, both incompatible
  normal forms, repeated-root and sparse hostile controls, projective
  finite-root gate, and source reconstruction.  Normal and optimized runs
  byte-match the frozen transcript.  Independent hostile audit is pending.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3806-binomial-cubic-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3807-trinomial-and-full-cubic-r-repairs-have-critical-points
related:
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
  - THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points
script: 04-computation/jc2_cubic_pseudoplane_quartic_r_repair_thm3813.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_quartic_r_repair_thm3813.out
script_sha256: 0f4058d30385a532470a2f2a2b776c9cabda26bceba2653ce985e41a69141235
output_sha256: 6efa9ae0837a42018f57bcf8fca235e75f76f37611166d94280f56f311f8dd1f
semantic_sha256: 32a2b2c1373a41c5f17e307d0b5ba6e898384fa7832815d0b86648ace8592704
hash_basis: raw LF bytes
---

# THM-3813 -- every quartic r-repair remains critical

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
characteristic zero.  On the `c=1` member of the THM-3785 cubic pseudo-plane
put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary `b_0,...,b_4 in k`, the regular function

```text
g(e)=b_0+b_1e+b_2e^2+b_3e^3+b_4e^4,
A=e^2-z/3+r g(e)                                      (2)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

When `b_4=0`, the complete degree-at-most-three result is the union of
THM-3799, THM-3805, THM-3806, and THM-3807.  It remains to prove the genuine
quartic case

```text
b_4 != 0.                                               (3)
```

## 1. Universal residual and the boundary criterion

Put

```text
u=re,                         K=1+2u,
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (4)
```

The Hamiltonian components computed from `(1)` are

```text
{A,r}=r^2-9z^2(2e+rg'),
{A,z}=3gr^2-3(1+2re)(2e+rg'),
{A,e}=9gz^2-(1+2re).                                  (5)
```

After `r=u/e`, one has `P/e^2={A,z}/3`; `Q` is the compatibility equation
for

```text
z^2=K/(9g),                     z^3=u(1+u)/e.          (6)
```

Keeping independent symbols `G,D` in place of `g,g'`, exact elimination in
`u` gives

```text
Res_u(P,Q)=G^3e^4 H_univ(e,G,D).                       (7)
```

After substituting the quartic profile, write the residual factor as `H(e)`.
Its degree and leading coefficient are

```text
deg H=22,                       LC(H)=19131876 b_4^5.  (8)
```

Thus `(3)` makes `H` a nonzero polynomial of fixed degree twenty-two.

Suppose, for contradiction, that every root of `H` lies on the forbidden
boundary `V(e g)`.  Factoring over `k`, with arbitrary multiplicities,

```text
H=mu product_alpha (e-alpha)^(n_alpha)                 (9)
```

implies

```text
e g H'/H=sum_alpha n_alpha e g/(e-alpha) in k[e].     (10)
```

In particular,

```text
H divides e g H'.                                     (11)
```

This necessary implication is valid without assuming that `g` is squarefree:
if `e g` has order `m>=1` and `H` has order `n>=1` at a boundary point, then
`e g H'` has order at least `m+n-1>=n` there.

## 2. The intrinsic weight-zero coefficient packet

Divide `e g H'` by `H` over `Q(b_0,...,b_4)[e]`.  The quotient is

```text
q = 22b_4e^4 +(53b_3+4)e^3/3
    +(44b_2/3+(-4b_3^2-4b_3+8)/(9b_4))e^2
    +(13b_1-4b_2b_3/(3b_4)
      +(8b_3^3-24b_3+16)/(27b_4^2))e
    +38b_0/3+4b_1(1-b_3)/(3b_4)-8b_2^2/(9b_4)
    +(28b_2b_3^2-20b_2b_3-8b_2)/(27b_4^2)
    +(-16b_3^4+16b_3^3+48b_3^2-80b_3+32)/(81b_4^3).
                                                               (12)
```

The only parameter denominator divides `81b_4^3`, which is invertible under
`(3)`.  Define the cleared remainder

```text
R=81b_4^3(e g H'-qH)=sum_(j=0)^21 r_j e^j.             (13)
```

If `(11)` held, every `r_j` would vanish.  The large raw coefficient packet
has a much smaller intrinsic description.  Set

```text
C_0=b_0b_4^3,       C_1=b_1b_4^2,       C_2=b_2b_4,
C_3=b_3,            T=b_4^7.                            (14)
```

These are precisely the weight-zero parameters for the recurring mod-seven
scaling, and `(3)` gives `T!=0`.  Formally substitute

```text
b_0=C_0/s^3, b_1=C_1/s^2, b_2=C_2/s, b_3=C_3, b_4=s,
T=s^7.                                                   (15)
```

For `21>=j>=15`, the pullback of `r_j` is a nonzero rational scalar times

```text
s^(j-17) F_j(C_0,C_1,C_2,C_3,T),                       (16)
```

where `F_j` is the signed integer primitive part.  The powers in descending
order are exactly

```text
4,3,2,1,0,-1,-2.                                      (17)
```

Because `s=b_4!=0`, no equation is added or lost in `(15)--(17)`.  In
particular, `(11)` would force

```text
F_21=F_20=...=F_15=0.                                 (18)
```

This passage is evaluation in five explicit invariants, not a claim that a
quotient parametrization is birational.  It therefore retains sparse profiles,
repeated roots of `g`, and every specialization with `b_4!=0`.

## 3. Seven top equations are already inconsistent

Work in the exact polynomial ring

```text
Q[C_0,C_1,C_2,C_3,T]                                  (19)
```

with graded reverse lexicographic order on the displayed variable ordering,
and let

```text
I=(F_21,F_20,F_19,F_18,F_17).                         (20)
```

Exact Buchberger reduction gives a 28-element Groebner basis for `I`.  The
normal forms of the next two primitive equations are

```text
NF_I(F_16) = -T(2688C_2-15308821)/4,
NF_I(F_15) = 27T(378784C_2-2308816611)/8.              (21)
```

The companion reconstructs every `F_j` directly from `(7),(12)--(16)`,
checks that each differs from its raw coefficient by only a nonzero rational
scalar, computes the exact basis, and freezes both the seven-polynomial packet
and the basis by independent hashes.  Thus `(21)` is an exact characteristic-
zero ideal-membership calculation, not a numerical solve or bounded census.

At any putative common zero of `(18)`, the first five equations make each
polynomial congruent to its normal form.  Since `T!=0`, the last two equations
would force simultaneously

```text
C_2=15308821/2688,
C_2=2308816611/378784.                                (22)
```

But cross multiplication gives

```text
15308821*378784-2308816611*2688=-407362596704 !=0.    (23)
```

This contradicts `(18)`.  Therefore `(11)` is impossible, and `H` has a root
`eta` with

```text
eta g(eta) != 0.                                      (24)
```

## 4. The off-boundary root reconstructs a critical point

At `e=eta`, equations `(7)--(8)` make the `u`-resultant vanish.  The leading
coefficient of the quartic `Q` is

```text
LC_u(Q)=-729g(eta)^3 !=0.                             (25)
```

Hence the homogenized pair has no common root at infinity, even if the degree
of `P` drops.  There is a common finite root `u_0`.  Directly,

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16,                          (26)
```

so `u_0`, `1+u_0`, and `K_0=1+2u_0` are nonzero.  Define

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                     (27)
```

No square or cube root is chosen.  Exact reduction gives

```text
z_0^2-K_0/(9g(eta)) = -Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0 = u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0 = -Q/(eta^2K_0^2),
{A,z}|_0 = 3P/eta^2.                                  (28)
```

Thus `(r_0,z_0,eta)` lies on `Y` and the last two Hamiltonian components
vanish.  The surface Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                    (29)
```

and `K_0!=0` kill `{A,r}`.  Since `Y` is smooth and symplectic, this is an
actual critical point of `A`, excluding every regular Darboux mate.

## 5. Root-space interpretation and exact scope

The weight-zero compression also explains a misleading local resonance.  At
a nonzero root `alpha` of `g`, with `d=g'(alpha)`, the universal residual has

```text
H(alpha)=d alpha[-729alpha(d-4alpha^2)^3-2].           (30)
```

Put `w=d-4alpha^2`.  The nontrivial arm law is
`alpha=-2/(729w^3)`.  The apparent high-contact factor
`531441w^7+16` would then give

```text
d w^6=w^7+16/531441=0,                                (31)
```

so `d=0`: it is a repeated-root seam, not an exceptional simple arm.  The
global seven-coefficient calculation above includes that seam automatically;
no local squarefreeness assumption enters the proof.

Together with THM-3799, THM-3805, THM-3806, and THM-3807, this candidate
closes every pure `r g(e)` carrier with `deg g<=4`.  It does **not** address
degree at least five, mixed `z^2h(e)+r g(e)` corrections, another pseudo-plane
arm profile, or the existence of arbitrary planar Keller maps.  The exact
companion has 429 active gates; normal and optimized executions byte-match
the frozen transcript.  **QED, subject to independent hostile audit.**
