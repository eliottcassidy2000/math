---
id: THM-2760
title: "Exact-prefix even Faber flux gcd and smooth-boundary exclusion"
status: >
  PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  every m=4k-2, exact-prefix specialization with nonzero square-root residue
  omega gives gcd(Phi_m/q,Psi_m)=r^(k-1) after q=r omega^3 and localization
  at omega.  The residual factors are coprime because one explicit carrier
  is a Schur--Szego composition of two simple-negative Jacobi transforms.
  Thus no nondegenerate exact-prefix boundary point with q omega nonzero can
  satisfy both top fluxes.  On the projective root-zero chart omega=0,q!=0,
  the two fluxes vanish together exactly for k=2 mod 3, equivalently m=6 mod
  12, and the third response is then nonzero.  The omitted affine corner q=0
  is outside that projective chart.  This is a boundary theorem; it does not
  classify the full h=0
  intersection, derive another response chart, or close another degree.
source: jc-even-zero-flux-next-2026-07-28
audit: >
  all-degree-split-hostile-audit-2026-07-28 independently replayed the
  specialized recurrence and root-zero lane through k=30, checked the
  load-bearing multiplicity, hyperbolicity, and simple-B-root statements
  against the primary Kostov--Shapiro TeX source, and verified that the
  proof uses the strict multiplicity threshold correctly: PASS.
external_input: >
  Kostov--Shapiro, On Schur--Szego composition of polynomials,
  arXiv:math/0605377, Propositions 4--5 and Theorem 6(ii)
related:
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2752-all-even-zero-first-flux-response-regularization-closure
  - THM-2755-all-even-zero-flux-componentwise-global-regularity-closure
  - THM-2759-split-degree-six-componentwise-monicization-closure
script: 04-computation/jc2_exact_prefix_even_faber_flux_gcd_thm2760.py
output: 05-knowledge/results/jc2_exact_prefix_even_faber_flux_gcd_thm2760.out
script_sha256: dd2fc28c09f4917428572b508c9f73acfb8d46b8c64053a7f584ccec4f9ed6e3
output_sha256: 5f09b3cb9e89d45a8299698b4e78779ac5d623426937f71ad26f690fe4797152
hash_basis: LF-normalized bytes
---

# THM-2760 -- the exact-prefix even flux gcd is uniform

**PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The degree-22 resultant `-76608 omega^6` is one specialization of an
all-degree fact.  Once the exact-prefix square-root residue is nonzero, the
two top fluxes have only their forced power of the normalized cubic residue
in common.  The proof exposes the residual polynomial as a strict
Schur--Szego composition; this also explains the mod-three exceptional lane
at square-root residue zero.

## 1. Statement and scope

Let `k>=1`, put

```text
m=4k-2,                         alpha=m/4=k-1/2,
```

and define the Faber coefficients by

```text
(1+2dz^2+qz^3+(d^2-s)z^4)^alpha=sum_(n>=0)c_n z^n,

Phi_m=4c_(m+1),                 Psi_m=4c_(m+2),
R_m=4c_(m+3)+2d c_(m+1).                              (1)
```

At an exact-prefix boundary lift the two residue identities are

```text
d=-omega^2,                     s=q omega.             (2)
```

If `omega!=0`, write

```text
q=r omega^3,                    s=r omega^4.            (3)
```

Then in `C[omega,omega^(-1),r]`, up to nonzero rational and `omega` units,

```text
gcd(Phi_m/q,Psi_m)=r^(k-1).                            (4)
```

The odd index in `Phi_m` makes it universally divisible by `q`, so the first
equation in `(4)` is polynomial.  Consequently no point with `q omega!=0`
can satisfy both exact-prefix top flux equations.  On the projective chart
`omega=0,q!=0`, where `(2)` gives `d=s=0`, the two top fluxes vanish together
exactly when

```text
k=2 mod 3,                       equivalently m=6 mod 12. (5)
```

On this exceptional lane `R_m!=0`.  Thus `(5)` is a genuine response-active
root-zero candidate, not a common zero of all three responses.  If one also
sets `q=0`, all positive-degree coefficients of `(1)` vanish; that affine
corner is not a point of the projective residue chart and is not classified
by `(5)`.

The theorem does not assert that every top point has `q omega!=0`, classify
the other components of the top intersection, or show that the root-zero
candidate in `(5)` carries a physical source map.  It supplies the uniform
exact-prefix exclusion used at the five nondegenerate points of THM-2752 and
isolates the degree `6 mod 12` boundary that a higher-degree argument must
retain.

## 2. The forced `r^k` factor

After the change of variable `x=omega z`, equations `(2)`--`(3)` turn the
quartic inside `(1)` into

```text
B_r(x)=1-2x^2+r x^3+(1-r)x^4
      =(1-x^2)^2+r x^3(1-x).                          (6)
```

Define

```text
P_k(r)=r^(-k)[x^(4k-1)]B_r(x)^(k-1/2),
Q_k(r)=r^(-k)[x^(4k  )]B_r(x)^(k-1/2).                (7)
```

These are polynomials.  Indeed, expand by perturbation order `j`:

```text
B_r(x)^(k-1/2)
 =sum_(j>=0) binom(k-1/2,j) r^j x^(3j)
       (1-x)^(2k-1-j)(1+x)^(2k-1-2j).                 (8)
```

For `j<k`, both final exponents are nonnegative and the summand in `(8)` has
degree exactly `4k-2=m`.  It contributes to neither coefficient in `(7)`.
Hence both top fluxes are divisible by `r^k`.

For `ell>=0`, direct coefficient extraction gives, with coefficients outside
their displayed finite ranges omitted,

```text
p_ell=[r^ell]P_k
 =binom(k-1/2,k+ell)(-1)^(k-1-3ell)2^(k-1-3ell)
      binom(k-1-ell,2ell),                            (9)

q_ell=[r^ell]Q_k
 =binom(k-1/2,k+ell)(-1)^(k-3ell)2^(k-3ell-1)
      (2binom(k-1-ell,2ell-1)+binom(k-1-ell,2ell)).  (10)
```

Here `0<=ell<=floor((k-1)/3)` in `(9)` and
`0<=ell<=floor(k/3)` in `(10)`.  Elementary binomial cancellation yields

```text
(k-3ell)q_ell=-(k+ell)p_ell.                          (11)
```

At the terminal index `ell=k/3`, when it exists, both sides of `(11)` are
zero because `p_ell=0`.

## 3. One squarefreeness carrier controls the gcd

Let

```text
theta=r d/dr,                         S_k=P_k-3Q_k.    (12)
```

Equation `(11)` is the polynomial differential identity

```text
(k-3theta)Q_k=-(k+theta)P_k.                          (13)
```

Solving `(12)`--`(13)` gives

```text
P_k=(k-3theta)S_k/(4k),
Q_k=-(k+theta)S_k/(4k).                               (14)
```

The constant terms of `P_k` and `Q_k` are nonzero and opposite, so zero is
not a root of `S_k`.  At any `rho!=0`, equations `(14)` show

```text
P_k(rho)=Q_k(rho)=0
iff S_k(rho)=S_k'(rho)=0.                             (15)
```

Thus `P_k,Q_k` are coprime exactly when `S_k` is squarefree.

Put `n=floor(k/3)` and normalize `H_k=S_k/S_k(0)`.  Equations `(9)`--`(12)`
simplify to the positive factorial formula

```text
H_k(r)=sum_(ell=0)^n h_(k,ell)r^ell,

h_(k,ell)=
 k k! (k-1-ell)! /
 (32^ell ell! (k-3ell)! (k+ell)!).                   (16)
```

## 4. The mod-three Jacobi factorization

For `k=3n+j`, define `(a,b)` by

```text
j=0: (a,b)=(n-1/3,n-2/3),
j=1: (a,b)=(n+1/3,n-1/3),
j=2: (a,b)=(n+2/3,n+1/3).                            (17)
```

Write falling and rising factorials as `u_down_ell` and `u_up_ell`.  The
factorial identity

```text
k!/((k-3ell)!ell!)
 =binom(n,ell)27^ell a_down_ell b_down_ell            (18)
```

turns `(16)` into

```text
h_(k,ell)=binom(n,ell)(27/32)^ell
  (a_down_ell/(k-1)_down_ell)
  (b_down_ell/(k+1)_up_ell).                          (19)
```

Define two degree-`n` polynomials in binomial basis,

```text
A_k(x)=sum binom(n,ell) a_down_ell/(k-1)_down_ell x^ell,
B_k(x)=sum binom(n,ell) b_down_ell/(k+1)_up_ell x^ell. (20)
```

If `*` denotes degree-`n` Schur--Szego composition,

```text
(sum binom(n,ell)u_ell x^ell)
 *(sum binom(n,ell)v_ell x^ell)
 =sum binom(n,ell)u_ell v_ell x^ell,
```

then `(19)` says exactly

```text
H_k(r)=(A_k*B_k)((27/32)r).                           (21)
```

Both factors have simple negative roots.  For `A_k`, reversal gives

```text
x^n A_k(1/x)
 =unit * 2F1(-n,k-n;a-n+1;-x).                       (22)
```

The right side is a Jacobi polynomial with parameters

```text
alpha_A=a-n,                  beta_A=k-a-n-1.          (23)
```

For `B_k`, first write

```text
B_k(x)=2F1(-n,-b;k+1;x)
 =(1-x)^n 2F1(-n,k+1+b;k+1;x/(x-1)).                 (24)
```

The final hypergeometric polynomial has Jacobi parameters

```text
alpha_B=k,                    beta_B=b-n.              (25)
```

Every parameter in `(23)` and `(25)` is greater than `-1`.  Explicitly,
`a-n` is `-1/3,1/3,2/3`, `b-n` is `-2/3,-1/3,1/3`, and
`k-a-n-1` is respectively `n-2/3,n-1/3,n+1/3`.  The standard Jacobi zero
lemma follows directly from Rodrigues' formula: orthogonality against every
polynomial of degree less than `n` under the positive weight
`(1-x)^alpha(1+x)^beta` forces `n` sign changes in `(-1,1)`.  Hence all `n`
zeros are distinct and interior.  The transformations `(22)` and `(24)`
therefore put all zeros of `A_k,B_k` on the negative real axis, distinctly.

The remaining input is the precise Schur--Szego multiplicity theorem of
Kostov--Shapiro, not a generic preservation slogan.  For real hyperbolic
degree-`n` polynomials `U,V` with every root of `V` negative, their
composition is hyperbolic.  If roots of multiplicities `mu,nu` satisfy
`mu+nu>n`, their forced product root has multiplicity `mu+nu-n`; every other
root (a `B`-root in their terminology) is simple.  This is Propositions 4--5
and Theorem 6(ii) of
[Kostov--Shapiro, arXiv:math/0605377](https://arxiv.org/abs/math/0605377),
also published in *C. R. Math.* 343 (2006), 81--86.

Apply it to `(20)`.  When `n>=2`, all input multiplicities are one, so no
pair has sum greater than `n`; every output root is a simple `B`-root.  The
coefficients are positive, so the real nonzero roots are negative.  The case
`n=1` is immediate and `n=0` is constant.  Therefore `H_k`, hence `S_k`, has
only simple negative roots.  Equations `(15)` prove

```text
gcd(P_k,Q_k)=1.                                       (26)
```

## 5. Returning to the homogeneous fluxes

Coefficient scaling under `x=omega z` and `(7)` gives

```text
Phi_m   =4 omega^(4k-1) r^k P_k(r),
Psi_m   =4 omega^(4k  ) r^k Q_k(r).                   (27)
```

Since `q=r omega^3`, equations `(26)`--`(27)` prove `(4)`.  Before
localizing `omega`, if one retains the dimensional variable `q`, the raw gcd
has one additional factor `omega` exactly when `k=2 mod 3`; this factor is a
unit on the chart of `(4)`.  Confusing the dimensional and localized gcds is
the sharp hostile to the clean statement.

At a boundary lift with `q omega!=0`, both prefactors in `(27)` are nonzero.
The two top flux equations would force a common root of `P_k,Q_k`, contrary
to `(26)`.  This is the uniform nondegenerate exact-prefix exclusion.

## 6. The root-zero chart and the `6 mod 12` lane

When `omega=0` and `q!=0`, equations `(2)` give `d=s=0`, and the series in
`(1)` is

```text
(1+qz^3)^(k-1/2).                                    (28)
```

Its coefficient of `z^N` is zero exactly when `3` does not divide `N`.
The two flux indices are `4k-1` and `4k`, congruent modulo three to `k-1`
and `k`.  They are both nonmultiples of three exactly for `k=2 mod 3`, which
proves `(5)`.  On that lane the response index `4k+1` is divisible by three,
and

```text
R_m=4 binom(k-1/2,(4k+1)/3) q^((4k+1)/3) !=0.         (29)
```

The half-integral upper argument and `q!=0` make `(29)` nonzero.
Thus degrees `m=6,18,30,...` retain one response-active root-zero boundary
candidate.  The theorem neither removes that candidate nor turns its
nonzero response into a componentwise pole-count conclusion.

## 7. Exact reproduction and failure boundary

The independent audit checked the load-bearing external input directly in
the primary Kostov--Shapiro TeX source.  Proposition `th:mult` gives a
product-root only when the input multiplicity sum is strictly greater than
the degree and explicitly excludes equality; Proposition `prop:comp`
preserves hyperbolicity for the relevant negative-root factor; and Theorem
`th:ordpart(ii)` makes every remaining `B`-root simple.  The two degree-`n`
Jacobi factors used here have simple negative roots, so for `n>=2` they
create no `A`-roots and all output roots are simple; `n=1` is immediate.
The audit also independently replayed the specialized recurrence and the
root-zero congruence through `k=30`.  Those finite checks corroborate, but do
not replace, the all-degree Schur--Szego proof.

Run

```bash
python3 04-computation/jc2_exact_prefix_even_faber_flux_gcd_thm2760.py
python3 -O 04-computation/jc2_exact_prefix_even_faber_flux_gcd_thm2760.py
```

Both executions byte-match

```text
05-knowledge/results/jc2_exact_prefix_even_faber_flux_gcd_thm2760.out.
```

The companion checks `k=1..18` by two independent constructions: the Faber
coefficient recurrence and formulas `(9)`--`(10)`.  It verifies `(11)`--
`(21)`, exact negative simple-root counts for both Jacobi factors and their
carrier, gcd `(26)`, the localized equality `(4)`, the extra-`omega` hostile,
and the root-zero classification `(5)`--`(29)`.  These finite checks are
hostile controls, not the source of the all-`k` quantifier.

```text
script_sha256 = dd2fc28c09f4917428572b508c9f73acfb8d46b8c64053a7f584ccec4f9ed6e3
output_sha256 = 5f09b3cb9e89d45a8299698b4e78779ac5d623426937f71ad26f690fe4797152
hash_basis    = LF-normalized bytes
```

This theorem does not classify the full top boundary, prove nonvanishing of
the third response at every other top point, control lower-row banks, derive
the split chart in another degree, or prove `JC(2)` or `DC(2)`.  QED.
