---
id: THM-3210
title: "Factorial exterior double-cancellation ray and two-step invisible window"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  THM-3186's exit-time amplitude satisfies deg_v E_L = L-2 over Q(n,d), so the complete
  length-three cancellation set is the rational graph
  `v=(n+3)(d-n-2)/(2d[(4n+9)d-(n+3)(4n+7)])`, and THM-3186's published hostile
  is exactly its `(n,d)=(1,5)` point.  On the ray `d=n+4`,
  `v=(n+3)/(3(n+4)(2n+5))` BOTH `E_3` and `E_4` vanish identically in `n`,
  while `Delta`, `v`, `beta_n`, `c_(n+1)=2` and `E_5` stay nonzero.  The
  visibility profile at lengths `2,3,4,5` is therefore
  visible / invisible / invisible / visible: the invisible window is two steps
  wide and visibility is NOT monotone in length.  For `n=1,2,3` that ray is the
  unique nonzero rational double-cancellation locus, and no triple
  cancellation `E_3=E_4=E_5=0` exists.
audit: >
  The exact symbolic companion checks `deg_v E_L=L-2` for four base indices and
  five lengths, derives and matches the closed-form length-three locus on nine
  base indices, confirms that THM-3186's hostile is a point of it, verifies
  `E_2=2`, `E_3=E_4=0` and the closed form of `E_5` symbolically in `n`,
  replays the double cancellation numerically together with strict positivity
  of `Delta`, `v` and `beta_n` for `n=1..40`, computes the shared resultant
  root to exclude triple cancellation, verifies rational uniqueness of the ray,
  and carries a positive control that moving `d` off the ray restores
  length-four visibility.  A separate companion imports neither SymPy nor
  primary code: custom Q[N,D,V] arithmetic reproduces the general locus and
  ray, exact Sylvester determinants reproduce all six resultants, Euclidean
  gcds exclude triple cancellation, and finite-field Frobenius tests prove the
  three quartics irreducible modulo 17, 19, and 13.  Normal/-O/stored replays
  are byte-identical and AST audits find no assert node.
source: death-star-imw-pi-n-s2-2026-08-02
depends_on:
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
related:
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
  - THM-3204-parabolic-continuant-single-gate-and-jacobi-smith-obstruction
script: 04-computation/factorial_double_cancellation_ray_thm3210.py
output: 05-knowledge/results/factorial_double_cancellation_ray_thm3210.out
script_sha256: d76de542b5904825deae4d556b150eec97311d44f7246ec7b1f394f2676157ba
output_sha256: f052824666dd572a569ca07d67aa1ec1489b370c963055bc35b5d958993916d2
independent_script: 04-computation/factorial_double_cancellation_ray_independent_audit_thm3210.py
independent_output: 05-knowledge/results/factorial_double_cancellation_ray_independent_audit_thm3210.out
independent_script_sha256: d20c2af92875cd76a0380205bbdcc864efd3207228c0a6405d4d9f7f463f4c45
independent_output_sha256: 4e88fb2df7f81dbfddc3cb86a1d52769a7a01d478acca03de4570fad1da17f71
hash_basis: LF-normalized bytes
---

# THM-3210 -- factorial exterior double-cancellation ray and two-step invisible window

**PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

THM-3186 proved that the signed exit-time amplitude

```text
E_L=sum_(j=1)^(L-1) c_(n+j) C_(j-1) prod_(h=j+1)^(L-1) u_(n+h)         (1)
```

can vanish while every local datum stays a unit, and exhibited one hostile,
`n=1, d=5, v=4/105`, at which `V_2 != 0` but `V_3 = 0`.  That witness is not
isolated and its window is not one step.

## 1. The amplitude has generic degree `L-2` in `v`

With THM-3183's factorial state
`a_i=2(i+1)(2i+1)v`, `b_i=i(i+1)Delta`, `c_i=d-i-1`, `Delta=1-4dv`,
`u_i=-b_i`, `alpha_i=a_i d`, `beta_i=b_i d`, the continuant recursion
`C_r=alpha_(n+r)C_(r-1)+d beta_(n+r)C_(r-2)` gives, in
`Q[n,d][v]`,

```text
deg_v E_L = L-2,        L>=2.                                          (2)
```

Here is an all-`L` proof rather than a finite degree fit.  The coefficient of
`v^r` in `C_r` comes only from the all-monomer term
`prod_(i=1)^r alpha_(n+i)`.  Hence the coefficient of `v^(L-2)` contributed by
exit time `j` has the form

```text
d^(L-2)(d-n-j-1) K_(L,j)(n),

K_(L,j)(n)
 =prod_(i=1)^(j-1) 2(n+i+1)(2n+2i+1)
  prod_(h=j+1)^(L-1) 4(n+h)(n+h+1).                       (2a)
```

For `n>=1`, every `K_(L,j)(n)` is positive.  Therefore the coefficient of
`d^(L-1)v^(L-2)` in `E_L` is the nonzero sum
`sum_(j=1)^(L-1)K_(L,j)(n)`.  This proves `(2)` over `Q(n,d)`; a special
numerical value of `d` may lower the degree and is not part of the generic
claim.

Thus `E_2` is a nonzero constant (the entry pivot), `E_3` is **affine** in
`v`, `E_4` is a conic, and so on.  The primary companion additionally checks
`n=1,2,3,4` and `L=2,...,6` exactly.

## 2. The complete length-three cancellation locus

Because `E_3=c_(n+1)u_(n+2)+c_(n+2)alpha_(n+1)` is affine in `v`,

```text
E_3=(n+2)[-(d-n-2)(n+3)Delta+2(2n+3)vd(d-n-3)],                        (3)
```

and its unique root is the rational graph

```text
v = (n+3)(d-n-2) / ( 2d[ (4n+9)d - (n+3)(4n+7) ] ).                    (4)
```

This is the **complete** length-three cancellation set, not a sample.  Two
degenerations are visible in `(4)`:

* `d=n+2` gives `v=0`; there `c_(n+1)=0`, so THM-3186's *entry* wall already
  fires and the length-two pivot `V_2=-c_(n+1)beta_n` vanishes.  This is not a
  signed cancellation.
* `(4n+9)d=(n+3)(4n+7)` makes the denominator vanish, so on that single ray
  `E_3` is a nonzero constant in `v`: **no length-three cancellation exists at
  all**.

THM-3186's hostile is the `(n,d)=(1,5)` point of `(4)`:
`v=(4)(2)/(2*5*(13*5-44))=8/210=4/105`.

## 3. The double-cancellation ray

**Theorem.** Put

```text
d=n+4,            v=(n+3)/(3(n+4)(2n+5)).                              (5)
```

Then, identically in `n`,

```text
Delta=(2n+3)/(3(2n+5)),      beta_n=n(n+1)(n+4)(2n+3)/(3(2n+5)),
c_(n+1)=2,                                                             (6)

E_2=2,     E_3=0,     E_4=0,                                           (7)

E_5=-4(n+2)(n+3)^2(n+4)(2n+3)(10n^3+101n^2+336n+366)/(27(2n+5)^2).     (8)
```

For every `n>=1` each factor in `(6)` and `(8)` is strictly positive, so
`v != 0`, `Delta != 0` and `2vDelta` is invertible: the ray satisfies
THM-3183's standing hypothesis, and `beta_n != 0` excludes the entry wall.
Since `V_L=-beta_n E_L`, the visibility profile is

```text
length      2          3            4            5
V_L      visible   invisible    invisible     visible.                 (9)
```

**Two consequences.**

1. **The invisible window is two steps wide.** THM-3186 established one
   cancelling length; `(7)` gives two consecutive ones, uniformly in the base
   index `n`.
2. **Visibility is not monotone in length.** A source can be visible, become
   invisible for two steps, and become visible again.  No "first visible
   depth" or "once visible, always visible" argument survives, and a search
   that certifies visibility by testing one or two lengths past the graph
   distance is unsound.

The `n=1` member of `(5)` is exactly THM-3186's published hostile.  That
witness was therefore already a double cancellation; only `V_3` was tested
there.

## 4. Uniqueness and the absence of a triple point

**FINITE-EXACT (`n=1,2,3`).** Let `Res_v` denote the resultant in `v`.

```text
Res_v(E_3,E_4) = const * d^2 * (d-(n+4)) * Q_n(d),                    (10)
```

with `Q_n` an irreducible quartic over `Q` (`Q_1=13d^4-96d^3+288d^2-672d+1008`,
`Q_2=17d^4-160d^3+600d^2-1800d+3600`,
`Q_3=7d^4-80d^3+360d^2-1320d+3300`).  Hence `d=n+4` is the **unique nonzero
rational** double-cancellation ray for those base indices.  Moreover

```text
gcd( Res_v(E_3,E_4), Res_v(E_3,E_5) ) = const * d^2,                  (11)
```

whose only root `d=0` is excluded, so **no triple cancellation
`E_3=E_4=E_5=0` occurs** for `n=1,2,3`.  Whether a three-step invisible window
exists for some larger `n`, or over an extension field, is **OPEN**.

The independent audit proves the three quartics irreducible after reduction
modulo `17,19,13`, respectively, using the degree-four Frobenius criterion;
it also reconstructs `(10)--(11)` from hand-written Sylvester determinants.

## 5. What this changes and what it does not

It sharpens THM-3186's projection no-go: the failure of the support graph and
the local Smith profiles to determine projected visibility is not a
codimension-large accident but a rational ray present at every base index, and
the failure persists for two consecutive lengths.  Any future identification of
THM-3183's offset-six PRS pivots with an oriented continuant product must
therefore carry an amplitude that survives `(5)`, not merely a graph-distance
or first-nonzero-length statistic.

It does **not** prove the empirical `floor(s/2)` Euclidean-depth staircase, any
fixed-offset closure, or any `NC(2)`, `GMC(2)`, `JC(2)` or `LRC(14)`
consequence, and `(10)`--`(11)` are finite-exact at `n=1,2,3` only.

Run

```text
python 04-computation/factorial_double_cancellation_ray_thm3210.py
python -O 04-computation/factorial_double_cancellation_ray_thm3210.py
python 04-computation/factorial_double_cancellation_ray_independent_audit_thm3210.py
python -O 04-computation/factorial_double_cancellation_ray_independent_audit_thm3210.py
```

and compare LF-normalized bytes with the two declared outputs.  Exact symbolic,
sparse-polynomial, rational, and finite-field arithmetic only; no floating
point, random sampling, imported executable, or assertion-sensitive test.

**QED** for sections 1--3; section 4 is finite-exact in its stated universe.
