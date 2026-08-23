---
id: THM-3726
title: "Automorphic Cohn constant SL2 orbit Broughton classification"
status: >
  PROVED + VERIFIED-EXACT.  Let R=[a,b;c,d] be any constant SL2 matrix.
  A constant combination of either row of M0R with the other can be closed
  only on the two-sheet locus (2a-d)^2=1, with an explicit typed linear
  compatibility condition.  Every such closed row is exactly L+L^2S in
  Jacobian-one linear source coordinates, and the complementary curl is 6S
  or -6S.  Hence no final polynomial left shear completes it.  This
  classifies constant exposed parameters over the whole constant right-SL2
  orbit; nonconstant exposed parameters and nonconstant right matrices remain
  outside the theorem.
source: root / 2026-08-22
audit: >
  SELF-AUDITED; independent hostile audit requested.  The exact companion
  checks both raw closure systems and their discriminant resultant, universal
  parametrizations of both sheets and orientations, both Broughton potentials
  and debts, all 180 integer SL2 matrices in the [-4,4] box, and finite
  Broughton cokernels through degree 12.  Normal and optimized runs
  byte-match the frozen transcript.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
related:
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3723-automorphic-cohn-c3-two-right-resonance
  - THM-3725-automorphic-cohn-hyperbolic-two-right-resonance
script: 04-computation/jc2_automorphic_cohn_constant_sl2_orbit_thm3726.py
output: 05-knowledge/results/jc2_automorphic_cohn_constant_sl2_orbit_thm3726.out
script_sha256: af6430f5c767b620451f197bc42686a6b20f663b7698702d9d55dabfaef1daa1
output_sha256: 2d17226da5525d039960232693cb923b2fe19e12ddb462a4dcb0c34ba3ced882
semantic_sha256: a9efb6a1de61faaad9d40407dcc091f89647b392029001c7fda25d2acf2905c7
hash_basis: raw LF bytes
---

# THM-3726 -- the whole constant right orbit has one Broughton closure locus

**PROVED + VERIFIED-EXACT.**  The `C_3` resonance of THM-3723 and the
hyperbolic/Pell resonance of THM-3725 look group-theoretically different.
They are nevertheless two elementary-factor charts on one geometric object:
the constant `SL_2` orbit's two-sheet row-closure locus.

Work over a characteristic-zero field `k`.  Put

```text
M_0=[4T^2,2XT-1;1+2XT,X^2],
R=[a,b;c,d] in SL_2(k),
N=M_0R,                                                det N=1.    (1)
```

Write the rows of `N` as `alpha,beta`, and let `h in k` be constant.

The lower combination `beta+h alpha` is closed exactly when

```text
ch=d-a,                         (4a-d)h=b.              (2)
```

The upper combination `alpha+h beta` is closed exactly when

```text
(a-d)h=-c,                      bh=4a-d.                (3)
```

Either system forces

```text
epsilon=2a-d in {+1,-1}.                                (4)
```

Conversely, `(4)` together with the corresponding typed system `(2)` or
`(3)` gives a closed row.  The square condition alone is not sufficient at a
degenerate zero-coefficient boundary; the linear system is the exact claim.

Every survivor has a polynomial potential

```text
Q=L+L^2S,                         J_(X,T)(L,S)=1,       (5)
```

and its complementary row has curl `6S` in the lower orientation and `-6S`
in the upper orientation.  Therefore, for arbitrary `f,g in k[X,T]`, neither

```text
E_+(f)E_-(h)M_0R,                 E_-(g)E_+(h)M_0R     (6)
```

is a Jacobian matrix whenever the exposed parameter `h` is constant.

## 1. The closure discriminant

Direct differentiation gives

```text
curl(beta+h alpha)
 =2{[a+ch-d]X+[(4a-d)h-b]T},                           (7)

curl(alpha+h beta)
 =2{[(a-d)h+c]X+[4a-d-bh]T}.                          (8)
```

Thus `(2)--(3)` are necessary and sufficient.  Compatibility of either
two-equation system is measured by the same expression:

```text
cb-(4a-d)(d-a)
 =(a-d)(4a-d)+bc
 =(2a-d)^2-(ad-bc).                                   (9)
```

Since `ad-bc=1`, every closure lies on `(4)`.  Equation `(9)` also explains
why the typed systems must be retained: if both coefficients multiplying
`h` vanish, scalar compatibility still has to be checked.

## 2. Universal lower-sheet parametrization

Suppose `(2)` holds and put `epsilon=2a-d`.  The complete parametrization is

```text
a=epsilon+ch,
d=2a-epsilon,
b=(2a+epsilon)h,                    epsilon^2=1.        (10)
```

Indeed `(10)` gives `ad-bc=epsilon^2` and is equivalent to `(2)` on
`SL_2`.  Define

```text
L=epsilon(X+2hT),
S=[cX+(2a+epsilon)T]/3.                                (11)
```

Using `(10)`, one checks

```text
J_(X,T)(L,S)=1,
d(L+L^2S)=beta+h alpha.                                (12)
```

The complementary curl is

```text
curl(alpha)=2[cX+(4a-d)T]=6S.                          (13)
```

A final upper left shear would have to solve

```text
J(f,L+L^2S)=6S,                                        (14)
```

which is impossible by THM-3716.

## 3. Universal upper-sheet parametrization

For `(3)`, use the equivalent parametrization

```text
d=2a-epsilon,
c=(a-epsilon)h,
bh=2a+epsilon,                      epsilon^2=1.        (15)
```

Set

```text
L=epsilon(hX+2T),
S=[(a-epsilon)X+bT]/3.                                (16)
```

Then

```text
J_(X,T)(L,S)=1,
d(L+L^2S)=alpha+h beta,                                (17)
```

and the complementary debt is

```text
curl(beta)=2[(a-d)X-bT]=-6S.                           (18)
```

The nonzero scalar form of THM-3716 again excludes the final shear.

## 4. Group-theoretic synthesis and honest boundary

The discriminant coordinate

```text
epsilon=2a-d                                             (19)
```

is the missing invariant behind the earlier factor charts:

- the THM-3723 `C_3` matrix has `(a,d)=(0,1)` and `epsilon=-1`;
- the THM-3725 trace-four Pell matrix has `(a,d)=(1,3)` and
  `epsilon=-1`;
- the identity and triangular single-shear boundaries occupy the
  `epsilon=+1` sheet.

Thus elliptic, hyperbolic, and triangular right words can expose different
row formulas while producing the same Broughton component.  The group type
of a factorization is not the Hamiltonian-cokernel type of its closed row.

The theorem classifies **constant** exposed combinations for every constant
`R`, hence constant right words of arbitrary elementary length after they are
multiplied out.  It does not prove that a nonconstant polynomial `h` cannot
close a row for a general constant `R`; THM-3721, THM-3723, and THM-3725 do
prove that stronger fact in their respective cells.  Nor does this theorem
classify nonconstant right matrices.  Any counterexample remaining in the
constant right orbit must therefore use a genuinely nonconstant exposed left
parameter outside those completed cells.

Reproduce the exact checks with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_constant_sl2_orbit_thm3726.py
python3 -B -O 04-computation/jc2_automorphic_cohn_constant_sl2_orbit_thm3726.py
```

The exact integer grid contains 180 `SL_2` matrices and 74 closed orientation
rows; it guards all zero-entry boundaries.  Equations `(7)--(18)` prove the
classification uniformly.  **QED.**
