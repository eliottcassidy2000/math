---
id: THM-3870
title: "Vertical-axis polynomial profiles force a projective companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  If the
  depressed-cubic discriminant contains the vertical axis, its polynomial
  profile differs at a finite A-adic order from the unique nonpolynomial
  binomial zero profile.  The first mismatch forces a residual component
  through finite-base C-infinity whose affine normalization omits another,
  distinct point above A-infinity.
source: root / vertical first-mismatch completion of THM-3866, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  The 862-gate audit independently derives
  the axis classification and positive-sign formal comparator, contracts
  completed divisibility back to k[A,C], and checks N=0, T=0, every degree
  regime, arbitrary A-dependent leading coefficients, shared-q_N
  cancellation, multiplicities, factor selection, dominance and the two
  distinct normalization punctures.  Its normal, optimized and frozen
  transcripts byte-match.  The separate 260-gate primary companion checks
  the universal response and a 30-row N=0,...,5 regime grid.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3866-all-polynomial-graph-branches-force-projective-companion
  - THM-3871-first-nongraph-triangular-parabola-companion
script: 04-computation/jc2_vertical_axis_first_mismatch_companion_thm3870.py
output: 05-knowledge/results/jc2_vertical_axis_first_mismatch_companion_thm3870.out
script_sha256: d4bbbee11cee32ecc93e057d36a9e0e3f5e959d2ec0fe43ea43cfeb81c29a381
output_sha256: 70f16eec8f8b25e956d73fec6df761e96d7d9ee07f81363678d9f1f3131194b3
semantic_sha256: 32e524ea2a1b86e16afe4160574b6e186b65b39b60a1c03a6b71ec51b943451d
independent_audit_script: 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py
independent_audit_output: 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
independent_audit_script_sha256: c5f6d3f8350d4d1c41cb910c9ce891cf140d822172bf4448cc99707202a89012
independent_audit_output_sha256: 0a9c9ea793027f79e578271b5f53f396c705161fffe90c6f74cd11a8f426410a
independent_audit_semantic_sha256: 48aa2b48290aa53595bf9bbea0f0a7015fde2ea13515543cf62b6a941b857079
hash_basis: raw LF bytes
---

# THM-3870 -- the vertical axis forces a projective companion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  For
`b=b(A,C) in k[A,C]`, put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.
```

If `A` divides `Delta_b`, then after removing the vertical factor the
reduced branch packet contains an irreducible component `Xi` distinct from
`A=0` whose affine normalization omits at least two distinct points of its
smooth projective normalization.  One lies above
`A=0,C=infinity`, and another lies above `A=infinity`.

Thus the vertical-axis lane is closed inside this depressed-cubic carrier.
The theorem says nothing about a non-graph, non-axis `A1` component, another
carrier, a Keller atlas, or `JC(2)`.

## 1. Classification, comparator, and first mismatch

For

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b,
```

specialization gives

```text
Delta_b(0,C)=9C^2-54b(0,C).
```

Thus `A | Delta_b` if and only if `b(0,C)=C^2/6`, equivalently, uniquely,

```text
b=C^2/6+A Q(A,C),                 Q in k[A,C].                 (A1)
```

There is no hidden sufficiency direction here: substitution in (A1) gives
the explicit polynomial branch formula

```text
Delta_b/A = -C^3-54Q-54ACQ-(3/4)AC^4-9A^2C^2Q-27A^3Q^2.       (A2)
```

Put

```text
P=1+(2/3)AC,                 u=1+AC+A^2b.
```

The exact identity `A^2 Delta_b=27(P^3-u^2)` shows that a zero profile in
`k[C][[A]]` with the prescribed constant term `u=1 mod A` must have
`u=P^(3/2)`.  The other sign has constant term `-1` and is impossible.
Consequently

```text
b_*=(P^(3/2)-1-AC)/A^2
   =C^2/6+A Q_*,
Q_*=sum_(j>=0) q_j(C)A^j,
q_j=gamma_j C^(j+3),
gamma_j=binom(3/2,j+3)(2/3)^(j+3) !=0.                        (A3)
```

Uniqueness can also be seen directly in (A2): after lower coefficients have
been fixed, the coefficient `q_j(C)` first occurs with the unit coefficient
`-54`.  The independent checker derives this recursion through eight levels
without assuming that `q_j` is a monomial, then recovers (A3).

Since every `q_j` is nonzero, `Q_*` has infinitely many nonzero
`A`-coefficients.  A polynomial `Q` therefore has a unique least finite
coefficient mismatch `N`, and

```text
Q=Q_<(N)+A^N T,
t=T(0,C) != q_N.                                                (A4)
```

This includes `N=0` and `T=0`.  Divisibility of `Q-Q_<(N)` by `A^N` is just
the equality of all earlier coefficients, so `T` is polynomial and unique.

## 2. Completed calculation and polynomial contraction

Let

```text
b_N=C^2/6+A Q_<(N),       u_N=1+AC+A^2b_N.
```

The exact response to any increment `x` is

```text
Delta_(b_N+x)-Delta_(b_N)=-54u_N x-27A^2x^2.                  (A5)
```

In `k[C][[A]]`, write

```text
b_*=b_N+A^(N+1)H_N,            H_N(0,C)=q_N(C).
```

Equations `Delta_(b_*)=0` and (A5) show that `Delta_(b_N)` belongs to
`A^(N+1)k[C][[A]]`.  Since it lies in `k[A,C]`, coefficientwise contraction
gives

```text
k[A,C] intersection A^(N+1)k[C][[A]]=A^(N+1)k[A,C],
Delta_(b_N)=A^(N+1)R_N,               R_N in k[A,C].          (A6)
```

After division and specialization,

```text
R_N(0,C)=54q_N(C).                                             (A7)
```

No analytic convergence or illegal division in the completion occurs.

For a uniform statement set `rho_N=1/6` when `N=0`, and
`rho_N=gamma_(N-1)` when `N>=1`.  The top term of `b_N` is then
`rho_N A^N C^(N+2)`.  The square term in `Delta` is the unique term of top
`C`-degree, giving

```text
deg_C R_N=2N+4,       lc_C(R_N)=-27rho_N^2 A^(N+1),
deg_C u_N=N+2,        lc_C(u_N)=rho_N A^(N+2).                 (A8)
```

This includes the hostile boundary

```text
R_0=-C^3-(3/4)AC^4.
```

## 3. Actual mismatch and every degree regime

Substituting `b=b_N+A^(N+1)T` in (A5) yields

```text
Delta_b=A^(N+1)S,
S=R_N-54u_NT-27A^(N+3)T^2,
S(0,C)=54(q_N-t) !=0.                                        (A9)
```

In particular, the residual has no factor `A`.  If `T=0`, equations
(A7)--(A8) immediately give generic degree `2N+4`, special degree `N+3`,
and a strict drop.

Suppose `T!=0`, put `d=deg_C T`, and let `tau(A)=lc_C(T)`.  The three total
degrees in (A9) are `2N+4`, at most `N+2+d`, and `2d`.

- If `d<N+2`, the `R_N` term uniquely wins and `D=2N+4`.
- If `d>N+2`, the quadratic term uniquely wins and `D=2d`.
- If `d=N+2`, the common leading coefficient is exactly

  ```text
  -27A^(N+1)(rho_N+A tau)^2.                                  (A10)
  ```

  Its bracket has nonzero constant term `rho_N`, so it is not the zero
  polynomial.  It may vanish at isolated values of `A`; that does not change
  the degree over `k[A]`.

The special degree is

```text
m=deg_C(q_N-t) <= max(N+3,d).                                 (A11)
```

Equations (A10)--(A11) give `D>m` in every regime.  This remains true when
`t` shares the leading term, several leading terms, or all but one term of
`q_N`.  In particular, at `d=N+3` one only needs `m<=d`, while `D=2d`.
If the highest `C`-coefficient of `T` is divisible by `A`, then
`deg_C t<d`, which only strengthens (A11).

The audit checker covers `N=0,...,5`, `T=0`, degrees below, at, and above
`N+2`, a resonant square with an isolated base root, and the hostile families

```text
T=q_N+1,
T=q_N+C^(N+1)+1,
T=q_N+1+A C^(N+5).
```

Thus it explicitly exercises both shared-`q_N` cancellation and an
`A`-dependent generic leading term invisible on the special fibre.

## 4. Factor selection with multiplicity

Homogenize `S` in `C` to its actual degree `D`, using `[C:Z]`.  If `S_0^h`
is the degree-`m` homogenization of the nonzero `S(0,C)`, then

```text
mathcal S(0;C,Z)=S_0^h(C,Z) Z^(D-m),
S_0^h(1,0) !=0.                                                (A12)
```

Hence the projective residual meets `P_0=(0,[1:0])`.  Factor in the UFD

```text
S=c product_i G_i^(e_i).
```

Because `C`-degrees add and each factor is homogenized at its own actual
degree, homogenization respects this product with every multiplicity.
Equation (A12) therefore selects at least one distinct reduced factor `G_i`
whose closure reaches `P_0`.

The selected factor is not `A`: equation (A9) says `S(0,C)` is nonzero.
Nor can it be a `C`-independent vertical factor through `A=0`, since such an
irreducible factor is associated to `A`.  It is not a component supported on
`Z=0`, because an affine factor homogenized at its actual `C`-degree is not
divisible by `Z`.  Thus the selected reduced curve is distinct from the axis
and dominates the `A`-line.

Multiplicity causes no loss: it changes only how much of the positive
`Z`-order in (A12) is allocated to the selected factor.  The checker factors
representative actual residuals with multiplicity.  It also uses the hostile
synthetic packet

```text
(AC+1)^2(C-A),
```

whose generic/special `C`-degrees are `3/1`; the reduced factor `AC+1`
reaches `P_0` with multiplicity two and is visibly not vertical.

## 5. Two distinct missing normalization points

Let `Xi` be a selected reduced factor, close it in
`P1_A x P1_C`, and normalize the closure.  A point `p_0` lies above `P_0`.
The pullback of the rational coordinate `C/Z` has a pole at `p_0`, so this
point is outside the affine normalization.

Dominance makes the induced projective morphism to `P1_A` nonconstant, hence
surjective.  Choose `p_infinity` above `A=infinity`.  The coordinate `A` has
a pole there, so it too is outside the affine normalization.  Their images
in `P1_A` are respectively zero and infinity; consequently
`p_0 != p_infinity`.  This proves two distinct missing normalization points,
not merely two boundary incidences that might coalesce after normalization.

The factor `AC+1` is an exact hostile model: its normalization is `P1`, and
its affine part misses precisely the displayed endpoints
`(A=0,C=infinity)` and `(A=infinity,C=0)`.

## 6. Exact and independent controls

Run from the repository root:

```bash
python3 -B 04-computation/jc2_vertical_axis_first_mismatch_companion_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_companion_thm3870.out
python3 -B -O 04-computation/jc2_vertical_axis_first_mismatch_companion_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_companion_thm3870.out
python3 -B 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
python3 -B -O 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
```

The 260-gate primary companion checks the axis classification, formal
coefficients, universal response, `N=0,...,5`, `T=0`, 30 mismatch rows and
the resonant square.  The structurally separate 862-gate hostile audit adds
`A`-dependent leading coefficients, cancellation with `q_N`, actual-degree
factorization with multiplicity, a nonreduced synthetic packet, dominance and
distinct projective endpoints.  Normal and optimized executions of both
companions byte-match their frozen LF outputs.  **QED.**
