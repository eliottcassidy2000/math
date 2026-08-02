---
id: THM-3125
title: "Monomial-ray first-window factorial closure in three variables"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every nonzero exponent vector v in N_0^3 and every d>=1, the first
  three factorial moments detect every polynomial
  A+B(X^v)^d+C(X^v)^(2d).  Gauss multiplication identifies the factorial
  weights on this monomial ray, up to a geometric gauge, with a finite
  nonempty product of positive Gamma shapes; current PROVED THM-3107 then
  closes the initial width-three cell.  This is a genuine sparse subclass of
  FC(3), but only a first-window result and not the ambient SFC(3).
audit: >
  An independent audit derived Gauss multiplication by partitioning the
  factors 1..qn into residue classes, checked the geometric gauge and exact
  hypotheses of current proved THM-3107, and tested the transfer and both
  orientation invariants on out-of-box rays with 49 and 55 Gamma layers.
  Normal, optimized, and stored transcripts and declared hashes agree.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-3107-finite-layer-product-gamma-consecutive-width-three-orientation
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
script: 04-computation/factorial_monomial_ray_fc3_thm3125.py
output: 05-knowledge/results/factorial_monomial_ray_fc3_thm3125.out
script_sha256: cbe521dfd4455f5c508c6ce01127f004d084e48206517b54bf5b9d8e0caef4be
output_sha256: 88b4f48e6fbce6659080907e2840d1b82aac0bd4b79e6a10829d6114bd24d56f
hash_basis: LF-normalized bytes
---

# THM-3125 -- monomial-ray first-window factorial closure in three variables

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let the three-variable factorial functional be

```text
L_3 : C[X_1,X_2,X_3] -> C,
L_3(X_1^alpha_1 X_2^alpha_2 X_3^alpha_3)
   =alpha_1! alpha_2! alpha_3!.                              (1)
```

This theorem embeds the all-finite-layer product-Gamma orientation of
THM-3107 into `FC(3)`.  It gives a genuine ambient three-variable subclass,
as opposed to the univariate three-slot computations elsewhere in the
factorial lane.

## 1. Statement

Fix

```text
v=(v_1,v_2,v_3) in N_0^3 \ {(0,0,0)},
s=X^v=X_1^v_1 X_2^v_2 X_3^v_3,
d>=1.                                                        (2)
```

For arbitrary `A,B,C in C`, put

```text
F=A+B s^d+C s^(2d).                                         (3)
```

Then

```text
L_3(F)=L_3(F^2)=L_3(F^3)=0        implies       A=B=C=0.     (4)
```

In particular, if `L_3(F^m)=0` for every `m>=1`, then `F=0`.  Thus every
polynomial in the family `(3)` satisfies the conclusion of the original
ambient three-variable Factorial Conjecture `FC(3)`.

## 2. Exact Gauss-multiplication transfer

Set `t=s^d` and restrict `L_3` to the polynomial algebra `C[t]`.  Its moment
sequence is

```text
W_n=L_3(t^n)=product_(i=1)^3 (d v_i n)!,                    (5)
```

where a zero coordinate contributes `0!=1`.  For every positive integer
`q`, Gauss multiplication gives the exact identity

```text
(qn)!=q^(qn) product_(k=1)^q (k/q)_n,                       (6)
```

with `(theta)_n=theta(theta+1)...(theta+n-1)`.  Apply `(6)` to

```text
q_i=d v_i                    for every i with v_i>0.         (7)
```

Because `v!=0`, the resulting shape list is finite and nonempty, and every
shape `k/q_i` is strictly positive.  Therefore

```text
W_n=G^n V_n,
G=product_(i:v_i>0) (d v_i)^(d v_i)>0,
V_n=product_(i:v_i>0) product_(k=1)^(d v_i) (k/(d v_i))_n.   (8)
```

The sequence `(V_n)` is exactly a finite product of positive Gamma rising
factorials, the universe of THM-3107.

The geometric factor in `(8)` is not discarded informally.  Define

```text
L_V(t^n)=V_n,                       H(t)=A+Bt+Ct^2.           (9)
```

For every polynomial `P`, coefficientwise use of `(8)` gives

```text
L_3(P(s^d))=L_V(P(Gt)).                                     (10)
```

Since substitution `P(t)->P(Gt)` is an algebra homomorphism, for every
`m>=1`,

```text
L_3(F^m)=L_V(H(Gt)^m),
H(Gt)=A+(BG)t+(CG^2)t^2.                                   (11)
```

Current PROVED, VERIFIED-EXACT, independently hostile-audited THM-3107 says
that the initial support `{0,1,2}` is good for every finite nonempty product
of positive Gamma shapes.  Apply it to `(V_n)` and the polynomial `H(Gt)`.
The first three equalities in `(11)` force

```text
A=BG=CG^2=0.                                                (12)
```

Because `G>0`, `(12)` is equivalent to `A=B=C=0`, proving `(4)`.

## 3. Connection contract

The transfer has the following exact type signature.

```text
source:    finite product-Gamma initial support {0,1,2}
target:    three-variable factorial monomial ray {1,s^d,s^(2d)}
map:       t=s^d, followed by the Gauss gauge t -> Gt
preserves: the first-three common-zero predicate
sidecar:   v, d, the positive shape list, and G
loses:     off-ray supports, translated ray supports, and shifted windows.
```                                                                  (13)

No Gaussian functional is used.  The reason the transfer works is the exact
factorization of every factorial `(d v_i n)!` into positive rising-factorial
layers, not a radial heuristic or a moment asymptotic.

## 4. Exact controls

The companion verifies `(6)--(8)` with rational arithmetic on the complete
box

```text
0<=v_i<=3, v!=0;              1<=d<=4;              0<=n<=7,
```

comprising `2016` Gauss cells.  It also rebuilds THM-3107's two normalized
quadratic--cubic orientation invariants directly from the factorial ray
weights for all `63*4=252` pairs `(v,d)` in the same exponent and dilation
box.  Both invariants are strictly negative in every cell.  The transcript
hashes all exact rational values, rejects `v=0`, and contains a hostile in
which omitting one copy of the geometric factor is detected.

These bounded checks verify the specialization and catch indexing or gauge
errors.  The universal quantifiers in `(4)` come from the exact identities
`(6)--(11)` and the all-finite-layer theorem THM-3107, not from extrapolation
of the finite box.

## 5. Boundary and nonclaims

This theorem proves a first-window statement only.  It does not prove:

- every shifted three-moment window for `(3)`, hence not ambient `SFC(3)`;
- `FC(3)` for polynomials outside a single anchored monomial ray;
- a translated ray family `s^a(A+Bs^d+Cs^(2d))`, because the powers carry
  different Gamma tilts, exactly the translation boundary in THM-3107;
- arbitrary three-slot supports not in arithmetic progression on one ray;
- any result for signed or non-Gamma moment sequences after the positive
  layer factorization is lost;
- the false three-dimensional Gaussian moment conjecture, the Jacobian
  conjecture, or any implication between those conjectures and `FC(3)`.

The zero exponent vector is excluded because it collapses all three displayed
monomials to the constant `1`.  Zero coordinates inside a nonzero `v` are
allowed and simply contribute no Gamma layers.

## 6. Reproduction

Run

```text
python3 04-computation/factorial_monomial_ray_fc3_thm3125.py
python3 -O 04-computation/factorial_monomial_ray_fc3_thm3125.py
```

Both executions must byte-match
`05-knowledge/results/factorial_monomial_ray_fc3_thm3125.out` after LF
normalization.

**End of proof.**
