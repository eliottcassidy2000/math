---
id: THM-3570
title: "Universal Pell-conic target-graph factor compiler"
status: >
  PROVED + VERIFIED-EXACT.  For every nonzero phi in C[a,b], the fixed
  THM-1300 target-graph core cubic is reducible over C(a,b) if and only if
  there is q in C(a,b)^* with phi=4(q^2+bq+4a)/q^3.  The associated rational
  root is x=-2q/(q^2+2bq+12a).  The converse follows from the exact
  quadratic-in-phi discriminant, whose nonsquare part is the Pell conic
  W^2=(1+bx)^2-12ax^2; its two rational branches are identified by
  q -> 12a/q.  The equivalent search equation is
  phi*q^3-4q^2-4bq-16a=0.
source: kps-s188
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - THM-3565-resonant-linear-a-target-graph-factor-classification
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
companion: 04-computation/jacobian_universal_pell_conic_factor_compiler_kps_s188.py
output: 05-knowledge/results/jacobian_universal_pell_conic_factor_compiler_kps_s188.out
script_sha256: d633a3f3f6191aa48524d4519ecc61374a411980695ce2b28b0f97f0310ec8f1
output_sha256: fa142160b424be9d6432debef441a1d4d40694ea8d96c9b4178033d0ebb66acc
hash_basis: LF-normalized bytes
---

# THM-3570 -- universal Pell-conic target-graph factor compiler

**PROVED + VERIFIED-EXACT.**  Reducibility of every polynomial target graph
for the fixed THM-1300 map is governed by one rational Pell conic.  This
turns the sparse degree-resonance condition of THM-3564 into an exact
all-degree factor compiler.

Put `K=C(a,b)`.  The theorem concerns `0!=phi in C[a,b]`; the excluded row
`phi=0` has the separate visible factor `F3=x(2-3xy-x^2z)`.

## 1. Statement

THM-2473's target-graph core cubic is

```text
E_phi(x)=L_phi x^3+(4+3bphi)x+2phi,

L_phi=27a^2phi^2+18abphi+16a-b^3phi-b^2.                 (1)
```

For `phi!=0`, the following are equivalent:

1. `E_phi` is reducible over `K`;
2. there is `q in K*` satisfying

```text
phi=4(q^2+bq+4a)/q^3;                                    (2)
```

3. the cubic equation

```text
phi q^3-4q^2-4bq-16a=0                                   (3)
```

has a nonzero root in `K`.

For a parameter in `(2)`, one rational root of `(1)` is

```text
x=-2q/(q^2+2bq+12a).                                     (4)
```

The denominator in `(4)` cannot vanish for `q in K`: as a quadratic in
`q`, its discriminant is `4(b^2-12a)`, which has odd valuation at the
irreducible divisor `b^2-12a` and is nonsquare in `K`.

Thus polynomial target-graph factor searches reduce exactly to rational
solutions of `(3)`.  No degree bound or genericity hypothesis is present.

## 2. Solve the core equation for the target graph

Fix a putative root `0!=x in K` and regard `E_phi(x)=0` as a quadratic in
`phi`:

```text
27a^2x^3 phi^2
 +[(18ab-b^3)x^3+3bx+2]phi
 +(16a-b^2)x^3+4x=0.                                    (5)
```

Its discriminant factors exactly as

```text
Delta_phi
=-[12ax^2-b^2x^2-2bx-1]
  [12ax^2-b^2x^2+bx+2]^2.                               (6)
```

Set

```text
S=12ax^2-b^2x^2+bx+2.                                   (7)
```

The divisor `S=0` contains no `K`-rational `x`.  As a quadratic in `x`, its
discriminant is

```text
3(3b^2-32a),                                             (8)
```

which has odd valuation one at the irreducible divisor `3b^2-32a` and is
therefore not a square in `K`.  Hence `S!=0` at every rational root under
consideration.

Because `(5)` has the rational solution `phi`, its discriminant is a square.
After dividing out the nonzero square `S^2`, there is `W in K` with

```text
W^2=(1+bx)^2-12ax^2.                                     (9)
```

This is the promised Pell conic: it is the norm equation for the two
linear factors of `(1+bx)^2-12ax^2`, and it has the rational base point
`(x,W)=(0,1)`.

## 3. Rational parametrization and branch involution

Since the actual core root has `x!=0`, draw the line through the base point

```text
W=1+tx,                         t in K.                  (10)
```

Substitution in `(9)` gives

```text
x=2(b-t)/(t^2-b^2+12a).                                  (11)
```

The case `t=b` would give `x=0`, so put

```text
q=t-b in K*.                                              (12)
```

Then `(11)` becomes `(4)`.  Substitution of `(10)--(12)` in the two
quadratic-formula branches of `(5)` gives

```text
phi_1=q(q^2+3bq+36a)/(108a^2),

phi_2=4(q^2+bq+4a)/q^3.                                  (13)
```

These are not distinct families.  Direct simplification gives

```text
phi_1(q)=phi_2(12a/q).                                    (14)
```

Since `12a/q` is again a nonzero element of `K`, every rational root of
`(1)` yields a parameter satisfying `(2)`.  This proves `1 => 2`; equations
`(2)` and `(3)` are visibly equivalent.

## 4. Direct converse

Conversely, take any `q in K*`, define `phi` by `(2)`, and define `x` by
`(4)`.  Clearing denominators and expanding gives

```text
L_phi x^3+(4+3bphi)x+2phi=0.                             (15)
```

Thus the cubic has a rational root and is reducible.  This proves `2 => 1`
and completes the equivalence.

The argument also handles points where one particular conic chart looks
singular: the involution `(14)` exchanges the two charts.  The only omitted
core-root case is `x=0`, which by `(1)` is exactly the explicitly excluded
row `phi=0`.

## 5. Relation to the first resonance and the Pell tower

THM-3565's complete `deg_a(phi)=1` family is the specialization

```text
q=-2/h(b).                                                (16)
```

Indeed `(2)` and `(4)` become

```text
phi=-2h^3a+b h^2-2h,

x=h/(3ah^2-bh+1).                                        (17)
```

Thus THM-3565 is the polynomial-descent locus of the universal conic, and
THM-3568 closes every source-coordinate factor on that locus.

Equation `(9)` is also the exact shared mechanism with the Chebyshev--Pell
program reserved in THM-3566: iterating rational norm parameters can create
odd collision towers, while polynomial target graphs impose the additional
divisibility condition `(3)`.  This is a precise map of objects and lost
information, not an inference that the reserved tower is proved.

For the next resonant degrees `4,7,...`, the sharp search is no longer an
unstructured factorization of a large pullback.  It is:

```text
find q in C(a,b)* such that 4(q^2+bq+4a)/q^3 lies in C[a,b]
with the prescribed degree and collision value.          (18)
```

The theorem does not classify those polynomial-descent loci, determine
whether either source component is `A2`, or prove a planar Jacobian
counterexample.

## 6. Exact verification

Run

```bash
python3 04-computation/jacobian_universal_pell_conic_factor_compiler_kps_s188.py
python3 -O 04-computation/jacobian_universal_pell_conic_factor_compiler_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion verifies the
discriminant factorization `(6)`, exceptional discriminant `(8)`, conic
parametrization, both branches `(13)`, involution `(14)`, direct converse
`(15)`, cubic equation `(3)`, and the THM-3565 regression `(16)--(17)`.  Five
additional rational controls include parameters depending on both `a` and
`b`.

**QED.**
