---
id: THM-3859
title: "Marked-root polynomial graphs force a punctured companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every polynomial
  graph component of the
  depressed-cubic branch has an exact marked-root normal form.  When the
  transverse quotient depends only on the graph coordinate, the companion
  is a primitive quadratic: if irreducible it has a finite pole and an
  infinity place on the same normalization; if reducible it contains a
  canonical G_m component.  Thus the graph never occurs without a component
  missing at least two projective points.
source: jc_sparse_direct_search / post-THM-3852 marked-root graph lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit rederived the
  UFD square/cube classification and its divisibility at A=0, checked
  primitivity and the square-discriminant dichotomy over k(A), and followed
  both normalization-place arguments.  In particular, at the negative
  special-fibre root C has a genuine pole, while in the reducible case the
  factor with constant term 18 makes A a unit and has coordinate ring
  k[A,A^(-1)].  Repeated discriminant roots do not repeat a global factor
  because D(0)=81.  The companion verifies the cusp
  square/cube lift, universal factorization, explicit residual quadratic,
  discriminant flanks, special fibre, difference-of-squares factorization,
  the excluded axis boundary, and both nonconstant and THM-3852 collision
  controls.  Normal and optimized runs byte-match the frozen 27-gate
  transcript and both recorded hashes.
  A second independent hostile audit (root, 2026-08-23) rederived the
  classification without using the primary checker, parametrized every
  square-discriminant companion, classified the complete constant boundary
  q=2s (including repetition of the marked graph), checked the excluded axis,
  and separated the fixed-nonzero-product Laurent sheet from this polynomial
  graph grammar.  Its assertion-free 43-gate transcript byte-matches under
  normal, optimized, and frozen replay.
related:
  - THM-3842-nonlinear-cubic-tower-trace-shift-eightfold-base-change
  - THM-3850-nonconstant-cubic-profile-irreducible-branch-puncture-formula
  - THM-3852-affine-two-variable-cubic-profile-line-factor-companion-no-go
  - THM-3851-reciprocal-cubic-discriminant-toric-root-sheet
script: 04-computation/jc2_marked_root_polynomial_graph_companion_thm3859.py
output: 05-knowledge/results/jc2_marked_root_polynomial_graph_companion_thm3859.out
script_sha256: f9ecaa52173f8d1f941332f82132c0fb52ea16a71e173cee18ea516702c5b0a4
output_sha256: c598ffa7afbcde842ee0ec353d4f05ea9c3a711dcca4b8b8c27fbd809694cf33
semantic_sha256: a22900e32db67db55945db2bae117c0fb6b1141ca1eb1f296d89fdcdc1d9e442
hash_basis: raw LF bytes
independent_script: 04-computation/jc2_marked_root_polynomial_graph_companion_independent_audit_thm3859.py
independent_output: 05-knowledge/results/jc2_marked_root_polynomial_graph_companion_independent_audit_thm3859.out
independent_script_sha256: 473cf6ffddcb8a4632c0d1ab398d92b89ad9287782d4da4dc4cc4cde1a0a8d9a
independent_output_sha256: e8145be241a0750d34f68f29a9f3c92769a48ee63e0200e8f7ae124a88bdc6f6
independent_semantic_sha256: 7478af841a342c046c64b903fbfbfc55fd98d7a225721147f5dfa9f29fb46987
---

# THM-3859 -- marked-root graphs force a punctured companion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  For any
`b=b(A,C) in k[A,C]`, put

```text
p=3/2+AC,                         u=1+AC+A^2b,                    (1)

Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                         (2)
```

Then

```text
8p^3-27u^2=A^2 Delta_b.                                         (3)
```

The theorem has a complete classification part and an all-degree companion
obstruction part.

### Marked-root classification

The graph

```text
Gamma_g: C=g(A),                       g in k[A],                 (4)
```

is an irreducible component of `V(Delta_b)` if and only if there are unique
`s in k[A]` and some `Q in k[A,C]` such that

```text
F=C-6s(1+As),
b=2s^2(3+4As)+QF.                                                (5)
```

Here `F=0` is exactly `(4)`.  Thus `(5)` classifies all branch components
that are polynomial graphs over the `A`-axis, not merely a selected ansatz.

### One-variable transverse obstruction

Assume in addition that the transverse quotient in `(5)` lies in `k[A]`:

```text
Q=q(A).                                                          (6)
```

Then

```text
Delta_b=F H_(s,q),                                               (7)
```

where `H_(s,q)` is primitive and has no vertical component.  At least one
irreducible component of `V(H_(s,q))` has affine normalization obtained from
its smooth projective normalization by deleting at least two distinct
points.  More sharply:

- if `H_(s,q)` is irreducible, one deleted point lies over `A=0` and at
  least one lies over `A=infinity`;
- if `H_(s,q)` is reducible, its reduced packet contains a canonical
  component isomorphic to `G_m`, hence exactly `P1 minus {0,infinity}`.

Repeated roots of the discriminant may make fibres singular or make the two
components of the **quadratic companion `H_(s,q)`** meet nontransversely, but
they do not create a repeated irreducible factor inside `H_(s,q)`: its
generic discriminant is nonzero.  The two companion cases above are
exhaustive.  The total branch `Delta_b=F H_(s,q)` can still repeat the marked
graph when `F` also divides `H_(s,q)`; the zero control `(25)` does exactly
this and is not being excluded.

Consequently a polynomial graph branch in the one-variable transverse
grammar `(5)--(6)` never occurs without a nonpolynomial companion.  The
theorem makes no claim for `C`-dependent `Q`, nor for the vertical axis
`A=0`; those are the precise surviving escape directions.

## 1. The cusp UFD forces the marked-root law

Restrict `(3)` to the graph `(4)`, whose coordinate ring is the UFD `k[A]`.
Set

```text
P=1+(2/3)Ag,              U=1+Ag+A^2b(A,g).                     (8)
```

Then `P^3=U^2`.  The polynomial `P` is nonzero: `P=0` would require
`Ag=-3/2`, impossible in `k[A]`.  Unique factorization, with scalar roots
absorbed because `k` is algebraically closed, gives a unique `z in k[A]`
such that

```text
P=z^2,                         U=z^3.                            (9)
```

Equations `(8)--(9)` imply

```text
Ag=(3/2)(z-1)(z+1),
A^2b(A,g)=(1/2)(z-1)^2(2z+1).                                  (10)
```

The second identity says its right side is divisible by `A^2`.  If
`z(0)=-1`, that right side has nonzero value `-2`, a contradiction.  Hence
`z(0)=1`, so there is a unique `s in k[A]` with

```text
z=1+2As.                                                        (11)
```

Substitution in `(10)` and cancellation of the nonzero polynomial `A`
give

```text
g=6s(1+As),
b(A,g)=2s^2(3+4As).                                             (12)
```

Since the monic polynomial `F=C-g(A)` generates the graph ideal, `(12)` is
equivalent to the division statement `(5)`.  Conversely `(11)--(12)` give
`P=z^2,U=z^3`, so `(3)` proves `F | Delta_b`.  This establishes the iff and
the uniqueness assertion.

The same calculation works on the normalization of any non-axis polynomial
branch, but `s` need not descend to `k[A]` unless the branch is a graph.
This is exactly why `(4)` is retained in the theorem statement.

## 2. Exact companion and its two discriminant flanks

Under `(5)--(6)`, direct division of `(2)` by `F` gives

```text
H=8AC^2+BC+E,                                                    (13)

B=-27A^2q^2+48A^2s^2-54Aq+48As+9,                              (14)

E=162A^3q^2s^2-432A^3qs^3+288A^3s^4
 +162A^2q^2s-324A^2qs^2+144A^2s^3
 +18As^2-54q+54s.                                               (15)
```

Its discriminant as a quadratic in `C` factors as

```text
D=B^2-32AE
 =27(Aq+4As+3)(3Aq-4As+1)^3.                                   (16)
```

In particular,

```text
H(0,C)=9C-54q(0)+54s(0),
B(0)=9,                         D(0)=81.                         (17)
```

The coefficients of `H` are primitive because its leading coefficient is
`8A` while `B(0)=9`.  Also `A` does not divide `H`, and no `A-a` with
`a!=0` can divide it because the `C^2` coefficient is `8a`.  Thus `H` has
no vertical component.  Since `D(0)!=0`, its generic quadratic has distinct
roots; a repeated global factor is impossible.

## 3. Irreducible companion: the two places are forced

Suppose `D` is not a square in `k(A)`.  By primitivity, `H` is irreducible.
Its function field is the quadratic extension described by

```text
w=16AC+B,                       w^2=D,
C=(w-B)/(16A).                                                   (18)
```

Because `D(0)=81`, the smooth projective normalization has two distinct
unramified points over `A=0`, with `w=9` and `w=-9`.  At the second point,

```text
w-B=-18+O(A),                                                     (19)
```

so `C` has a pole.  This point is absent from the affine normalization.
The other point is the unique finite point visible in `(17)`; its `C` value
is `6(q(0)-s(0))`.

The finite morphism from the smooth projective normalization to `P1_A` has
at least one point over `A=infinity`.  The coordinate `A` has a pole at every
such point, so all of them are absent from the affine normalization.  A
point over zero and a point over infinity are distinct.  Since the
normalization is connected in the irreducible case, this proves the required
two-puncture lower bound on one component.

## 4. Reducible companion: the bad component is exactly `G_m`

Suppose `D` is a square in `k(A)`.  Since `k[A]` is a UFD and `D` is a
polynomial, after changing sign there is `W in k[A]` with

```text
W^2=D,                          W(0)=9.                          (20)
```

Write

```text
B=9+AB1,                       W=9+AV.                           (21)
```

The difference-of-squares identity

```text
32AH=(16AC+B-W)(16AC+B+W)                                      (22)
```

has first factor divisible by `A`.  Cancelling that factor gives the exact
polynomial factorization

```text
32H=(16C+B1-V) J,

J=18+A(16C+B1+V).                                               (23)
```

Both factors are primitive and linear in `C`, hence irreducible.  They are
distinct because `D` is not the zero polynomial.  The first is a polynomial
graph.  On the second, `A` is a unit and

```text
C=-(18/A+B1+V)/16,
```

so

```text
k[A,C]/(J) ~= k[A,A^(-1)].                                      (24)
```

Thus `J=0` is canonically `G_m`; its smooth projective normalization is
`P1_A`, and its affine normalization deletes exactly the distinct points
`A=0` and `A=infinity`.  If the sign in `(20)` is reversed, the two factors
in `(22)` swap roles and the same component is recovered.  Zeros of `W` may
make the components meet, but do not change `(24)`.

## 5. Boundaries and controls

The zero choice `s=q=0` gives

```text
F=C,                b=0,
Delta_0=C^2(8AC+9).                                              (25)
```

Here the graph itself repeats, while its companion contains the predicted
`G_m` hyperbola.  The genuinely nonconstant choice `s=A,q=0` gives

```text
F=C-6A(1+A^2),
D=27(4A^2+3)(1-4A^2)^3,                                        (26)
```

an irreducible-companion control with the two forced places of Section 3.

The unique doubled-line boundary of THM-3852 is recovered exactly by

```text
s=-1/6,                         q=-1/3.                          (27)
```

Indeed, with `L=A-6C-6` and
`K=A^2-24AC-6A-27`,

```text
F=-L/6,
b=(A-18C-9)/54,
H=LK/18,                        D=(A-3)^4.                       (28)
```

Thus the smooth two-place conic from THM-3852 is precisely the canonical
`G_m` factor `(24)` in the square-discriminant marked-root chart.

The entire constant-profile boundary is just as rigid.  Put `s=a` and
`q=d`, with `a,d in k`.  Then

```text
D=27(3+A(d+4a))(1+A(3d-4a))^3.                                (29)
```

The two odd-multiplicity linear flanks have the same root exactly when
`d=2a`; the case where both are constant is `a=d=0` and is already included.
Thus the constant companion is reducible exactly when `q=2s`.  On this
complete square boundary,

```text
H=(6Aa^2-C+6a)(12A^2a^2-8AC+12Aa-9).                         (30)
```

The first factor is `-F`, so the total branch repeats the marked graph, and
the second is the canonical `G_m` component.  Constants with `d!=2a` lie in
the irreducible case.

The excluded axis boundary is also explicit: `A=0` is a branch component
exactly when

```text
b(0,C)=C^2/6.                                                    (31)
```

It is not a graph over `A` and is not covered by this theorem.  Likewise a
general graph profile has `Q=Q(A,C)` in `(5)`; allowing genuine `C`
dependence is the next open transverse operation.

This boundary is disjoint from the fixed-nonzero-product toric root sheet of
THM-3851.  A root of `T^3-UT^2+VT-c` with `c!=0` satisfies

```text
V=cT^(-1)-T^2+UT,                                               (32)
```

so its sheet is a Laurent graph with source `A1 x G_m`, not a polynomial
graph with source `A1`.  At `c=0` the pole disappears only because the cubic
becomes reducible.  Hence the toric sheet neither evades this theorem nor
fills its genuinely `C`-dependent transverse direction.

## 6. Exact replay

Run

```bash
python3 04-computation/jc2_marked_root_polynomial_graph_companion_thm3859.py
python3 -O 04-computation/jc2_marked_root_polynomial_graph_companion_thm3859.py
python3 04-computation/jc2_marked_root_polynomial_graph_companion_independent_audit_thm3859.py
python3 -O 04-computation/jc2_marked_root_polynomial_graph_companion_independent_audit_thm3859.py
```

The first two commands must byte-match
`05-knowledge/results/jc2_marked_root_polynomial_graph_companion_thm3859.out`.
The latter two must byte-match
`05-knowledge/results/jc2_marked_root_polynomial_graph_companion_independent_audit_thm3859.out`.
The assertion-free companions perform respectively 27 and 43 exact
structural and hostile gates.  Their raw-LF SHA-256 hashes are recorded in
the metadata.
