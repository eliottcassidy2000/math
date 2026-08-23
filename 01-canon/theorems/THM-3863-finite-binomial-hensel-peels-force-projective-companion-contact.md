---
id: THM-3863
title: "Finite binomial Hensel peels force projective companion contact"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  The unique formal binomial lift of the cusp
  admits a polynomial marked-graph quotient coefficient by coefficient, but
  every finite truncation leaves an exact residual C-degree drop N over
  A=0.  A reduced component through that projective contact dominates the
  A-line and therefore loses both its point over finite A=0,C=infinity and a
  distinct point over A=infinity.
source: jc_sparse_direct_search / THM-3859--THM-3860 Newton-edge synthesis, 2026-08-23
audit: >
  SELF-AUDITED exact proof candidate.  The companion constructs the quotient
  recurrence and replays N=1 through N=6, checking every Euclidean remainder,
  coefficient degree and leading term, A-adic factorization, residual degree
  drop, exact projective Z^N contact, and the normalized N=2 control.  Normal
  and optimized runs byte-match the frozen 96-gate transcript.  Independent
  hostile audit is still required before promotion.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
script: 04-computation/jc2_finite_binomial_hensel_peel_thm3863.py
output: 05-knowledge/results/jc2_finite_binomial_hensel_peel_thm3863.out
script_sha256: 00b487568b34de58d606662816829e17b0d87088d7a3797e67e2d089652ee902
output_sha256: 4aba045eafd1c8b37890bbb427a8adcdaca4a5ebcea384b5abb75c99492875a8
semantic_sha256: 3128c19fcb8a5843975d52dc7202abc68a247d38af2d8ad1f143589371e2437f
hash_basis: raw LF bytes
---

# THM-3863 -- every finite binomial peel leaves projective contact

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
P=1+(2/3)AC,

b_*={P^(3/2)-1-AC}/A^2
   =sum_(j>=0) beta_j A^j C^(j+2),                              (1)

beta_j=binom(3/2,j+2)(2/3)^(j+2).                              (2)
```

Here `P^(3/2)` is the unique `A`-adic binomial series with constant term
one.  Every `beta_j` is nonzero.

Fix any `sigma in k` and define the polynomial graph and its marked-root
profile

```text
F=C-6sigma(1+A sigma),
b0=2sigma^2(3+4A sigma).                                       (3)
```

There is a unique formal quotient

```text
Q_*=(b_*-b0)/F=sum_(j>=0) q_j(C)A^j       in k[C][[A]],         (4)
```

and its coefficients satisfy

```text
q_j in k[C],                   deg q_j=j+1,
lc_C(q_j)=beta_j.                                                (5)
```

For every integer `N>=1`, truncate one step before `q_N`:

```text
Q_N=sum_(j=0)^(N-1) q_j(C)A^j,
b_N=b0+FQ_N.                                                     (6)
```

Then the depressed-cubic branch polynomial

```text
Delta_N=-27A^2b_N^2+8AC^3-54ACb_N+9C^2-54b_N                 (7)
```

has the exact factorization

```text
Delta_N=A^N F R_N,                 R_N in k[A,C],               (8)
```

with sharp degrees

```text
deg_C R_N=2N+1,
R_N(0,C)=54q_N(C),             deg_C R_N(0,C)=N+1.              (9)
```

Let `X_N=V(R_N)`.  At least one irreducible component of
`(X_N)_red` has affine normalization obtained from its smooth projective
normalization by deleting at least two distinct points.  One lies over the
finite base value `A=0` and `C=infinity`; the other lies over
`A=infinity`.  This conclusion is valid even if `R_N` is reducible, has
repeated factors, or several branches share the finite-base projective
contact.

Equivalently, every finite polynomial Hensel peel transfers rather than
removes the companion debt.  Only the infinite, nonpolynomial binomial
series `(1)` makes the branch identity vanish to every `A`-adic order.

This theorem concerns the canonical finite truncations `(6)`.  It does not
classify arbitrary `C`-dependent transverse quotients in THM-3859.

## 1. The binomial series is the unique formal zero-branch lift

For any profile `b`, set

```text
u=1+AC+A^2b.                                                     (10)
```

The universal identity is

```text
A^2 Delta_b=27(P^3-u^2).                                        (11)
```

In the complete ring `k[C][[A]]`, one has `u=1 mod A`.  Therefore
`Delta_b=0` formally if and only if

```text
u=P^(3/2),                                                       (12)
```

the unique square root of `P^3` congruent to one.  Solving `(10)` gives
exactly `b=b_*` in `(1)`.  Formula `(2)` is the binomial coefficient
expansion.  None of its coefficients vanish in characteristic zero, so the
unique formal solution is not a polynomial.

This is a uniqueness statement in the completed ring, not a claim that a
formal series has algebraized.

## 2. The marked graph divides the formal lift coefficientwise

On `F=0`,

```text
P=1+4A sigma(1+A sigma)=(1+2A sigma)^2.                          (13)
```

The chosen binomial branch therefore gives

```text
P^(3/2)=(1+2A sigma)^3,
b_*|_(F=0)=6sigma^2+8A sigma^3=b0.                              (14)
```

Because `F` is monic in `C`, `(14)` gives the unique formal quotient `(4)`
with polynomial coefficients.

There is also a completely elementary recursion.  Put `q_(-1)=0` and

```text
d_j=beta_j C^(j+2)
    -[j=0]6sigma^2-[j=1]8sigma^3.                               (15)
```

Equating the coefficient of `A^j` in `FQ_*=b_*-b0` yields

```text
(C-6sigma)q_j-6sigma^2q_(j-1)=d_j.                              (16)
```

The numerator obtained by solving `(16)` vanishes at `C=6sigma`; this is
also immediate from `(14)`.  Hence `q_j` is polynomial.  Its top term comes
only from `beta_j C^(j+2)`, proving `(5)` inductively.

## 3. Exact A-adic order of every truncation

Equations `(4)` and `(6)` give

```text
b_*-b_N=A^N F(q_N+Aq_(N+1)+...).                               (17)
```

Let `u_*=P^(3/2)` and let `u_N` be `(10)` for `b_N`.  From `(17)`,

```text
u_*-u_N=A^(N+2)F(q_N+O(A)),
u_*+u_N=2+O(A).                                                  (18)
```

The marked-root identity `(13)--(14)` gives `F | Delta_N`, while `(18)`
gives `A^N | Delta_N`.  Since `A` and the monic graph polynomial `F` are
coprime in `k[A,C]`, their product divides.  Substitution in `(11)` then
also gives

```text
R_N(0,C)=54q_N(C).                                               (19)
```

The factor has exact `A`-adic order `N`, because `q_N` is nonzero.

For the generic `C`-degree, `(5)--(6)` give

```text
deg_C b_N=N+1,
lc_C(b_N)=beta_(N-1)A^(N-1).                                   (20)
```

The unique highest `C`-degree term of `(7)` comes from
`-27A^2b_N^2`.  After division by `A^N F`,

```text
deg_C R_N=2N+1,
lc_C(R_N)=-27 beta_(N-1)^2 A^N.                                (21)
```

Together `(5)`, `(19)`, and `(21)` prove all of `(9)`.

## 4. The finite-base projective contact forces two missing points

Homogenize `R_N` only in the `C` coordinate to degree `2N+1`, using
`[C:Z]` on `P1_C`; call the result `mathcal R_N(A;C,Z)`.  Let
`q_N^h(C,Z)` be the degree-`N+1` homogenization of `q_N`.  Equation `(9)`
gives the exact special fibre

```text
mathcal R_N(0;C,Z)=54 q_N^h(C,Z) Z^N.                           (22)
```

At

```text
P_0=(A=0,[C:Z]=[1:0]),                                         (23)
```

one has `q_N^h(1,0)=beta_N!=0`.  Thus the total projective branch has
contact of exactly `N` with `C=infinity` over `A=0`.  This is a contact
multiplicity; it is not asserted to consist of `N` distinct normalization
points.

Factor in the affine UFD

```text
R_N=c product_i G_i^(e_i),                                     (24)
```

with distinct irreducible `G_i`.  The `C`-homogenization of `(24)` is the
product of the individual homogenizations at their actual `C`-degrees, so
`(22)--(23)` put `P_0` on the projective closure of at least one reduced
component `Gamma=V(G_i)`.

That component is not vertical.  Indeed `(19)` is a nonzero polynomial, so
`A` does not divide `R_N`; an irreducible component supported over `A=0`
would have defining factor `A`.  Nor is its projective closure the boundary
`Z=0`, because `G_i` is an affine factor at its actual `C`-degree.  Hence
the projection of `Gamma` to the `A`-line is nonconstant.

Take the closure of `Gamma` in `P1_A x P1_C` and its smooth projective
normalization

```text
nu: Gamma_tilde -> Gamma_bar.                                  (25)
```

Surjectivity of normalization supplies at least one point `p_0` over
`P_0`.  The pullback of `C/Z` has a pole at `p_0`, since `C!=0` and `Z=0`
there.  Therefore `p_0` is absent from the affine normalization of
`Gamma`.

Because `Gamma` dominates the `A`-line, its normalized projective map to
`P1_A` is nonconstant and hence surjective.  It has at least one point
`p_infinity` over `A=infinity`.  The function `A` has a pole there, so this
point is also absent from the affine normalization.  Finally,

```text
A(p_0)=0,                         A(p_infinity)=infinity,        (26)
```

so the two points are distinct.  This proves the componentwise conclusion
without assuming irreducibility or squarefreeness of the total residual.
Multiplicities `e_i` and collisions among branches affect the contact ledger
but not the reduced component selected above.

## 5. Newton-edge synthesis and the first concrete peel

THM-3859 begins with the quadratic Newton edge

```text
8A C^2+B C+E,                    B(0) in k*,                     (27)
```

which sends one sheet to `C=infinity` over `A=0`.  Taking `q_0` in `(6)`
peels that sheet into an `A=0` factor, but `(22)` shows that the next
residual still has projective contact.  Repeating finitely never ends the
debt.

The same valuation mechanism appears in THM-3860's vertical rational lift:
its equation `M(S-f)=s` has positive divisor order on `M` and unit generic
value on `s`, forcing `S` to have negative order.  Here the successive
Newton edges are the polynomial shadows of the unique square root `(1)`.
This is a structural comparison, not a dependency on THM-3860's provisional
global theorem.

For an exact low-degree control, take `N=2` and `sigma=-1/6`.  Then

```text
b_2=(2A^2C^2-2A^2C-A^2-12AC^3+108C^2)/648,                    (28)
```

and

```text
Delta_(b_2)=A^2 F R_2,
deg_C R_2=5,
R_2(0,C)=(3C^3-3C^2+C+1)/12.                                  (29)
```

Thus the first two vertical layers are genuinely paid by a polynomial
profile, yet the residual degree drops from five to three and retains exact
contact two at finite-base infinity.  This is a positive polynomial
near-miss and a hostile control for any claim that a finite jet already
algebraizes the binomial lift.

## 6. Exact replay

Run

```bash
python3 04-computation/jc2_finite_binomial_hensel_peel_thm3863.py
python3 -O 04-computation/jc2_finite_binomial_hensel_peel_thm3863.py
```

Both commands must byte-match
`05-knowledge/results/jc2_finite_binomial_hensel_peel_thm3863.out`.  The
assertion-free companion verifies `N=1,...,6` through 96 exact gates.  The
all-`N` proof is the coefficient recursion and degree argument above; the
finite replay is a hostile control, not its logical substitute.  Raw-LF
SHA-256 hashes are recorded in the metadata.
