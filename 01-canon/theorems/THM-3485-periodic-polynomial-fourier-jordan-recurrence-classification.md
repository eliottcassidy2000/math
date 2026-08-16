---
id: THM-3485
title: "Periodic-polynomial Fourier/Jordan recurrence classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Over any
  characteristic-zero field, the minimal shift
  polynomial of a period-p polynomial word is the product over Fourier
  colours (x-zeta^j)^(deg Q_j+1), omitting zero colours.  Its exact recurrence
  defect from the naive p(d+1) bound is the sum of the missing top Jordan
  rungs.  Over Q the product descends by cyclotomic order.  A shared leading
  lane coefficient removes at least one top rung from every nontrivial colour.
source: codex-2026-08-15-periodic-polynomial-recurrence
audit: >
  independent proof and replay audit verified Fourier/Jordan minimality,
  Galois descent, reduced generating denominator, defect and shared-layer
  formulas, prime/composite-period boundaries, rational cyclotomic grouping,
  the THM-3484 specialization, zero-word boundary, modular Hankel-rank
  minimality, security, and byte-identical normal/optimized/stored replay
depends_on: []
related:
  - THM-3484-ternary-weighted-determinant-minimal-recurrence-and-cubic-fourier-degree-drop
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
script: 04-computation/periodic_polynomial_fourier_jordan_recurrence_thm3485.py
output: 05-knowledge/results/periodic_polynomial_fourier_jordan_recurrence_thm3485.out
script_sha256: 2d6441a7c875ded74b7ff6b83509108da450f6bb18553d37652269258936a6a0
output_sha256: 85859d4b1585c6b72f5224f5022b258230fa4d630317fa078ac943e1f4f43789
semantic_sha256: 834dfd2be5aed2f4d92b6f7fe742bc0df3a3008e61a0e1bc3fd39eef2dac54c0
hash_basis: LF-normalized bytes
---

# THM-3485 -- Fourier colours classify periodic-polynomial recurrences

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem isolates the general mechanism behind THM-3484's order drop.
It is elementary, but keeping the exact zero-colour and descent boundaries
explicit prevents a finite-fit recurrence from being mistaken for a minimal
one.

## 1. Periodic polynomial words and their colours

Let `K` be a characteristic-zero field, let `p>=1`, and choose polynomials

```text
P_0(t),...,P_(p-1)(t) in K[t].                       (1)
```

Define the sequence

```text
a_n=P_r(n) when n==r (mod p),       n>=0.             (2)
```

In a splitting field `L` of `x^p-1` over `K`, fix a primitive `p`th root
`zeta` and put

```text
Q_j(t)=1/p sum_(r=0)^(p-1) zeta^(-jr)P_r(t),
                                                     0<=j<p.       (3)
```

For `Q_j!=0`, write `d_j=deg Q_j`; zero colours have no degree.  The finite
Fourier indicator identity gives, for every `n>=0`,

```text
a_n=sum_(j=0)^(p-1) zeta^(jn) Q_j(n).                (4)
```

Fourier inversion also shows that all `Q_j` vanish exactly when all `P_r`
vanish.

## 2. Exact minimal shift polynomial

Let `E` denote forward shift, `(Ea)_n=a_(n+1)`.  The theorem is

```text
boxed:
chi_a(x)=product_(j:Q_j!=0) (x-zeta^j)^(d_j+1).      (5)
```

This is the minimal annihilator of `(2)` under shifts—equivalently, the
minimal polynomial on its cyclic shift module—over `L` and over `K`.
Equivalently, the reduced ordinary generating denominator is

```text
D_a(z)=product_(j:Q_j!=0) (1-zeta^j z)^(d_j+1).      (6)
```

### Proof

For a nonzero polynomial `R` of degree `e`, direct subtraction gives

```text
(E-lambda)(lambda^n R(n))
 =lambda^(n+1)(R(n+1)-R(n)).                          (7)
```

In characteristic zero the finite difference on the right has degree exactly
`e-1`.  Hence `lambda^nR(n)` is killed by `(E-lambda)^(e+1)` and by no smaller
power.

For distinct nonzero `lambda`, the spaces of polynomial-exponential
sequences `lambda^nR(n)` are generalized eigenspaces of `E` with pairwise
coprime minimal polynomials.  Their sum is direct: applying all other primary
factors isolates any chosen eigenspace, where those factors are invertible.
Therefore the minimal polynomial of the sum `(4)` is the least common
multiple of the primary factors, which is their product `(5)`.

Every automorphism of `L/K` permutes the roots `zeta^j` and sends the
corresponding `Q_j` to the permuted Fourier colour.  It preserves each degree
and hence fixes `(5)`.  Thus `chi_a` lies in `K[x]`.  A smaller polynomial over
`K` would remain a smaller annihilator after extension to `L`, contradicting
the already proved minimality there.  This also proves `(6)`, since distinct
colours give distinct poles and `(7)` gives their exact orders.  QED.

The proof applies unchanged to eventual recurrences: a nonzero polynomial in
`n` cannot vanish on an infinite tail, so deleting a finite prefix removes no
primary factor.

## 3. Exact defect from the naive bound

Assume in this section that the word is nonzero, and let

```text
d=max_(r:P_r!=0) deg P_r.                            (8)
```

The usual quasi-polynomial bound is

```text
chi_a divides (x^p-1)^(d+1),
order(chi_a)<=p(d+1).                                (9)
```

For each colour define its missing-rung count

```text
delta_j = d-deg Q_j       if Q_j!=0,
delta_j = d+1             if Q_j=0.                  (10)
```

Then `(5)` gives the exact identity

```text
boxed:
p(d+1)-deg chi_a=sum_(j=0)^(p-1) delta_j.             (11)
```

Thus every recurrence compression of a periodic polynomial word is exactly a
Fourier-colour degree loss.  A missing colour deletes its entire Jordan
block; cancellation of one leading layer removes one top rung.

More explicitly, let

```text
c_m=( [t^m]P_0,...,[t^m]P_(p-1) ).                   (12)
```

Then `d_j` is the largest `m` for which the `j`th finite Fourier coefficient
of `c_m` is nonzero.  This is a complete layer-by-layer classifier, not only a
leading-term test.

## 4. Shared leading layers

Suppose the degree-`d` lane coefficients are one common nonzero scalar:

```text
c_d=(c,c,...,c),        c!=0.                         (13)
```

The trivial colour retains degree `d`, while every nontrivial colour loses
that layer.  Hence

```text
deg chi_a <= p(d+1)-(p-1).                           (14)
```

Equality holds exactly when every nontrivial colour is nonzero in degree
`d-1`.

More generally, for `1<=h<=d`, if the top `h` vectors

```text
c_d,c_(d-1),...,c_(d-h+1)                            (15)
```

are constant across lanes, then every nontrivial colour has degree at most
`d-h`, and at least `h(p-1)` orders disappear.  For rational lanes and prime
`p`, if `c_(d-h)` is nonconstant, Galois transitivity makes all nontrivial
Fourier coefficients at that layer nonzero.  The drop is then exactly
`h(p-1)`.

The sharp opposite boundary is also useful: if every coefficient vector
`c_m` is constant, all nontrivial colours vanish and `(2)` is just one
ordinary polynomial sequence, with minimal polynomial `(x-1)^(d+1)`.

## 5. Rational cyclotomic form

For `K=Q`, colours whose roots have the same multiplicative order form one
Galois orbit and have the same exponent.  For each `m|p`, let `e_m` be
`d_j+1` for any `j` such that `zeta^j` has order `m`, or zero if that orbit is
absent.  Then `(5)` descends explicitly as

```text
boxed:
chi_a(x)=product_(m|p) Phi_m(x)^e_m.                 (16)
```

Composite periods can therefore lose selected cyclotomic orders entirely.
The `p=4` parity-only hostile in the companion has colours of orders one and
two but no order-four colour; its recurrence order is six rather than the
naive sixteen.  This is why “one common leading coefficient” is not the whole
classification.

## 6. THM-3484 is the first nontrivial cubic instance

For the three determinant lanes of THM-3484,

```text
(d_0,d_1,d_2)=(7,6,6).                               (17)
```

Equations `(5)` and `(16)` give

```text
chi(x)=(x-1)^8 Phi_3(x)^7
      =(x-1)^8(x^2+x+1)^7,                           (18)
```

of order `22=24-2`.  The two missing orders are exactly one top Jordan rung
in each nontrivial cubic colour.  THM-3484 additionally supplies the explicit
determinant word, recurrence numerator, shortened-factor hostiles, and
nonzero Hankel determinant; none of those instance-specific calculations is
needed for the general proof above.

## 7. Relation to jets, tournaments, and harmonic subsets

Equation `(7)` identifies polynomial degree with generalized-eigenvector or
Jordan depth.  Fourier cancellation removes top rungs before the recurrence
is interlaced.  THM-3476 has a parallel grammar: powers of its evaluation
polynomial `P` are exactly transverse Hasse-jet order.  The two theorems are
not instances of one another.  THM-3476 works in characteristic two with a
formal evaluation kernel and proves that physical first-live jet depth is
unbounded; the present theorem works in characteristic zero with a fixed
finite period and classifies a finite shift module.

The tournament-sized-four input in THM-3484 enters only through the two `K4`
tree polynomials that create its degree-seven lanes.  Formula `(5)` sees the
resulting periodic word, not the graph or a tournament orientation.

Likewise each residue class modulo `p` has natural and harmonic coefficient
`1/p`, but an arbitrary subset of the natural numbers need not have either.
Periodicity is the sidecar that makes both the recurrence and harmonic
coefficient finite-state; `(5)` does not classify nonperiodic harmonic
subseries.

## 8. Deterministic exact companion and scope

Run

```bash
python -B 04-computation/periodic_polynomial_fourier_jordan_recurrence_thm3485.py
python -B -O 04-computation/periodic_polynomial_fourier_jordan_recurrence_thm3485.py
```

The standard-library companion constructs cyclotomic polynomials exactly,
detects Fourier-colour degrees by reduction at primitive roots, builds `(16)`,
and compares it with rational Berlekamp--Massey reconstruction.  It checks the
THM-3484 specialization, five structured hostiles, 35 generic packets, and
common-leading packets for every `2<=p<=9`.

For the identically zero word, the convention is `chi=D=1` and recurrence
order zero; the defect identity `(11)` is deliberately stated only for a
nonzero word.  The companion's ambient degree-zero hostile records naive
bound seven and defect seven rather than assigning a polynomial degree to
zero.

This theorem proves no LRC current, bispectrum, endpoint
factorization, Rule 30 complexity statement, tournament classification,
Jacobian statement, or LRC(14) conclusion.
