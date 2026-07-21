---
id: THM-2022
title: "Frobenius amplification of the lowest balanced Wick face proves NC2 and GMC(2)"
status: >
  PROVED. For every finite exact support in C[Z,W], a complex torus point of
  the Gaussian moment nullcone descends to an algebraic torus point. The
  lowest balanced face supplies a nonzero Laurent constant term Q by the
  one-variable Duistermaat--van der Kallen theorem. At a suitable good prime
  p, divide the moment of order p*m0 by the common factorial (p*A0)!.
  Kummer and Lucas identify the surviving residue layer with the p-fold
  dilation of the face channels, and Frobenius makes its residue exactly
  Q^p. Thus no support whose charge convex hull contains zero can be null;
  the nullcone consists exactly of the two strict one-sided charge loci.
  The earlier exposed-two-vertex gap>1 candidate reserved under this ID is
  subsumed and no longer needed.
source: codex-2026-07-21-NC2-followup
supersedes_reservation: "exposed two-vertex factorial face with gap greater than one"
depends_on:
  - THM-1630  # one-variable Duistermaat--van der Kallen constant-term theorem
  - THM-1540  # NC2 implies GMC(2)
related:
  - THM-1645
  - THM-2019
  - THM-2020
  - THM-2033
  - THM-2040
  - THM-2041
  - HYP-8800
  - HYP-8765
script: 04-computation/gmc2_frobenius_lowest_face_codex_20260721.py
output: 05-knowledge/results/gmc2_frobenius_lowest_face_codex_20260721.out
---

# THM-2022 -- Frobenius amplification of the lowest balanced Wick face

Let `Z` be a circular complex Gaussian and put `W=Zbar`, so

```text
E[Z^A W^B] = A! if A=B, and 0 otherwise.
```

For a nonzero polynomial, combine equal monomials and delete zero
coefficients:

```text
P = sum_(i=1)^k c_i Z^(a_i) W^(b_i),       c_i != 0,
q_i = a_i-b_i.
```

The theorem proves the two-dimensional nullcone conjecture NC2:

> **Theorem.** If `E[P^m]=0` for every `m>=1`, then either every `q_i>0`
> or every `q_i<0`. Conversely, either strict one-sided condition makes all
> positive moments vanish. Consequently GMC(2) holds.

The proof does not isolate first-return circuits. It preserves the whole
lowest balanced face, including every balanced-channel collision on that
face, and amplifies
its nonzero constant term by Frobenius.

## 1. Exact Wick channels

For `m>=1`, let

```text
R_m = {r in N^k : |r|=m and sum_i q_i r_i=0},
A(r) = sum_i a_i r_i = sum_i b_i r_i       (r in R_m).
```

Direct multinomial expansion and Wick balance give

```text
M_m(c):=E[P^m]
 = sum_(r in R_m) binom(m;r) A(r)! c^r,                 (1)
binom(m;r)=m!/prod_i r_i!.
```

The crucial two-dimensional feature is that a balanced channel has one
nonnegative radial exponent `A(r)` and hence one monotone Wick factorial
`A(r)!`.

## 2. Reduction to algebraic coefficients on a fixed exact support

Fix the monomial support and regard the moments in (1) as polynomials over
`Q` in the variables `c_i`. Put

```text
I=(M_1,M_2,...) subset Q[c_1,...,c_k].
```

If the complex null locus of `I` meets the coefficient torus, then the
localized algebra

```text
R=Q[c_1,...,c_k,(c_1...c_k)^(-1)]/I
```

is nonzero. (Noetherianity first replaces the infinite moment list by a
finite generating sublist.) A maximal-ideal residue field of `R` is a field
finitely generated as a `Q`-algebra, hence a finite extension of `Q` by
Zariski's lemma. Therefore a complex torus null point would imply a torus
null point

```text
c=(c_1,...,c_k) in (K^*)^k                              (2)
```

for some number field `K`. It is enough to exclude (2).

## 3. The lowest balanced face retains a nonzero toral address

Assume the support is not strictly one-sided. Equivalently,

```text
0 in conv{q_1,...,q_k}.
```

Consider the rational linear program

```text
delta = min {sum_i a_i x_i : x_i>=0, sum_i x_i=1,
                                  sum_i q_i x_i=0}.          (3)
```

Linear-programming duality supplies `lambda in Q` such that

```text
a_i-lambda q_i >= delta                                    (4)
```

for every `i`. Let `F` be the equality set in (4). An optimal
point in (3) is supported on `F`, so

```text
0 in conv{q_i:i in F}.                                     (5)
```

Define the face Laurent polynomial

```text
f_F(u)=sum_(i in F) c_i u^(q_i).                            (6)
```

There is no hidden projection collision in (6). If two points on `F` had
the same charge, (4) at equality would give the same `a`; then
`b=a-q` would also agree, contradicting distinctness of the exact monomial
support.

By (5) and the one-variable Duistermaat--van der Kallen theorem
(`THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2.md`),
the constant terms of all positive powers of `f_F` cannot vanish. Choose
`m0>=1` such that

```text
Q:=CT_u(f_F^m0) != 0.                                      (7)
```

Every multiplicity vector contributing to (7) is balanced, has length
`m0`, and is supported on `F`. For all of them, (4) gives the same integer

```text
A0=A(s)=delta*m0 in N.                                     (8)
```

This also covers endpoint cases. If the optimal face meets charge zero only
at a neutral monomial, then its nonzero coefficient already gives a nonzero
constant term (and one may take `m0=1`).

## 4. Kummer isolates the dilated face layer

Choose a rational prime `p` satisfying

```text
p>m0,                                                       (9)
```

outside the finite set of primes that ramify in `K` or at which one of the
algebraic numbers `c_i` or `Q` is not a unit. Fix a prime ideal
`pfrak` of `K` above `p`. Put

```text
M=p*m0.
```

For every `r in R_M`, (3)-(4) imply

```text
A(r)>=delta*M=p*A0,                                        (10)
```

with equality exactly when `r` is supported on `F`.

Kummer's carry formula for a multinomial coefficient says

```text
v_p(binom(M;r)) = total base-p carry count (with multiplicity). (11)
```

Because `p>m0`, the base-`p` digits of `M` are `(0,m0)`. Thus
there are no carries in (11) if and only if every units digit of every
`r_i` is zero, equivalently

```text
r=p*s,       |s|=m0.                                       (12)
```

Balance of `r` is equivalent to balance of `s`.

By (10), every channel factorial in (1) is divisible by `(p*A0)!`.
After dividing the entire moment by that integer, there are exactly three
cases modulo `pfrak`.

1. If `r` is not divisible componentwise by `p`, (11) makes the multinomial
   coefficient divisible by `p`, so the normalized term vanishes.
2. If `r=p*s` but `s` is not supported on `F`, then `A(s)>A0`.
   Both are integers, so `A(s)>=A0+1`; the factorial quotient
   `(p*A(s))!/(p*A0)!` contains the factor `p*(A0+1)` and vanishes.
3. If `r=p*s` and `s` is a balanced face channel, then
   `A(r)=p*A0`, so the factorial quotient is exactly one.

Therefore the complete residue layer of `M_M/(p*A0)!` is precisely

```text
{p*s : |s|=m0, q dot s=0, supp(s) subset F}.                (13)
```

It may contain many channels. No unique-channel assertion is used.

## 5. Lucas and Frobenius prevent cancellation inside the layer

For a channel `p*s` in (13), the multinomial Lucas congruence gives

```text
binom(p*m0;p*s_1,...,p*s_k)=binom(m0;s_1,...,s_k) mod p.     (14)
```

Writing bars for residues in the characteristic-`p` residue field, (7),
(13), and (14) yield the exact normalized identity

```text
M_(p*m0)(c)/(p*A0)!
 = sum_s binom(m0;s) c^(p*s)                    mod pfrak
 = (sum_s binom(m0;s) c^s)^p                    mod pfrak
 = Qbar^p                                       mod pfrak. (15)
```

The middle equality is Frobenius: coefficients from `F_p` are fixed by the
`p`-th power map. By the choice of `p`, `Qbar` is nonzero. Thus (15) is
nonzero, and so

```text
M_(p*m0)(c) != 0,
```

contradicting the assumed algebraic torus null point. Section 2 then excludes
every complex torus null point on every support whose charge convex hull
contains zero.

It follows that a null polynomial has all charges strictly positive or all
strictly negative. Conversely, charges add under multiplication, so a strict
one-sided polynomial has no charge-zero monomial in any positive power and
all its moments vanish. This proves NC2.

The implication recorded in
`THM-1540-the-two-dimensional-nullcone-conjecture.md` now gives GMC(2): if,
say, every charge of `P` is at
least one, then for a fixed polynomial `Q_0`, every charge of `Q_0 P^m` is
nonzero once `m` exceeds the finite negative-charge range of `Q_0`; hence
`E[Q_0P^m]=0` eventually. The negative-charge case is identical.

## 6. Why the transfer works, and why it is two-dimensional

The mechanism was found by importing two seemingly unrelated repository
ideas: retain a side channel under Frobenius twist, and attach the base-`p`
carry state to a tropical wall. Here the lower balanced face is the exact
wall, `Q` is the retained side channel, and Kummer says every non-dilated
allocation crosses a carry wall. Frobenius preserves the *whole* face sum
`Q`, so colliding primitive circuits never need to be separated.

This does not prove an analogous statement in higher Gaussian dimension.
For one complex variable, balance leaves the single scalar Wick factor
`A(r)!`, monotone in the lower-face height. With several complex variables
the Wick factor is a product of coordinate factorials; lowering one scalar
functional need not increase every coordinate valuation. That is precisely
where the above minimum-layer argument stops.

For Tournament Analysis, take balanced channels as vertices and compare
their divisibility after the normalization `(p*A0)!`, with lexicographic
order only as a tie
path. The selected quotient preserves the minimum residue layer but
forgets its residue sum. Formula (15), not transitivity or arbitrary
tie-breaking, restores that missing coordinate. The challenged assumption
is therefore explicit: noncancellation does not require a dominant channel;
an entire tied face can survive as one Frobenius power.

## 7. Scope and supersession

The proof uses arbitrary finite monomial support and arbitrary complex
coefficients. It includes arbitrary radial polynomials, arbitrary many
charges, arbitrary many primitive return circuits, neutral terms, and all
degree-resonance bands. It therefore closes the residual left by THM-2017,
THM-2018, THM-2019, HYP-8765, and HYP-8766, without invalidating their finer
effective and asymptotic statements.

The former THM-2022 candidate required an exposed two-vertex face with
factorial gap greater than one. That archimedean estimate is unnecessary:
the lowest balanced face always exists, and the carry/Frobenius argument
handles a many-vertex face with no metric gap beyond strict face separation.

The constant-term input is J.J. Duistermaat and W. van der Kallen,
"Constant terms in powers of a Laurent polynomial," *Indagationes
Mathematicae* (N.S.) 9(2) (1998), 221--231, Theorem 2 and Remark 3.
