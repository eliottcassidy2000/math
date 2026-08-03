---
id: THM-3213
title: "Tournament normalized cyclic diagonal and fast moving-jet transform"
status: >
  PROVED + VERIFIED-EXACT.  For a fixed tournament seed A of order a, the
  complete factorial path-cover row of A^{join r} is one binomial transform
  computable by a single truncated polynomial product in
  O(M(ar)+ar log r) arithmetic, improving the quadratic moving-jet evaluator.
  For the transitive seed T_r=K1^{join r}, every fixed C3-substitution
  coordinate has an exact trivariate diagonal formula.  Smooth-point analysis
  gives its full leading constant after factorial-cube normalization; dividing
  once more by binom(r-1,d-1) produces a lambda^r/r asymptotic and proves the
  twice-normalized sequence is not C-finite for every fixed d.  A prime-shift
  congruence puts a fresh p^3 in the singly normalized denominator for every
  fixed d and all but finitely many primes.  More generally, for every fixed
  q-block quotient tournament Q, F_(Q[T_r,...,T_r])(d)/(r!)^q has fresh p^q
  denominators and is not C-finite.  P-recursiveness remains open.
audit: >
  The exact transform companion checks 262 coordinates for K1, C3, and the
  same-H mask-40/mask-76 hostile.  The exact modular companion reconstructs
  the C3 quotient kernel and checks primes 5,7,11,13,17,19.  A second exact
  modular companion checks the prime-shift identity through d=12, the d=4
  constant 1944, the p-index Fermat-quotient near miss, and its finite scout;
  an independent vertex/set-partition companion reconstructs the full U_2
  path-cover profile and 1944 without the quotient kernel.  The all-depth and
  quotient-wide proofs are algebraic below.  A floating 80-digit companion
  checks the analytically proved asymptotic for d=1..4 and r=40..640 and is
  explicitly labelled NUMERICAL CONTROL ONLY.  All five
  normal/-O/stored transcripts agree and all scripts have no assert node.
source: root/frontier-synthesis-cont-2026-08-02
depends_on:
  - THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence
  - THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction
  - THM-3121-path-cover-walk-content-substitution-kernel
related:
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
fast_script: 04-computation/tournament_fast_binomial_transform_thm3213.py
fast_output: 05-knowledge/results/tournament_fast_binomial_transform_thm3213.out
fast_script_sha256: 9034eb8a23f9325378053656163aee1609c92cac3e66e21e145c2a553f8c7a9d
fast_output_sha256: c5308414097a741931418e4063d864182ddccce5084adc0b423ecd97959f79c2
prime_script: 04-computation/tournament_normalized_prime_denominator_thm3213.py
prime_output: 05-knowledge/results/tournament_normalized_prime_denominator_thm3213.out
prime_script_sha256: a8c512ad551e08329d23a0c9f36d36bf4053e1240a172895872e60755d26f4db
prime_output_sha256: e3b9c1b8de55b12296b3b7956ab57b2fc3a4ea6ca98184a62b57a72614f0f9c1
prime_shift_script: 04-computation/tournament_prime_shift_all_depth_thm3213.py
prime_shift_output: 05-knowledge/results/tournament_prime_shift_all_depth_thm3213.out
prime_shift_script_sha256: 9e6dd5129fea3c069e218b425e6376a156a5503df08b7c08572f90918aa942d1
prime_shift_output_sha256: d3622e0237ec1d1f8c09bfa2499922277050159c9164f5bceaa20e5c85b07ce2
d4_vertex_script: 04-computation/tournament_d4_prime_shift_vertex_audit_thm3213.py
d4_vertex_output: 05-knowledge/results/tournament_d4_prime_shift_vertex_audit_thm3213.out
d4_vertex_script_sha256: 49d2532ebbc9d9562fb8c28cf26faa7a94eed30519ec6d9b4ae0709d88d9d6ba
d4_vertex_output_sha256: 699d51258140dad404b6db9529d73f6c12a02d21a065c301f2a6af56236c933f
asymptotic_control_script: 04-computation/tournament_normalized_diagonal_asymptotic_control_thm3213.py
asymptotic_control_output: 05-knowledge/results/tournament_normalized_diagonal_asymptotic_control_thm3213.out
asymptotic_control_script_sha256: e298e38884709b9bdb6f18fa808413571b0682168217635f4e795a55c62cab34
asymptotic_control_output_sha256: 65b6b105a8c7343689e72bb4999d159a19c68b0db9bd927a449ab82cbfaebd90
hash_basis: LF-normalized bytes
---

# THM-3213 -- normalized cyclic diagonals and a fast moving-jet transform

**PROVED + VERIFIED-EXACT.**

## 1. A single-product transform for the complete moving jet

Let `A` be a fixed nonempty tournament of order `a`, and put

```text
B_r=A^(join r),                         N=ar.             (1)
```

For a tournament `S`, let `p_S(c)` count spanning covers by `c` unordered
directed paths and let

```text
F_S(c)=c! p_S(c),
Q_S(t)=sum_c p_S(c)(t)_(falling c).                       (2)
```

Write `F_r(c)=F_(B_r)(c)`.  THM-3181 gives

```text
Q_(B_r)(t)=Q_A(t)^r,                                     (3)
```

and falling-factorial inversion gives

```text
F_r(c)=sum_(j=0)^c (-1)^(c-j) binom(c,j) Q_A(j)^r.        (4)
```

Equation `(4)` is a binomial transform.  Dividing by `c!` turns it into one
ordinary convolution:

```text
F_r(c)/c!
 =[z^c]
   (sum_(j=0)^N Q_A(j)^r z^j/j!)
   (sum_(ell=0)^N (-z)^ell/ell!).                         (5)
```

Thus the **entire** row `(F_r(0),...,F_r(N))` is obtained by one truncated
polynomial product followed by factorial scaling.

Let `M(N)` denote the arithmetic cost of multiplying two degree-`N`
polynomials over a characteristic-zero field.  For fixed seed `A`, evaluate
the fixed-degree polynomial `Q_A` at `0,...,N` in `O(N)` arithmetic and use
binary powering for the `N+1` exponents.  Equation `(5)` then costs

```text
O(M(N)+N log r) arithmetic,              O(N) storage.    (6)
```

This is a unit-cost arithmetic bound, not a bit-complexity claim.  It improves
THM-3202's in-place `O(N^2)` forward-difference table.  After `(5)`, any fixed
`C_3` output depth costs only `O_d(N)` additional arithmetic by THM-3202's
bounded diagonal band; all depths `d<=D` cost `O_D(N)` for fixed `D`.  No
claim is made for all cyclic output depths growing with `N`.

The exact companion checks `(5)` against direct differences on `262` profile
coordinates: `A=K_1` through level 12, `A=C_3` through level 7, and two
order-five same-H hostiles through level 4.  The hostile C3-lift Hamiltonian
counts remain

```text
178036299,                         193215375,             (7)
```

so the speedup does not turn scalar `H` into a sufficient sidecar.  It
computes the complete moving jet more efficiently; it does not compress away
the data that cyclic substitution consumes.

## 2. Exact transitive-seed diagonal

Now specialize to the transitive tournament

```text
T_r=K_1^(join r),                  U_r=C_3[T_r,T_r,T_r].  (8)
```

Put

```text
u(x)=exp(x)-1,
P(X,Y,Z)=X+Y+Z+XY+XZ+YZ+3XYZ,
W(X,Y,Z)=P(X,Y,Z)/(1-XYZ).                                (9)
```

For the transitive tournament, a path cover is a set partition whose blocks
inherit their unique increasing paths.  Hence

```text
F_(T_r)(c)=c! S(r,c)=r![x^r]u(x)^c.                      (10)
```

Insert `(10)` into THM-3202's quotient-walk substitution formula.  For every
fixed `d>=1` and every `r>=1`,

```text
A_(r,d):=F_(U_r)(d)/(r!)^3
 =[x^r y^r z^r] W(u(x),u(y),u(z))^d.                     (11)
```

This is an exact closed diagonal formula.  It is not a fitted recurrence and
does not require expanding the tournaments on `3r` vertices.

## 3. Smooth diagonal asymptotic with the leading constant

Let

```text
lambda=log(2).                                            (12)
```

For every fixed `d>=1`,

```text
A_(r,d) ~ K_d r^(d-2) lambda^(-3r),                      (13)

K_d=
 9^d /
 ((d-1)! (2lambda)^(d-1) 4pi sqrt(3) lambda(1-lambda)).  (14)
```

Here is a self-contained smooth-point proof.  The singular denominator in
`(11)` is

```text
H(x,y,z)=1-u(x)u(y)u(z).                                 (15)
```

Because `u` has strictly positive coefficients,

```text
|u(x)|<=u(|x|),
```

with equality on `|x|=lambda` only at `x=lambda`.  Therefore

```text
c=(lambda,lambda,lambda)                                 (16)
```

is the unique minimal point in the positive diagonal direction:
`u(lambda)=1`, so `H(c)=0`.  It is smooth, since `H_z(c)=-2`, and the three
logarithmic gradient coordinates agree by symmetry.  The numerator at `c` is

```text
P(1,1,1)^d=9^d.                                          (17)
```

Eliminate `z` on the local singular surface:

```text
z=h(x,y)=log(1+1/(u(x)u(y))),             h(lambda,lambda)=lambda. (18)
```

Near `z=h`,

```text
H=(-hH_z)(1-z/h)+O((1-z/h)^2),            -hH_z=2lambda at c. (19)
```

Extracting the order-`d` pole in `z` contributes

```text
r^(d-1)/((d-1)!(2lambda)^d)               (20)
```

times the remaining bivariate saddle and the numerator `9^d`.

For logarithmic local variables `x=lambda exp(s)`, `y=lambda exp(t)`, the
phase is

```text
Phi(s,t)=log(h(lambda exp(s),lambda exp(t))/lambda)+s+t.  (21)
```

Its first derivatives vanish at the origin.  Direct differentiation gives

```text
-Hess Phi(0,0)
 =(1-lambda) [[2,1],[1,2]],
det(-Hess Phi)=3(1-lambda)^2.                             (22)
```

On the Cauchy torus `s,t` are imaginary, so `(22)` is the positive Gaussian
decay matrix.  The two-dimensional Laplace integral contributes

```text
1/(2pi r sqrt(3)(1-lambda)).                              (23)
```

and the exponential factor is `lambda^(-3r)`.  Multiplying `(17)`, `(20)`,
and `(23)` yields `(13)--(14)`.  Strict minimality from the equality case of
`|u(x)|<=u(|x|)` makes the complement of the local saddle exponentially
smaller.  This proves the asymptotic; the 80-digit companion is only a
numerical control of its constant.

## 4. Twice-normalized non-C-finiteness

For `r>=d`, define

```text
B_(r,d)=A_(r,d)/binom(r-1,d-1).                           (24)
```

Since `binom(r-1,d-1)~r^(d-1)/(d-1)!`, equations `(13)--(14)` give

```text
B_(r,d) ~ C_d lambda^(-3r)/r,                            (25)

C_d=9^d/((2lambda)^(d-1)4pi sqrt(3)lambda(1-lambda))>0.  (26)
```

Every complex C-finite sequence is a finite exponential polynomial
`sum_j P_j(r)alpha_j^r`, with each `P_j` an ordinary polynomial.  At the
maximal root modulus, such a sum cannot have a nonzero `alpha^r/r` leading
term: after division by that modulus, a nonzero finite unit-circle
exponential polynomial cannot tend to zero, while deleting all maximal terms
would lower the exponential growth rate.  Therefore

```text
{B_(r,d)}_(r>=d) is not C-finite for every fixed d>=1.    (27)
```

Equivalently, its ordinary generating function in `r` is not rational.  This
is an asymptotic proof, not a finite recurrence fit.

## 5. Prime-shift new-denominator obstruction at every depth

The singly normalized sequence `A_(r,d)` needs a separate argument because
its leading form `r^(d-2)lambda^(-3r)` is compatible with a repeated
characteristic root when `d>=2`.

### 5.1 The zero-shift control `d=1,2,3`

Let `p>3` be prime.  Falling-factorial inversion and Fermat's congruence give

```text
F_(T_p)(1)=1 mod p,
F_(T_p)(c)=Delta^c(t^p)|_(t=0)=Delta^c(t)|_(t=0)=0 mod p
                                                for c>=2. (28)
```

Reduce THM-3202's substitution formula modulo `p`.  Every factor-profile cell
except `(1,1,1)` vanishes, so

```text
F_(U_p)(d)=[XYZ]W(X,Y,Z)^d mod p.                         (29)
```

The coefficient in `(29)` is read directly from `P^d`:

```text
[XYZ]W^d = 3,6,6,0                 for d=1,2,3,d>=4.      (30)
```

For `d=2`, the six terms are the ordered pairs `X*YZ`, `Y*XZ`, `Z*XY`; for
`d=3`, they are the six permutations of `X*Y*Z`.  Total degree excludes
`d>=4`.  Hence for `d=1,2,3`,

```text
p does not divide F_(U_p)(d).                             (31)
```

Since `v_p(p!)=1`, the reduced denominator of

```text
A_(p,d)=F_(U_p)(d)/(p!)^3                                (32)
```

contains exactly `p^3`.  Every prime `p>3` is therefore a new denominator
prime at its own index.

A rational-valued C-finite sequence can be taken to satisfy a recurrence over
`Q`: a complex linear dependence among its rational shift columns has a
rational null vector after row reduction.  Clearing the recurrence and
initial denominators then confines every term denominator to one finite set
of primes.  Equations `(31)--(32)` contradict that confinement.  Thus

```text
{A_(r,d)} is not C-finite for d=1,2,3.                    (33)
```

For `d>=4`, `(30)` makes the **zero-shift** test silent.  The repair is to
move the prime index while keeping its residue profile fixed.

### 5.2 The all-depth prime shift

Fix any `d>=1` and put

```text
m=ceil(d/3),                    s=m-1.                    (33a)
```

For every prime `p` and every `c`, Fermat's congruence gives

```text
F_(T_(p+s))(c)
 =Delta^c(j^(p+s))|_(j=0)
 =Delta^c(j^(s+1))|_(j=0)
 =F_(T_(s+1))(c)                                  mod p. (33b)
```

The equality includes `c>s+1`: both sides there are zero modulo `p`, because
the right finite difference exceeds the degree.  THM-3202 expresses
`F_(U_r)(d)` as an integer cubic form in the complete transitive profile.
Therefore `(33b)` implies

```text
F_(U_(p+s))(d)=F_(U_m)(d)                         mod p. (33c)
```

Set

```text
K_d=F_(U_m)(d).                                           (33d)
```

This integer is strictly positive.  The tournament `U_m` has `3m>=d`
vertices; every tournament has a Hamilton path, and cutting such a path into
`d` nonempty consecutive paths gives a `d`-path cover.

Now take any prime `p>s` not dividing `K_d` and put `r=p+s`.  Then
`p<=r<2p`, so `v_p(r!)=1`, while `(33c)` makes the numerator a `p`-adic unit.
The reduced denominator of

```text
A_(p+s,d)=F_(U_(p+s))(d)/((p+s)!)^3                       (33e)
```

therefore contains exactly `p^3`.  Only finitely many primes are excluded by
`p<=s` or `p|K_d`.  Infinitely many fresh denominator primes prove

```text
{A_(r,d)} is not C-finite for every fixed d>=1.           (33f)
```

For `d=4`, the cheapest shift is `s=1` and

```text
K_4=F_(U_2)(4)=1944=2^3 3^5.                             (33g)
```

Thus **every** prime `p>3` contributes a new `p^3` denominator at `r=p+1`.
For `p>5`, the first layer at the unshifted index `r=p` is instead the
Fermat-quotient combination

```text
12(24q_p(2)-21q_p(3)+5q_p(5)) mod p.                     (33h)
```

Here `q_p(a)=(a^(p-1)-1)/p mod p` for `p` not dividing `a`.  The combination
vanishes at `p=13`; exact controls give
`v_5(F_(U_5)(4))=3` and `v_13(F_(U_13)(4))=2`.  The finite scout finds no
other zero first layer through `10^6`, but this is only finite evidence.  The
unshifted test is therefore the corrected near miss, not the uniform proof.

### 5.3 Quotient-wide form

The same argument is not special to `C_3`.  Fix any tournament `Q` on `q>=1`
vertices and put

```text
V_r=Q[T_r,...,T_r].                                      (33i)
```

THM-3121 writes `F_(V_r)(d)` as an integer polynomial in `q` copies of the
complete factor profile.  Hence `(33b)` gives

```text
F_(V_(p+s))(d)=F_(V_(s+1))(d)                     mod p. (33j)
```

Choose `m=ceil(d/q)`, `s=m-1`, and
`K_(Q,d)=F_(V_m)(d)>0`; positivity again follows by cutting a Hamilton path
on the `qm>=d` vertices.  For every prime `p>s` not dividing `K_(Q,d)`, the
reduced denominator of

```text
F_(V_(p+s))(d)/((p+s)!)^q                                (33k)
```

contains exactly `p^q`.  Thus every fixed quotient and fixed positive depth
has a factorial-normalized sequence which is not C-finite.  This
quotient-wide statement concerns recurrence class only; it supplies no
uniform algorithm for growing quotient or output depth.

## 6. Recurrence taxonomy and operational boundary

The three normalizations

```text
F_(U_r)(d),             A_(r,d)=F_(U_r)(d)/(r!)^3,
B_(r,d)=A_(r,d)/binom(r-1,d-1)                            (34)
```

differ by nonzero hypergeometric terms for fixed `d`.  Termwise
multiplication by a hypergeometric sequence and by its reciprocal preserves
P-recursiveness.  Consequently the three sequences in `(34)` are
P-recursive either all together or not at all.  Their P-recursive status is
**OPEN**.

By contrast, their C-finite behavior differs sharply:

| sequence | fixed-`d` status |
|---|---|
| raw `F_(U_r)(d)` | not C-finite for every `d`, by THM-3202's factorial growth |
| `A_(r,d)` | not C-finite for every `d`, by the prime shift `(33a)--(33f)` |
| `B_(r,d)` | not C-finite for every `d`, by `(25)` |

This separates three notions that finite sequence fitting tends to conflate:
fast exact evaluation, a fixed constant-coefficient closed form, and a
polynomial-coefficient recurrence.  Equation `(5)` supplies the first;
equations `(27),(33f)` refute the second at every fixed depth; the third
remains open.

## Exact reproduction

Run

```text
python3 04-computation/tournament_fast_binomial_transform_thm3213.py
python3 -O 04-computation/tournament_fast_binomial_transform_thm3213.py
python3 04-computation/tournament_normalized_prime_denominator_thm3213.py
python3 -O 04-computation/tournament_normalized_prime_denominator_thm3213.py
python3 04-computation/tournament_prime_shift_all_depth_thm3213.py
python3 -O 04-computation/tournament_prime_shift_all_depth_thm3213.py
python3 04-computation/tournament_d4_prime_shift_vertex_audit_thm3213.py
python3 -O 04-computation/tournament_d4_prime_shift_vertex_audit_thm3213.py
python3 04-computation/tournament_normalized_diagonal_asymptotic_control_thm3213.py
python3 -O 04-computation/tournament_normalized_diagonal_asymptotic_control_thm3213.py
```

Each pair must reproduce its declared output byte for byte.  The first four
companions use exact integer/rational arithmetic.  The fifth is explicitly a
floating numerical control and is not a proof dependency for `(13)`.  All
checks use explicit exceptions, so optimized mode retains them.  QED.
