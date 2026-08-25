---
id: THM-4034
title: "Exceptional-quartic global conductor of degree 178"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the
  THM-3703 exceptional-quartic restriction algebra S inside its normalization
  K[x], the global conductor is cK[x], where c is the unique monic
  degree-178 conductor polynomial.  It factors exactly as
  x^2(x^2-1)^2 h_172 with h_172 squarefree and coprime to the retained
  factor.  Intrinsically, c is the monic gcd of the three divided-difference
  resultants for (B,C), (B,E), and (C,E).  The numerical-semigroup conductor
  170 is not this ideal conductor: x^170 is not even in S.  No singularity
  classification, higher stable lift, Keller map, or JC(2) claim is made.
source: jc2-double-zero-rebuild-20260824 / conductor reconstruction, 2026-08-24
audit: >
  PASS -- an independent quotient-matrix reconstruction recovered the sharp
  178-column lower boundary and the same monic kernel at the good split
  fibre.  A separate divided-difference calculation recovered resultant
  degrees 578,488,338, the degree-178 triple gcd, the hostile degree-186
  two-resultant gcd in that reduction, the retained multiplicities, and the
  x^170 failure.  The exact PID argument, trace-dual resultant lemma,
  good-reduction directions, field typing, hashes, and strict scope were
  hostile-audited.  Final normal and optimized executions byte-match the
  frozen transcript.
depends_on:
  - THM-3703-russell-cylinder-exceptional-quartic-sagbi-module
related:
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3688-russell-cylinder-exceptional-quartic-actual-j1-j2-lift
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
  - THM-4039-exceptional-quartic-j3-lift-and-adjacent-gate-rigidity
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.out
script_sha256: f44baddb7d7c8c4d204d19952dadabfa4a52235e667295555481591a4f66ea11
output_sha256: 6854db44bc4959f4b6b1edc298b066e060533de426e5840957b26000cf1afd8b
hash_basis: raw LF bytes
---

# THM-4034 -- the global conductor begins in degree 178

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The numerical
semigroup boundary `170` in THM-3703 does not identify the conductor ideal of
the actual restriction algebra.  The latter has a unique monic generator of
degree `178`, and all three coordinate pairs are needed for its intrinsic
divided-difference description.

All rings below have characteristic zero.

## 1. The conductor condition as a finite exact kernel

Retain the irreducible quartic field and restriction algebra of THM-3703,

```text
K=Q[alpha]/(F_6(alpha)),
S=K[B,C,E] subset K[x],
S=direct_sum_(r=0)^17 K[Z]p_r,             deg Z=18.    (1)
```

THM-3703 proves that `K[x]` is the normalization of `S` and that the `89`
leading-degree gaps of `S` end at `169`.  Its monic Apéry basis gives an
exact descending normal form

```text
rho:K[x] -> span_K{x^g:g is one of the 89 gaps},
ker rho=S.                                                (2)
```

Because `Z` is monic of degree `18`, ordinary division gives

```text
K[x]=direct_sum_(j=0)^17 K[Z]x^j.                        (3)
```

Let

```text
mathfrak c=(S:K[x])={f in K[x]:fK[x] subset S}.           (4)
```

Equations `(2)` and `(3)` turn membership in `(4)` into the finite condition

```text
f in mathfrak c
 iff rho(x^j f)=0 for every 0<=j<18.                      (5)
```

For polynomials of degree at most `177`, `(5)` is a `K`-matrix with

```text
18*89=1602 rows,                    178 columns.          (6)
```

The companion builds every entry over `K` and reduces it at the good split
place

```text
p=137,                         alpha=44 mod 137.          (7)
```

All denominators and the leading coefficient of `F_6` are units there.  The
selected `178`-square has

```text
det=1 mod 137,
selected-row hash=
be7e4888affc2676b671a0791a05794b2b7229ec997a2312247210a3eb50d145. (8)
```

This is the sound positive-minor direction: the exact square over `K` is
nonzero, so `(6)` has full column rank in characteristic zero.  Consequently

```text
mathfrak c intersect K[x]_(<=177)=0.                      (9)
```

No modular rank failure is used.

## 2. Exact reconstruction and the principal conductor

Write a monic candidate as

```text
c=x^178+sum_(i=0)^177 u_i x^i.                           (10)
```

The square selected in `(8)` determines the `178` unknown field coefficients
uniquely.  Expanding in the ordered power basis
`(1,alpha,alpha^2,alpha^3)` gives a `712`-square over `Q`, which the companion
solves exactly.  It then discards the selected-square shortcut and checks all
`1602` equations in `(5)` for the reconstructed `c`.  Equivalently, it checks

```text
c,xc,...,x^17c in S.                                    (11)
```

The canonical ascending-degree, ascending-power-basis coefficient hash is

```text
sha256(c)=8e19b56da62b319108eacd59884e440ac902d1a886e6fafa0ac7d9fb1573da99. (12)
```

Equations `(3)` and `(11)` imply `cK[x] subset S`, so `c` lies in the
conductor.  Conversely, `(4)` is an ideal of the PID `K[x]`; write it as
`dK[x]` with `d` monic.  Since `c` lies in this ideal, `d` divides `c`.
Equation `(9)` gives `deg d>=178`, while `(10)` gives `deg c=178`.  Hence

```text
boxed: (S:K[x])=cK[x],                    deg c=178.      (13)
```

In particular

```text
dim_K K[x]/cK[x]=178,
dim_K S/cK[x]=178-dim_K K[x]/S=178-89=89.                (14)
```

The equality `178=2*89` is recorded here as exact conductor arithmetic; no
extra local singularity or Gorenstein classification is inferred from it.

## 3. Retained factor and the false degree-170 shortcut

Exact division gives

```text
c=x^2(x^2-1)^2 h_172,                    deg h_172=172,  (15)
```

where `h_172(-1)h_172(0)h_172(1)` is nonzero.  Its coefficient hash is

```text
sha256(h_172)=939fa31b2b9efcda5c935873ee782cc5139f1dac3c21b68fe8dff1c6ed04c54b. (16)
```

At the good place `(7)`, the reduction of `h_172` is squarefree and coprime
to `x(x^2-1)`.  The corresponding resultants are therefore nonzero in
characteristic zero, proving the same squarefree and coprime statements over
`K`.  Thus the three retained normalization points occur in the conductor
divisor with exact multiplicity two, and all other roots of `c` are simple.
This multiplicity statement alone does not classify the image singularities.

The numerical-semigroup conductor `170` from THM-3703 means that every
leading degree at least `170` occurs in `S`.  It does not say that every
monomial of that degree belongs to `S`.  The exact hostile normal form gives

```text
rho(x^170) !=0,                         so x^170 notin S. (17)
```

Therefore the tempting formula `mathfrak c=x^170K[x]` is false at its first
generator, not merely unproved.

## 4. Intrinsic three-resultant formula

For polynomials `f,g in K[x]`, define

```text
Delta_f(x,y)=(f(x)-f(y))/(x-y),
R_(f,g)(x)=Res_y(Delta_f(x,y),Delta_g(x,y)).              (18)
```

Use `(18)` for the three pairs `(B,C)`, `(B,E)`, and `(C,E)`.  Their nonzero
good reductions prove that each pair is primitive: a nonprimitive pair would
leave a common non-diagonal divided-difference factor and make its resultant
zero.  Thus `K(f,g)=K(x)`.  Here is the conductor mechanism.  Scaling `f` and
`g` by units of `K` is harmless, so take them monic.  Put `R=K[f]`,
`A_0=R[g]`, and `N=K[x]`.  These are finite free `R`-modules of the same
rank.  If `H` is the monic minimal polynomial of `g` over `K(f)`, their trace
duals are

```text
N^*=(1/f')N,                         A_0^*=(1/H'(g))A_0. (19)
```

Restriction of trace functionals gives `N^* subset A_0^*`; multiplying `(19)`
through gives

```text
(H'(g)/f')N subset A_0.                                  (20)
```

The norm identity for divided differences identifies `H'(g)/f'`, up to a
nonzero element of `K`, with `R_(f,g)`.  Thus every resultant in `(18)` lies
in `(A_0:N)`, hence in `(S:K[x])=cK[x]`.  Therefore `c` divides all three.

For degrees `(30,21)`, `(30,18)`, and `(21,18)`, the leading-form bounds at
infinity are respectively

```text
(30-1)(21-1)-(gcd(30,21)-1)=578,
(30-1)(18-1)-(gcd(30,18)-1)=488,
(21-1)(18-1)-(gcd(21,18)-1)=338.                         (21)
```

The good reduction `(7)` attains all three bounds, so the exact leading
coefficients are units at that place.  The reduced monic triple gcd has
degree `178` and equals the reduction of `(12)`.  Hence the monic exact triple
gcd is integral at this place and reduces to a common divisor of those three
polynomials, so it has degree at most `178`.  Since it is divisible by the
degree-`178` polynomial `c`, equality follows:

```text
boxed: c=monic gcd(R_(B,C),R_(B,E),R_(C,E)).             (22)
```

All three resultants are load-bearing in this proof.  As a hostile control,
the first two alone have gcd degree `186` in the fibre `(7)`, with extra
factor

```text
x^2(x^2-1)^3.                                           (23)
```

Equation `(23)` is explicitly a mod-`137` hostile statement, not an exact
two-resultant factorization claim.

## 5. Scope and reproduction

This theorem computes the global conductor ideal and its divisor in the
normalization.  It does not classify the singularities of `Spec S`, prove a
complete-intersection or Gorenstein presentation, change the stagewise
Jacobian-image theorem, construct a higher stable lift or global pair, or
prove a Keller map or counterexample to `JC(2)`.

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py
```

Normal and optimized exact runs byte-match the stored LF transcript.  The
script contains no Python `assert` statements.  **QED.**
