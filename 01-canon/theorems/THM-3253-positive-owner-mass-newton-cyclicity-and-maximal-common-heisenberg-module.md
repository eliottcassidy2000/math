---
id: THM-3253
title: "Positive owner-mass all-gauge cyclicity and maximal common Heisenberg module"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every integer dilation in THM-3246's (3,5;1,2) lane, compactify its
  168 strictly positive cell masses by zero on the THM-3234 Singer plane.
  The resulting 13 by 13 matrix is nonsingular in all 8,064 primitive Singer
  gauges.  Exact degree-26 Newton certificates handle the natural 672-gauge
  family; exact Z-polynomial factorization plus finite-field root-free
  certificates handle all twelve multiplier classes.  A separate shifted
  Newton certificate independently proves all 84 reduced classes.  Placing
  any mass matrix on one central slice gives a nonnegative packet whose H_13
  orbit spans the sharp 2,041-dimensional common submodule from THM-3250.
  The plane relocation remains abstract and is not a canonical endpoint current.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins THM-3234, promoted THM-3246 and
  promoted THM-3250; rebuilds all 168 numerator quadratics and their
  positivity; verifies the deterministic Singer, scalar, anti-diagonal,
  Frobenius and reflection reductions; computes exact Bareiss determinants;
  and checks every finite head and all 27 Newton coefficients in each tail.
  It additionally reconstructs 84 exact integer determinant polynomials,
  certifies their exact powers of g, compares them with independent
  finite-field interpolation, and enumerates a root-free prime certificate
  for every residue.  It tests the 672 Newton gauges at four independent
  dilations and verifies the charged/neutral dimension arithmetic.  Normal,
  optimized and stored transcript replay and the LF hashes are required.
  Independent reconstruction additionally checked all 8,064 gauges at g=1,
  all symmetry identities, the exceptional exact-order lane, the finite-field
  certificates, normal and optimized replay, and the declared LF hashes.
  A second independently developed companion computes every reduced
  determinant by fraction-free elimination, checks a 57-step exact head and
  one-sign Newton tails from the sharp common base 58, and matches 64,512
  two-prime hostile determinants plus separate domain-determinant controls.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate
related:
  - THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity
  - THM-3247-heisenberg-central-fourier-decomposition-and-canonical-current-cyclicity
script: 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
output: 05-knowledge/results/lrc_positive_owner_mass_newton_cyclicity_thm3253.out
script_sha256: 89aa2a399848ae52e8dd18de9967c7ea2940c04521434ad99407f7be96bdd700
output_sha256: a96010c22126d391bf490f8535dcb3b93809f63e8b705fd23c6223a962bdae58
secondary_script: 04-computation/lrc_positive_owner_mass_uniform_newton_audit_thm3253.py
secondary_output: 05-knowledge/results/lrc_positive_owner_mass_uniform_newton_audit_thm3253.out
secondary_script_sha256: a57b12a6f4ccd2cc35007896f567656b11ec95748bda408380373730874e2719
secondary_output_sha256: 908dd123de44e169cca80f4349e6c976937e129ba5925c0da36f83f72e3eeb87
hash_basis: LF-normalized bytes
---

# THM-3253 -- positive owner-mass all-gauge cyclicity and maximal common Heisenberg module

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3252 shows that the signed second-corrector word clears THM-3250's
charged determinant gate.  The determinant is not merely an asymptotic
shadow.  The actual positive overlap masses from THM-3246 are already cyclic
at every integer dilation in every primitive Singer gauge.

The resulting central-slice packet is nonnegative.  It does not give a
positive intertwiner to the blowup carrier: instead, its orbit realizes
exactly the largest source submodule which can occur on both sides of the
neutral-character mismatch.

## 1. Positive all-dilation owner masses

In THM-3246's lane

```text
(P,Q;e,f)=(3,5;1,2),
D_g=(504g-1)(840g-2),                                  (1)
```

write

```text
N_(g,j)=D_g I_(g,j),             0<=j<168.             (2)
```

For every integer `g>=1`, THM-3246 gives

| cells `j` | `N_(g,j)` |
|---|---|
| `0<=j<=5` or `162<=j<=167` | `12096g^2-1032g+2` |
| `6<=j<=23` or `144<=j<=161` | `12096g^2-24g` |
| `24<=j<=71` | `(16044-168j)g^2+48g` |
| `72<=j<=95` | `4032g^2+96g` |
| `96<=j<=143` | `(168j-12012)g^2+48g` |

Every displayed quadratic is strictly positive for `g>=1`: its value at
one and its derivative on `[1,infinity)` are positive.  Also

```text
N_(g,167-j)=N_(g,j).                                   (3)
```

Since `D_g>0`, the 168 masses `I_(g,j)` are strictly positive.

## 2. All primitive Singer gauges

Use THM-3234's deterministic plane

```text
F_169=F_13[u]/(u^2-2),          alpha=1+2u.             (4)
```

For every primitive multiplier and every phase,

```text
a in (Z/168Z)^*,
b in Z/168Z,                                              (5)
```

define a mass matrix on `F_13^2` by

```text
A_g^(a,b)(0,0)=0,
A_g^(a,b)(alpha^(b+aj))=I_(g,j).                        (6)
```

Thus `(6)` is nonnegative, is strictly positive on the punctured plane, and
has only the abstract completion point as a zero.  The theorem is

```text
det A_g^(a,b) != 0                                      (7)
```

for every integer `g>=1` and all

```text
phi(168)*168=48*168=8064                                (8)
```

gauges in `(5)`.

## 3. Phase and generator reductions

Let `B_g^(a,b)` be the cleared integer matrix obtained from `(6)` with
`I_(g,j)` replaced by `N_(g,j)`.  Put

```text
F_(a,b)(g)=det B_g^(a,b).                               (9)
```

Three exact symmetries reduce the 8,064 cases to 84 determinant polynomials.

First,

```text
alpha^14=6 in F_13^*.                                   (10)
```

Changing `b` by 14 therefore applies the same scalar permutation to rows
and columns.  Its two determinant signs square, so

```text
F_(a,b+14)(g)=F_(a,b)(g).                               (11)
```

Second,

```text
alpha^7=9u,
(x,y) |-> alpha^7(x+yu)=(5y,9x).                       (12)
```

This is transpose followed by multiplication permutations on the two axes.
Their signs are the Legendre symbols of 5 and 9, namely `-1` and `+1`.
Hence

```text
F_(a,b+7)(g)=-F_(a,b)(g).                              (13)
```

It remains to handle `b=0,...,6` for every generator class.

Finally, Frobenius is

```text
(x,y)->(x,-y),                  alpha^j->alpha^(13j).   (14)
```

It is a column permutation, so it identifies multiplier `a` with multiplier
`13a`, with the phase correspondingly multiplied by 13.  Reflection `(3)`
gives

```text
A_g^(-a,b)=A_g^(a,b+a)                                  (15)
```

after relabelling `j` by `167-j`.  Thus the 48 unit multipliers split into
twelve free orbits under

```text
H={1,13,155,167}={+/-1,+/-13}.                         (16)
```

One may take the representatives

```text
1,5,11,17,19,23,29,31,43,47,59,71.                    (17)
```

Together with `b=0,...,6`, this leaves `12*7=84` cases.

## 4. Exact Newton certificates for the natural gauges

Abbreviate `F_b=F_(1,b)`.  Every entry of `B_g^(1,b)` is quadratic in `g`,
so

```text
degree F_b <= 26.                                       (18)
```

For a polynomial of degree at most 26 and an integer base `r`, Newton's
identity is

```text
F_b(r+n)=sum_(k=0)^26 Delta^k F_b(r) binom(n,k)          (19)
```

for every integer `n>=0`.  The binomial coefficients are nonnegative.

Exact fraction-free determinants give:

1. For every `b=0,1,2,3,4`, all 27 numbers

   ```text
   Delta^k F_b(1),       0<=k<=26,                      (20)
   ```

   are strictly negative.  Their joint digest is

   ```text
   1aeb6d4070908447584ff0fee52c2ab7e7f0d2287f766bf69b962ef6f2815e16. (21)
   ```

   Hence `F_b(g)<0` for every integer `g>=1`.

2. For `b=5`, direct evaluation gives

   ```text
   F_5(g)<0,                 1<=g<=17.                  (22)
   ```

   All 27 coefficients `Delta^k F_5(18)` are strictly positive, so

   ```text
   F_5(g)>0,                 g>=18.                     (23)
   ```

   The exact head and tail digests are

   ```text
   b228dcdcca8a65cc11b61894e5cf07921238596ff03680873d88fe427d444274,
   f7a36db8a325e9bc576343143455035dab7e8bc1b16341159406c7427b0d52af. (24)
   ```

3. For `b=6`,

   ```text
   F_6(1)<0,                                               (25)
   ```

   while all 27 coefficients `Delta^k F_6(2)` are strictly positive.
   Their digest is

   ```text
   5fb6da1e11014bec9e2e97daa588f982ff318d6a437ebfc08da15141c2da40b4. (26)
   ```

   Thus `F_6(g)>0` for every `g>=2`.

Equations `(20)--(26)` prove `F_b(g)!=0` for the seven representatives.
The symmetries in Section 3 prove `(7)` for the natural 672 gauges with
`a in H`.  This sign-refined subtheorem is independently hostile-audited.
Since

```text
det A_g^(a,b)=det B_g^(a,b)/D_g^13,                     (27)
```

the positive rational mass matrices have the same nonvanishing.

The sign changes between `g=17,18` in phase 5 and between `g=1,2` in phase
6 explain why a single base-one one-sign Newton certificate would be false.
The shifted tails are load-bearing, not cosmetic.

## 5. Exact finite-field certificates for all generators

For each of the 84 pairs in `(17)` with `0<=b<=6`, exact interpolation from
the 27 Bareiss values at `g=0,...,26` gives

```text
F_(a,b)(g)=sum_(i=0)^26 c_i g^i in Z[g].               (28)
```

The conversion from the Newton basis to the monomial basis is performed over
`Q` and every denominator is checked to be one.  If

```text
r_(a,b)=min{i:c_i!=0},
Q_(a,b)(g)=F_(a,b)(g)/g^r_(a,b),                       (29)
```

then `Q_(a,b)` lies in `Z[g]`, has nonzero constant term, and the exact order
census is

```text
r=4:2,  r=5:18,  r=6:47,  r=7:15,  r=8:2.             (30)
```

For every pair there is a prime

```text
p in {29,31,37,41,43,47,53,59,61}                     (31)
```

for which the reduction of `Q_(a,b)` has no root in `F_p`.  Enumerating all
residues gives the first-certificate census

```text
p=29:34, 31:21, 37:9, 41:5, 43:7,
  47:4, 53:2, 59:1, 61:1.                              (32)
```

The complete 84-row `(a,b,p,r)` table has digest

```text
34771ce3677e1fa9cb1932324b54aa743e6e9b97304805e65459071d3b9fc132. (33)
```

The companion independently interpolates the same polynomial directly in
`F_p` and compares every coefficient with the reduction of `(28)`.  It also
checks `(28)` at `g=80`, outside the interpolation window.

Suppose `F_(a,b)(g)=0` for an integer `g>=1`.  Since `g^r` is nonzero as an
integer, `(29)` gives `Q_(a,b)(g)=0`.  Reduction modulo its certificate prime
contradicts `(31)`.  This includes `p|g`, because residue zero was enumerated
and `Q_(a,b)(0)!=0 mod p`.  Thus all 84 representatives are nonsingular, and
Section 3 proves `(7)` for all 8,064 gauges.

The exact-factor step is load-bearing.  In the case `(a,b)=(17,3)`, reduction
modulo 31 makes the apparent order at zero rise from seven to eight because
31 divides the exact seventh coefficient.  Factoring only a modular order
would be invalid; the exact quotient instead receives root-free prime 59.

## 6. Independent uniform Newton proof for all gauges

The secondary companion proves the all-gauge theorem by a disjoint exact
route.  Use the twelve multiplier representatives in (17) and phases
b=0,...,6.  For each of the resulting 84 classes, put

~~~text
P_(a,b)(g)=det B_g^(a,b).                               (33a)
~~~

Every P_(a,b) has degree exactly 26.  Exact Bareiss elimination gives:

~~~text
P_(a,b)(g)!=0,                         1<=g<=57;         (33b)

Delta^k P_(a,b)(58), 0<=k<=26,
are all nonzero and have one common sign for fixed (a,b).                  (33c)
~~~

The first line comprises 4,788 exact head determinants.  The second comprises
2,268 exact Newton coefficients, with 44 positive-sign and 40 negative-sign
classes.  Their ordered digests are

~~~text
head:
97260646b6d140268a649f6f5d4e3a8c0d6af6849945e81b172a2170e310a63f,

tail:
31fe23e30c06e0ddde4bb19f19c5d44e5795f00d0456e6a5e71c846d88a9c008. (33d)
~~~

Newton's identity (19) makes every summand one-signed for g>=58, so
(33b)--(33c) independently prove nonvanishing for every positive integer g.
The head is load-bearing: sixteen classes change sign at adjacent positive
integers, and

~~~text
P_(43,1)(57) P_(43,1)(58)<0.                            (33e)
~~~

Thus no common one-sign Newton tail can start before 58.

As an indexing and elimination hostile, the companion also checks all 8,064
gauges at g in {1,57,58,1000} modulo each of 1000000007 and 1000000009.
All 64,512 determinants are nonzero.  A separate SymPy domain determinant
matches 18 boundary, identity, and sign-transition controls.  This uniform
Newton proof and Section
5's exact-factor/root-free proof share the mass table and Singer model but
not their all-dilation nonvanishing certificate.

## 7. A nonnegative packet attaining the common-module ceiling

Fix any case in `(7)` and abbreviate its mass matrix by `A`.  Let `K` be
THM-3250's splitting field.  On its exact-address carrier define the rational packet

```text
W_A=sum_(s,t in F_13) A_(s,t)[s,t,0].                  (34)
```

It is nonnegative.  In the central Fourier basis,

```text
[s,t,0]=(1/13)sum_(kappa in F_13) E_(s,t)^kappa.        (35)
```

Hence every charged block of `(34)` has coefficient matrix `A/13`.
By `(7)` and THM-3250 it is cyclic and spans all 169 dimensions.  The twelve
central idempotents lie in `K[H_13]`, so all block projections belong to the
orbit span.  The charged blocks therefore contribute

```text
12*169=2028.                                            (36)
```

In the neutral block `H_13` translates only the `s` coordinate and leaves
`t` as multiplicity.  For a Fourier row vector `ell_a` in `s`, invertibility
of `A` implies

```text
ell_a A != 0                         for every a.        (37)
```

Thus `(34)` has a nonzero component in each of the thirteen common neutral
characters `chi_(a,0)`.  An orbit in a one-dimensional isotypic character
contributes only one dimension, so its neutral orbit span is exactly 13.
Combining `(36)` and `(37)`,

```text
dim_K span(H_13.W_A)=2028+13=2041.                     (38)
```

THM-3250 proves that 2041 is the maximum rank of any equivariant map from
the exact-address carrier to the regular nonvertical blowup carrier.  The
module in `(38)` has all charged blocks plus one copy of every common neutral
character, so it is precisely a maximal common submodule.  The nonnegative
packet uses the neutral overlap; it does not contradict or eliminate the
remaining 156-dimensional neutral mismatch.

## 8. Scope

The entries in `(6)` are genuine positive THM-3246 cell masses, but placing
them at Singer-plane points and on the central slice `delta=0` is an abstract
relocation.  No physical LRC owner-to-plane map, canonical endpoint packet,
Boolean observable, positive equivariant map, Markov clutch, or compatibility
with the full affine/Singer action is constructed.  In particular `(34)` is
not the canonical current treated in THM-3247.

The theorem treats one ordered `(3,5;1,2)` lane and all 8,064 primitive
Singer gauges.
It proves no row exclusion, arbitrary-radial NC2 theorem, Gaussian Moment
Conjecture, or `LRC(14)` decrement.  Its exact gain is that positivity of the
owner weights and maximal charged/neutral cyclic rank coexist at every
dilation.  The remaining obstruction is realization and equivariant
transport, not matrix rank.

## 9. Exact companions

Run

```text
python3 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
python3 -O 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
python3 04-computation/lrc_positive_owner_mass_uniform_newton_audit_thm3253.py
python3 -O 04-computation/lrc_positive_owner_mass_uniform_newton_audit_thm3253.py
```

and compare LF-normalized bytes with the two declared outputs.  The
companions use exact integer and finite-field arithmetic only, with no
floating point, randomness, discovery cache, or optimization-sensitive
assertions.

QED.
