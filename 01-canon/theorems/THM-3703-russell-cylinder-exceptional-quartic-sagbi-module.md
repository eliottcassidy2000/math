---
id: THM-3703
title: "Russell-cylinder exceptional-quartic SAGBI and Apéry module"
status: >
  PROVED + VERIFIED-EXACT.  Over the irreducible quartic field of THM-3683,
  the exceptional-fold restriction algebra has a six-element monic SAGBI
  basis of degrees 18,21,30,71,83,124 and an explicit free rank-18 module
  presentation over its degree-18 generator.  Its leading-degree semigroup
  has genus 89 and conductor 170, so the filtered restriction space through
  degree 375 has dimension 287.  The 287 raw target restrictions used by the
  THM-3687/structured coupled selector form an exact K-basis of that filtered
  space; split fibres supply a nonzero-minor certificate, not an inference
  from modular failure.  The global conductor ideal is not computed here.
source: jc_zero_debt_lift / exceptional-quartic restriction-ring basis, 2026-08-22
depends_on:
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
related:
  - THM-3687-russell-cylinder-exceptional-quartic-actual-j0-lift
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py
auxiliary_script: 04-computation/jc2_russell_cylinder_exceptional_quartic_restriction_basis_modular_probe.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.out
auxiliary_output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_restriction_basis_modular_probe.out
script_sha256: 3022ce723c58af0c1b59283ac8d2eb7629fbd74a63eeb139f43aa0e1612d67d9
output_sha256: 41fdeabedb9113255bffa135b124be6d53df2d8acb1a8a977fd5066831407f3b
hash_basis: raw LF bytes
---

# THM-3703 -- the exceptional restriction algebra has a six-degree grammar

**PROVED + VERIFIED-EXACT.**  The four sixth-debt roots do not merely share
the same successful raw monomial counts.  They are four embeddings of one
restriction algebra with an exact finite SAGBI and Apéry-module description.
This replaces the 1,017-column cutoff packet by a canonical 18-residue
object and explains the recurrent rank `287` structurally.

All rings below have characteristic zero.

## 1. Field and restriction algebra

Use the irreducible polynomial and field from THM-3683/3687,

```text
F_6(T)=72783360T^4-77822208T^3-28419741T^2
                         +7849770T-1276420,
K=Q[alpha]/(F_6(alpha)).                                (1)
```

Let `Q_alpha` be THM-3683's exceptional fold and restrict the Russell
cylinder compiler

```text
D=1+x^2Q_alpha,
B=(D-1)(D+2)^2,       C=xD(D+2),       E=Q_alpha(D+3),
C^2E=B(B+4).                                            (2)
```

Thus

```text
S=K[B,C,E] subset K[x],        (deg B,deg C,deg E)=(30,21,18). (3)
```

At the retained triple `x=-1,0,1`, one has `(B,C,E)=(0,0,-3)`.
For a nonzero polynomial over `K`, write `monic(H)=H/lc(H)`.  Put

```text
Z =monic(E+3),             C_0=monic(C),       B_0=monic(B),

G_71 =monic(Z^4-C_0^2B_0),
G_83 =monic(C_0^4-Z^3B_0),
G_124=monic(G_71 Z^3-G_83 C_0^2).                     (4)
```

The subtractions in `(4)` are not guessed degree labels: in each line the
two summands are monic with the same leading degree.  Exact cancellation
gives

```text
deg(Z,C_0,B_0,G_71,G_83,G_124)=(18,21,30,71,83,124). (5)
```

Every object in `(4)` is the restriction of an actual target element: only
target-ring addition, multiplication, a constant shift, and multiplication
by units of `K` were used.

## 2. The free rank-18 Apéry module

In residue order modulo `18`, define

```text
p_0 =1,                       p_9 =C_0^3,
p_1 =C_0G_124,                p_10=G_71G_83,
p_2 =C_0G_71,                 p_11=G_83,
p_3 =C_0,                     p_12=B_0,
p_4 =G_83^2,                  p_13=C_0G_71G_83,
p_5 =B_0G_83,                 p_14=C_0G_83,
p_6 =C_0^2,                   p_15=C_0B_0,
p_7 =C_0G_83^2,               p_16=G_124,
p_8 =C_0B_0G_83,              p_17=G_71.              (6)
```

These elements are monic and their degrees are

```text
(a_0,...,a_17)=
(0,145,92,21,166,113,42,187,134,
 63,154,83,30,175,104,51,124,71),                    (7)
```

with `a_r congruent r (mod 18)`.  The exact identity is

```text
S = direct_sum_(r=0)^17 K[Z] p_r.                     (8)
```

Here is a short proof whose finite part is replayed literally.  Let the
right side of `(8)` be `M`.  Distinct summands have distinct leading degrees
modulo `18`, so the sum is direct.  Every `p_r` lies in `S`, hence `M subset
S`.  Multiplication by `Z` preserves `M` tautologically.  The companion
reduces all

```text
3*18=54 products       Zp_r, C_0p_r, B_0p_r            (9)
```

to zero by descending leading-degree reduction against `(6)`.  The largest
reduction takes `120` steps and all reductions together take `522`; every
field coefficient and step is hashed.  Thus `M` contains `1` and is stable
under the generators `Z,C_0,B_0` of `S`, giving `S subset M`.  This proves
`(8)`.

## 3. SAGBI, gaps, and normalization

Equation `(8)` makes the leading-degree semigroup exact:

```text
Gamma=<18,21,30,71,83,124>.                            (10)
```

Its Apéry set modulo `18` is exactly `(7)`.  The `89` gaps are

```text
1..17,
19,20,22..29,31..35,37,38,40,41,43..47,49,50,52,53,
55,56,58,59,61,62,64,65,67,68,70,73,74,76,77,79,80,
82,85,86,88,91,94,95,97,98,100,103,106,109,112,115,
116,118,121,127,130,133,136,139,148,151,157,169.      (11)
```

In particular,

```text
genus(Gamma)=89,          conductor(Gamma)=170.        (12)
```

The word *conductor* in `(12)` means the numerical-semigroup conductor:
every degree at least `170` occurs and `169` does not.  It is not a claim
that the global conductor ideal `(S:K[x])` is `x^170K[x]` or has degree
`170`; that ideal remains open.

Because `(8)` gives a monic filtered basis with one element in every degree
of `Gamma`,

```text
dim_K K[x]/S=89,             dim_K S_(<=375)=376-89=287. (13)
```

It also proves that the six elements in `(5)` are a SAGBI basis: their
leading monomials generate the full initial algebra displayed in `(10)`.

Finally `Z` is monic of degree `18`, so `x` is integral over `S`.  By
Lüroth, if `Frac(S)=K(h)`, then `deg h` divides the degree of every element
of `S`; but `gcd(18,71)=1`.  Hence `deg h=1` and

```text
Frac(S)=K(x),                    normalization(S)=K[x]. (14)
```

## 4. Why the raw coupled packet has exactly 287 live restrictions

The finite-field sidecar rebuilds the raw target packet

```text
M_375={b^i c^j e^k:30i+21j+18k<=375},   |M_375|=1017 (15)
```

at every root in the completely split primes `137` and `163`.  In all eight
fibres its restriction rank is `287`, with the same pivot-column hash and
the same degree-gap list `(11)`.  Moreover the `395+174+395` coupled selector
uses precisely those `287` distinct pivot restrictions.

The positive minor at `(p,alpha)=(137,44)` has a valid characteristic-zero
consequence.  Its entries are the good reduction of the exact `K`-matrix;
its nonzero determinant is therefore the reduction of a nonzero element of
`K`.  Thus the corresponding `287` exact raw restrictions are independent.
They lie in `S_(<=375)`, whose exact dimension is `287` by `(13)`, so they
form a `K`-basis of the entire filtered restriction space.

This use of a modular minor is one-directional and sound: nonzero reduction
proves a characteristic-zero determinant is nonzero.  No modular rank
failure is used to infer characteristic-zero infeasibility.

Consequently the coupled lift should be viewed on the exact free module
`(8)`, not as a 1,017-monomial mystery.  The remaining `174` transverse
directions are a gap-module/cokernel problem on this canonical object.

## 5. Scope and reproduction

This theorem proves the restriction-algebra structure uniformly over the
field `K`, hence under all four embeddings `K -> C`.  It does not prove a
global conductor-ideal formula, `J_1=J_2=0`, a continuation through higher
stable orders, a Keller pair, or a counterexample to `JC(2)`.

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py

python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_restriction_basis_modular_probe.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_restriction_basis_modular_probe.py
```

Normal and optimized exact runs byte-match the stored transcript.  The
split-prime sidecar is explicitly `FINITE_FIELD_ONLY`; the exact module proof
is over `Q` and `K`.  **QED.**
