---
id: THM-2406
title: "Degree-eighteen H4 weighted-pole deep-wall collapse"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the a_0=1 weighted chart of THM-2389, four necessary pole
  equations already force both 126D=25B^2 and 20BC+21W=0. The
  characteristic-zero certificate is over the complete slope splitting
  field Q(epsilon), epsilon^6+10epsilon^4+25epsilon^2+23=0, and is
  saturated only by the proved pole units a_3 A_3(1). After losslessly
  removing those unit factors, the four equations have degrees
  2,3,4,6 and a four-element triangular Groebner basis with a quadratic
  quotient. Both wall functions, and the four omitted synchronization
  equations, have exact zero normal form. Thus a genuine H_4 survivor
  would lie on THM-2345's already-closed complete common-root wall.
  Together with the preceding degree-eighteen closures, this closes the
  inherited genuine nonsplit, polynomial exact-square-prefix reduced
  degree-eighteen branch. It does not close terminal branches outside
  THM-2262's reduction or prove JC(2) or DC(2).
source: codex-2026-07-26-h4-weighted-pole-wall
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2389-degree-eighteen-h4-three-pole-hermite-jet-synchronization
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
  - THM-2376-degree-eighteen-h2-coprime-cube-locus-closure
  - THM-2386-degree-eighteen-h4-common-root-elimination
  - THM-2387-degree-eighteen-h4-elliptic-three-isogeny-atlas
script: 04-computation/jc2_degree18_h4_weighted_pole_wall_thm2406.py
output: 05-knowledge/results/jc2_degree18_h4_weighted_pole_wall_thm2406.out
script_sha256: 4f527ba8b810b7b8039008207dafded9a837f8584ea1c7104a8176690a0d0126
output_sha256: ec7dbb0955f302a1388aec8a9c5b162173b70c623f0ad72f6fb186fd4c94e193
hash_basis: working-tree bytes (LF)
---

# THM-2406 -- the H4 pole section collapses to the deep wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2389 replaces a putative degree-eighteen `H_4` identity by a
three-pole Hermite system. This theorem performs the remaining
characteristic-zero elimination. Its conclusion is

```text
genuine H_4 survivor
  => 126D-25B^2=0
  => contradiction by the complete-wall conclusion of THM-2345.      (1)
```

The computation also forces the deeper equation

```text
20BC+21W=0.                                                        (2)
```

No normalization by `B,C,D`, or `W` is made. In particular, the
coordinate axes and vanishing-parameter boundaries in (1)--(2) are not
discarded.

## 1. A small model of the complete slope splitting field

Use the integer-scaled spectral coordinate

```text
w=(1701/2)v.
```

At infinity its three slopes are the roots of

```text
g(alpha)=alpha^3+2940alpha+30184.                                  (3)
```

The cubic is irreducible over `Q`, and

```text
Disc(g)=-126247730112=-23*74088^2.                                 (4)
```

Thus its Galois group is `S_3` and its splitting field has degree six.
For an ordered root triple `(alpha_0,alpha_1,alpha_inf)`, put

```text
epsilon=(alpha_1-alpha_inf)/42.
```

Direct elimination gives the unexpectedly small polynomial

```text
h(epsilon)
 =epsilon^6+10epsilon^4+25epsilon^2+23.                             (5)
```

It is irreducible and separable. In

```text
K=Q[epsilon]/(h)
```

the three roots are

```text
alpha_0   =154/(3epsilon^2+5),

alpha_1   =-alpha_0/2+21epsilon,

alpha_inf =-alpha_0/2-21epsilon.                                   (6)
```

The denominator `3epsilon^2+5` is nonzero and has an explicitly checked
inverse in `K`. The companion checks (3)--(6) and this unit identity
exactly. Since `K` has degree six and
contains all three roots, it is the full splitting field. Every one of
the six slope orderings is therefore a Galois conjugate of (6).
Checking the single ordered system over `K` is faithful to all six
orderings; it is not a check of one numerical branch.

## 2. The root-free weighted chart

Write the three-pole map as in THM-2389:

```text
y=A_3(t)/[t(t-1)],

A_3=a_0+a_1t+a_2t^2+a_3t^3,

a_0 a_3 A_3(1)!=0.                                                  (7)
```

The weighted action rescales all `a_i` together. Since `a_0!=0`, set

```text
a_0=1,          x=a_1,          y_2=a_2,          z=a_3.             (8)
```

Here `y_2` is a coefficient, not the target coordinate `y`.
The two surviving pole units are

```text
z!=0,                    A_3(1)=1+x+y_2+z!=0.                       (9)
```

Use the Hermite polynomial from THM-2389:

```text
R_6=R_*+kappa[t(t-1)]^2.                                           (10)
```

In the integer-scaled spectral coordinate the homogenized identity is

```text
N
 =R_6^3+12P^h(A_3,t(t-1))R_6+56Q^h(A_3,t(t-1))=0,                  (11)
```

where

```text
P^h
 =245A_3^4+1890BA_3^2d^2+(-24300B^2+122472D)d^4,

Q^h
 =539A_3^6+11340BA_3^4d^2+183708CA_3^3d^3
  +(72900B^2-367416D)A_3^2d^4
  +(2361960BC+2480058W)A_3d^5,

d=t(t-1).                                                           (12)
```

The slope and missing-linear-jet conditions cancel degrees eighteen
and seventeen, so `deg_t N=16`.

## 3. Four necessary equations

The order-two equations at `t=0,1` reconstruct `kappa,B`. Their
determinant is

```text
nonzero element of K times A_3(1)^4.                                (13)
```

Thus (9) makes the reconstruction unique. The order-three,
order-four, and order-five equations at `t=0` then reconstruct
`C,D,W` through the nonzero pivots already proved in THM-2389.
The companion independently verifies that the determinant constant and
all three reconstruction pivots are nonzero elements of `K` and
multiplies each by its exact inverse.
After exact cancellation, all five reconstructed parameters are
polynomials in `(x,y_2,z)` over `K`; no further denominator remains.

Every full synchronized survivor necessarily satisfies the following
four equations:

```text
E_B      : order-two matching at infinity,

E_C      : order-three matching at t=1,

E_D      : order-four matching at t=1,

E_6      : the sixth-order lock at t=0.                              (14)
```

The raw cross-multiplied numerators contain powers of the pole units.
Exact division gives

```text
E_B =z^4       e_B,          deg(e_B)=2,  terms(e_B)=8,

E_C =A_3(1)^3 e_C,           deg(e_C)=3,  terms(e_C)=14,

E_D =A_3(1)^2 e_D,           deg(e_D)=4,  terms(e_D)=29,

E_6 =             e_6,       deg(e_6)=6,  terms(e_6)=55.            (15)
```

Dividing by the first three factors is lossless precisely because of
(9). This unit removal is load-bearing computationally: retaining the
factors creates a very large intermediate basis, while (15) exposes the
small actual ideal.

## 4. Exact localized elimination

Introduce one saturation variable `eta` and form

```text
I=(
   eta*z*A_3(1)-1,
   e_6,e_B,e_C,e_D
  )
  in K[eta,z,y_2,x].                                                 (16)
```

With grevlex variable order `(eta,z,y_2,x)`, the exact reduced basis has
four elements. Define

```text
R=-253+793epsilon^3+775epsilon+100epsilon^5.
```

Three of the basis elements have the compact form

```text
q(x)
 =x^2
  +[3(253+83epsilon^3+173epsilon+6epsilon^5)/253]x
  +3(253+249epsilon^3+519epsilon+18epsilon^5)/506,

z +(R/506)x+R/253,

y_2-(R/253)x-3R/506.                                                (17)
```

The fourth is monic linear in `eta`; its `eta` coefficient and inverse
are both exactly `1` in `K[x]/(q)`, and its remaining coefficient is
printed in the stored transcript. Thus the localized quotient is generated by
`1,x` and is at most quadratic over `K`. The basis is not the unit
ideal: it retains the boundary solutions that the wall theorem will
exclude.

Now reconstruct `(B,C,D,W)` and reduce the two invariant wall
functions by this exact basis. The normal forms are

```text
NF_I(126D-25B^2)=0,

NF_I(20BC+21W)=0.                                                   (18)
```

These are ideal-membership statements over the characteristic-zero
splitting field. They do not rely on point enumeration, radicality, or
finite-field lifting.

The calculation is stronger than needed. The four synchronization
equations omitted from (14) also have zero normal form:

```text
E_C,inf=E_D,inf=E_W,1=E_W,inf=0 in K[eta,z,y_2,x]/I.                (19)
```

Hence (14) is a lossless four-equation replacement for the full
THM-2389 synchronization system on the pole chart.

## 5. H4 is empty

Let a genuine `H_4` survivor exist. THM-2389 supplies (7), one of the
six slope orderings, the Hermite reconstruction, every equation in
(14), and the pole units (9). Normalize by (8) and apply the appropriate
Galois conjugate of (16)--(18). Then

```text
126D=25B^2.                                                         (20)
```

On (20), the structured covariants satisfy

```text
P(0)=0,                         Q(0)=0.                              (21)
```

Indeed the relevant coefficient identities are

```text
122472*25=24300*126,

367416*25=72900*126.                                                (22)
```

Equation (2) additionally kills the linear coefficient of `Q` because

```text
2480058*20=2361960*21.                                              (23)
```

The primary handoff is THM-2345: its complete-wall conclusion says that
(20), with `B=0` and every axis included, carries no degree-eighteen
Keller trajectory. This is the required contradiction. Alternatively,
(21) contradicts THM-2386's coprimality conclusion for a genuine
`H_4` survivor. Therefore the `H_4S_4^2` stratum is empty.

## 6. The inherited degree-eighteen corollary

Within THM-2262/2297's inherited genuine nonsplit, polynomial
exact-square-prefix reduction, the degree-eighteen branch ledger is:

| Degree-eighteen branch | Exact closure |
|---|---|
| support at most two in `(B,C,D,W)` | THM-2314, THM-2316, THM-2320, THM-2324, THM-2328 |
| higher-support `H_0S_6^2` | THM-2335 |
| higher-support `H_2S_5^2` | THM-2371 and THM-2376 |
| higher-support `H_4S_4^2` | this theorem, via THM-2389 and THM-2345 |

THM-2332 proves that these exhaust the remaining genus-zero square
classes in that reduction. Consequently the **inherited reduced
degree-eighteen Keller trajectory branch is empty**.

This is a degree-specific, reduction-specific closure. Terminal branches
outside THM-2262's hypotheses, split/even order-raising, the formal Weyl
lift, `JC(2)`, and `DC(2)` remain open.

## 7. Exact controls and the sharp lock boundary

Run

```bash
python3 04-computation/jc2_degree18_h4_weighted_pole_wall_thm2406.py
python3 -O 04-computation/jc2_degree18_h4_weighted_pole_wall_thm2406.py
```

Both modes are byte-identical to the stored output. The companion:

1. proves irreducibility and separability of (3) and (5), checks the
   discriminant and all root formulas;
2. rebuilds (10)--(14) from the structured covariants;
3. checks every exact unit division, the slope denominator, the
   determinant constant, all reconstruction pivots and their inverses,
   and the monic `eta` coefficient;
4. computes (16)--(19) over `K`;
5. specializes the exact basis at two independent split places:

   ```text
   p=59,  epsilon=1,  slopes=(34,4,21),
   q=x^2+18x+24, roots 14,27;

   p=101, epsilon=-5, slopes=(60,67,75),
   q=x^2+45x+14, no F_101 root;
   ```

6. searches the complete `F_59` pole chart before imposing `E_6` and
   finds the explicit hostile

   ```text
   (x,y_2,z)=(58,34,47),

   E_B=E_C=E_D=0,

   126D-25B^2=22,          E_6=2.                                  (24)
   ```

Thus the sixth-order lock is genuinely load-bearing: the first three
matching equations alone do not force the wall. The finite-place
checks are controls for the characteristic-zero computation, not its
logical basis.

## 8. Independent audit

The independent hostile audit reconstructed the full degree-six slope
field, checked the `a_0=1` weighted normalization and the exact
`a_3 A_3(1)` localization, and verified every reconstruction pivot and
denominator. It converted the four-element grevlex basis to lex order by
FGLM, reduced all six Buchberger `S`-pairs to zero, and independently
reduced all nineteen coefficients of the reconstructed numerator and
the four omitted synchronization equations. Both execution modes and
declared hashes match. Finally, it checked both handoffs: the wall
contradicts THM-2386 directly and lies in THM-2345's complete closed
wall. No mathematical defect was found.

The transcript labels explicitly report Boolean vanishing checks as
`*_vanishes=1`; they are not polynomial values.  Normal and optimized
execution reproduce the stored transcript byte-for-byte, and the LF
hashes in the frontmatter were independently verified.
