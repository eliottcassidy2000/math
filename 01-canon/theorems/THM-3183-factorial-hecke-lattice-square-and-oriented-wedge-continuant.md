---
id: THM-3183
title: "Factorial Hecke lattice square and oriented wedge continuant"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The scalar-moment and x-weighted factorial Gauss--Manin frames form an
  exact commuting square of integral lattices.  At the prime reset its
  cokernel lengths split as 2=1+1-0, and on exterior squares as 4=2+2-0.
  The scalar exterior transfer has an exact graph-distance-two continuant
  pivot with a complete factor wall.  A sharp same-Smith hostile proves
  that these data do not determine a projected Newton row or the empirical
  floor(s/2) Euclidean-depth staircase.
audit: >
  The exact symbolic companion checks the matrix intertwining, determinants,
  all four reset Smith profiles and their exterior-square profiles, the
  oriented two-step pivot, three same-Smith hostile tails, and the offset-six
  leading PRS pivots H and J.  Normal and optimized replay agree with the
  stored transcript and declared hashes.  An independent immutable audit
  rederived every matrix, Smith profile, length identity, pivot, wall, and
  hostile.  A separate pure-integer multinomial interpolation at p=0..20,
  with degree bounds 11 and 15, independently proves the H and J identities.
  MISTAKE-359 repairs an evidence-only float/truncated-window defect in the
  original top-jet helper: the maintained script now uses exact Rational
  normalization, constructs every required predecessor coefficient, and
  matches 66 direct integer A/B/R/S coefficients at three primes.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/factorial_hecke_lattice_wedge_continuant_thm3183.py
output: 05-knowledge/results/factorial_hecke_lattice_wedge_continuant_thm3183.out
script_sha256: 5da6dea97dd8cb51b332d4a133a448fc8d7f0f15a2c07463d721043cd976ac65
output_sha256: 2ba36ca1337b74b111300f1eac8fc45b422d499495e87d4a0080908d81630b19
hash_basis: LF-normalized bytes
---

# THM-3183 -- factorial Hecke lattice square and oriented wedge continuant

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3182 found two different integral Smith profiles for the same rational
quadratic-moment system.  The `x`-weighted frame has reset type `(1,p,p)`,
whereas the scalar companion has type `(1,1,p)`.  The discrepancy is not a
coordinate error.  It is one elementary modification of the output lattice.
The exact commuting square below records where that extra layer lives.

## 1. The two integral frames

Put

```text
q(x)=d-x+v x^2,                  Delta=1-4dv,
M_n=int_0^infinity q(x)^n e^(-x) dx,
X_n=int_0^infinity x q(x)^n e^(-x) dx,
D_n=d^n.                                                       (1)
```

For `n>=1`, use the two states

```text
z_n=(M_n,X_n,D_n)^t,            Y_n=(M_n,M_(n-1),D_n)^t.       (2)
```

THM-3182 gives `z_(n+1)=G_n z_n`, where

```text
G_n=
[-(n+1)                  2v(n+1)                  d       ]
[-(n+1)(2n+3+2d)        (n+1)(1+2v(2n+3))       d(2n+3)]
[0                       0                         d       ].   (3)
```

The scalar recurrence gives `Y_(n+1)=S_nY_n`, with

```text
S_n=[a_n b_n c_n; 1 0 0; 0 0 d],
a_n=2(n+1)(2n+1)v,
b_n=n(n+1)Delta,
c_n=d-n-1.                                                   (4)
```

Assume `2vDelta` is invertible and define

```text
P_n=
[1                                      0              0      ]
[(1+2(2n+1)v)/(2v)       nDelta/(2v)       -1/(2v)]
[0                                      0              1      ]. (5)
```

The integration-by-parts identity

```text
2vX_n=[1+2(2n+1)v]M_n+nDelta M_(n-1)-D_n                  (6)
```

says exactly that `z_n=P_nY_n`.  More strongly, direct multiplication gives
the rational Gauss--Manin intertwining

```text
P_(n+1) S_n = G_n P_n.                                     (7)
```

Thus `(7)` is a commuting square of integral lattices for `n>=1` whenever its displayed
denominators are units.  Its determinants are

```text
det P_n=nDelta/(2v),
det S_n=-dn(n+1)Delta,
det G_n=-d(n+1)^2Delta.                                    (8)
```

Equation `(7)` therefore conserves determinant length:

```text
val(det G_n)=val(det S_n)+val(det P_(n+1))-val(det P_n).    (9)
```

## 2. The prime-reset Hecke square

Let the base be the height-one DVR `Z_p[v]_(p)`, or an unramified DVR
extension in which `p` is a uniformizer.  Assume `p` is odd, put `d=p+s`,
and assume `svDelta` is a unit.  At `n=p-1`, the four maps in `(7)` have
Smith types

```text
P_(p-1): (1,1,1),             S_(p-1): (1,1,p),
P_p:     (1,1,p),             G_(p-1): (1,p,p).             (10)
```

For `P_(p-1)` this follows from the unit determinant in `(8)`.  The reduction
of `P_p` has rank two, has a unit `2x2` minor on the first and third
coordinates, and its determinant has valuation one.  The scalar reset
`S_(p-1)` likewise has rank two, contains the unit minor `-s`, and has
determinant valuation one.  THM-3182 proves the last profile.

If `ell(A)` denotes the length of the torsion cokernel of a full-rank integral
map, `(9)` becomes the exact Hecke balance

```text
2=1+1-0.                                                    (11)
```

For a rank-three map with invariant-factor valuations `(e_1,e_2,e_3)`, its
exterior square has valuations

```text
(e_1+e_2,e_1+e_3,e_2+e_3).                                (12)
```

Consequently the exterior-square profiles around the same square are

```text
Lambda^2 P_(p-1): (1,1,1),       Lambda^2 S_(p-1): (1,p,p),
Lambda^2 P_p:     (1,p,p),       Lambda^2 G_(p-1): (p,p,p^2), (13)
```

and their lengths obey

```text
4=2+2-0.                                                    (14)
```

The second height-one layer in the weighted reset has therefore not appeared
from a coordinate-free Smith invariant.  One unit of length belongs to the
scalar dynamics and one to the adjacent output-lattice modification.  The
commuting lattice square, rather than either isolated Smith tuple, is the
invariant object.

## 3. The oriented scalar wedge continuant

In the ordered scalar wedge basis

```text
w_01=M_n wedge M_(n-1),
w_02=M_n wedge D_n,
w_12=M_(n-1) wedge D_n,                                    (15)
```

the exterior transfer is

```text
Lambda^2 S_n=
[-b_n   -c_n       0   ]
[ 0      a_nd      b_nd]
[ 0       d         0  ].                                  (16)
```

The line spanned by `w_01` is invariant.  On the quotient by that line, the
transverse transfer is the Jacobi/continuant block

```text
K_n=[a_nd b_nd; d 0],             det K_n=-b_nd^2.          (17)
```

Orientation now adds information which Smith factors forget.  The source
`w_12` has no `w_01` component after one step, but after two steps

```text
(Lambda^2 S_(n+1))(Lambda^2 S_n)w_12
 =b_nd[-c_(n+1)w_01+a_(n+1)d w_02+d w_12].                 (18)
```

Its first visible coefficient is exactly

```text
-c_(n+1)b_nd=(n+2-d)n(n+1)Delta d.                         (19)
```

Hence the hidden-to-visible graph distance is exactly two whenever the five
factors in `(19)` are units.  The complete elementary wall is

```text
n(n+1)dDelta(n+2-d)=0.                                     (20)
```

This is the first rigorous continuant sidecar for the fixed-offset
Euclidean staircase.  It is oriented and index-sensitive; the Smith profile
alone cannot recover it.

## 4. Sharp same-Smith projection hostile

Let `pi` be a DVR uniformizer and take the fixed reset

```text
R=diag(1,pi,pi),             pi^(-1)Lambda^2R=diag(1,1,pi). (21)
```

All three unimodular tails

```text
U_infinity=I,             U_1=I+E_(0,2),
U_2=I+pi E_(0,2)                                             (22)
```

preserve the Smith type of `R`.  Nevertheless the `w_01` coefficient of

```text
(Lambda^2 U_i)(pi^(-1)Lambda^2R)w_12                       (23)
```

has respectively valuation

```text
infinity, 1, 2.                                             (24)
```

Thus even the complete reset Smith type and the named second-order source do
not determine when a fixed projected coordinate becomes visible.  The
oriented tail minors are indispensable.

## 5. Exact interface with the offset-six PRS

The distinction is already visible in THM-3176.  Use its offset-six notation
`A,B,R,S,T` and write `R_j=[v^j]R`, `S_j=[v^j]S`.  Exact top-jet calculation
gives

```text
R_(p+3)/(2p)!
 =-8 prod_(i=1)^5(p+i)(2p+1)(2p+3) H,                     (25)

S_(p+2)/(2p)!
 = 4 prod_(i=1)^5(p+i)(2p+7) J,                            (26)
```

where

```text
H=24p+109,
J=256p^4-27648p^3-365600p^2-1528800p-2096649.             (27)
```

Thus, in THM-3176's `p>=197` range, `H` and `J` are precisely the nontrivial
arithmetic wall factors of the leading coefficient-degree Schur-complement
pivots.  The next polynomial `K` of THM-3176 enters the inherited connection coefficient

```text
P_0=4(p+6)(2p+1)K                                          (28)
```

of the third pseudoquotient.

By contrast, substituting `d=p+6` and `n=p+j` into the bare two-step pivot
`(19)` gives only

```text
(j-4)(p+j)(p+j+1)(p+6)(1-4(p+6)v).                        (29)
```

The arithmetic walls `H,J,K` therefore live in the coefficient-degree
Schur complements and their connection data, not in the reset Smith tuple or
the bare time-index pivot.  Proving the empirical law

```text
first separating Euclidean depth=floor(s/2)                (30)
```

still requires a closed identification between successive PRS pivots and an
oriented continuant product, together with its arithmetic wall atlas.

## 6. Scope and exact evidence

The theorem proves an integral lattice square, an exterior-length identity,
one exact two-step continuant pivot, and a projection no-go.  It does not
prove `(30)`, arbitrary fixed-offset closure, `NC(2)`, `GMC(2)`, or an LRC
consequence.

Run

```text
python 04-computation/factorial_hecke_lattice_wedge_continuant_thm3183.py
python -O 04-computation/factorial_hecke_lattice_wedge_continuant_thm3183.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
exact symbolic and integer arithmetic only.  It checks `(7)--(14)`, the
complete exterior matrix `(16)`, the two-step formula `(18)--(20)`, all three
hostile tails `(21)--(24)`, the offset-six leading identities `(25)--(27)`,
and the bare pivot `(29)`.  Equation `(28)` is inherited verbatim from
THM-3176 and is not an independent companion check here.

**QED.**
