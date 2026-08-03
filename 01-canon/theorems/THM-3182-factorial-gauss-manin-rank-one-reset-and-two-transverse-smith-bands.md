---
id: THM-3182
title: "Factorial Gauss--Manin rank-one reset and two transverse Smith bands"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The quadratic factorial moments and their first x-weighted companions form
  an exact three-state Gauss--Manin transfer.  At d=p+s each p-period has a
  rank-one Frobenius reset, the entire state descends from
  degree p+a to degree a, and its generic height-one Smith type in the
  x-weighted integral lattice is (1,p,p).  This exposes two transverse
  p-adic directions in that specified lattice but does not by itself prove
  the observed floor(s/2) Euclidean-depth staircase.
audit: >
  Two independent immutable audits rederived the transfer and scalar
  recurrence, rank-one reset and all-a descent, both Smith types, complete
  exterior-square first and second layers, scalar-companion gauge/index
  boundary, discriminant wall, and p=2/odd-p scopes.  The repaired companion
  checks the actual scalar unit minor -s.  Fresh normal and optimized replay
  agree with the 13-line stored transcript and declared LF hashes.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on: []
related:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
script: 04-computation/factorial_gauss_manin_rank_one_reset_thm3182.py
output: 05-knowledge/results/factorial_gauss_manin_rank_one_reset_thm3182.out
script_sha256: 0eb5eca0ae0664ed56005c13638f2c00f4376eb0798c05eba3f35173b630491e
output_sha256: b63165b6a95dbc124e041f871f3ee0ddf3d9f2e7195c8e4d6df029574658e546
hash_basis: LF-normalized bytes
---

# THM-3182 -- factorial Gauss--Manin rank-one reset and two transverse Smith bands

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The factorial Frobenius projector has a stronger state-space form.  Its
first singular transfer forgets two directions modulo `p`, but it remembers
both directions to first `p`-adic order in the `x`-weighted integral lattice.
This gives a precise candidate source for the two positive-height bands in
the fixed-offset Newton calculations, with a load-bearing lattice-gauge
boundary recorded in Section 5.

## 1. Three-state integration-by-parts transfer

Let

```text
q(x)=d-x+v x^2,
M_n^[d](v)=int_0^infinity q(x)^n e^(-x) dx,
X_n^[d](v)=int_0^infinity x q(x)^n e^(-x) dx,
D_n^[d]=d^n.                                                (1)
```

These integrals mean the exact factorial functional, so all entries lie in
`Z[d,v]`.  Put

```text
z_n^[d]=(M_n^[d],X_n^[d],D_n^[d])^t.                        (2)
```

For every `n>=0`,

```text
z_(n+1)^[d]=T_n^[d](v) z_n^[d],                             (3)

T_n^[d]=
[ -(n+1)                    2v(n+1)                  d       ]
[ -(n+1)(2n+3+2d)          (n+1)(1+2v(2n+3))       d(2n+3)]
[  0                        0                         d       ]. (4)
```

Its determinant is

```text
det T_n^[d]=-d(n+1)^2(1-4dv).                               (5)
```

Thus the quadratic discriminant `1-4dv` is exactly the singular wall of the
Gauss--Manin transport.

### Proof

Since `q'=-1+2vx`, integration of `(q^(n+1)e^(-x))'` gives

```text
M_(n+1)=d^(n+1)-(n+1)M_n+2v(n+1)X_n.                       (6)
```

The boundary term at zero is `-d^(n+1)`.  Next integrate
`(x q^(n+1)e^(-x))'`.  Both endpoint terms vanish, and

```text
xq'=2q-2d+x.                                                (7)
```

Consequently

```text
X_(n+1)
 =(2n+3)M_(n+1)-2d(n+1)M_n+(n+1)X_n.                       (8)
```

Substituting `(6)` into `(8)` gives the second row of `(4)`; the third row is
tautological.  Expanding the determinant of the upper `2x2` block gives
`(n+1)^2(4dv-1)`, proving `(5)`.

Eliminating `X_n` from two consecutive copies of `(6)` recovers the scalar
recurrence

```text
M_(n+1)
 =d^n(d-n-1)+2(n+1)(2n+1)vM_n
   +n(n+1)(1-4dv)M_(n-1),                                  (9)
```

so `(3)` is a genuine state lift of the quadratic factorial recurrence, not
an unrelated auxiliary construction.

## 2. Rank-one prime reset and full state descent

Let `p` be prime, let `s>=1` be an integer with `p` not dividing `s`, and put

```text
d=p+s.                                                       (10)
```

At the first singular index `n=p-1`, reduction of `(4)` gives

```text
T_(p-1)^[p+s]
 ==[0 0 s; 0 0 s; 0 0 s]                         (mod p).  (11)
```

Fermat gives `D_(p-1)^[p+s]==1`, hence `(11)` resets every incoming state to

```text
z_p^[p+s]==s(1,1,1)^t=s z_0^[s]                   (mod p).  (12)
```

For every `j>=0`, coefficientwise in `v`,

```text
T_(p+j)^[p+s](v)==T_j^[s](v)                      (mod p).  (13)
```

Induction from `(12)` therefore proves the full state-level Frobenius
descent

```text
z_(p+a)^[p+s](v)==s z_a^[s](v)                    (mod p)   (14)
```

for every `a>=0`.  The first coordinate of `(14)` is THM-3148's fixed-offset
residual-polynomial congruence.  The second coordinate is new: the first
`x`-weighted response descends through the same reset.

The condition `p` not dividing `s` is load-bearing.  If `p|s`, then `d` is
not a unit, the matrix in `(11)` is zero rather than rank one, and Fermat does
not give `(12)`.

## 3. The two transverse Smith directions

Work first over the height-one DVR

```text
R=Z_p[v]_(p),                                                (15)
```

whose residue field is `F_p(v)`.  Equivalently, specialize `v` in any
unramified DVR extension subject to

```text
p does not divide d(1-4dv).                                 (16)
```

The Smith type of the reset matrix is

```text
Smith_R(T_(p-1)^[d])=(1,p,p)                                (17)
```

up to unit factors.

Indeed, `d` is a unit entry, so the first determinantal divisor is a unit.
The reduction `(11)` has rank one, so every `2x2` minor is divisible by `p`.
But the minor on rows `(1,3)` and columns `(1,3)` is exactly

```text
-pd,                                                         (18)
```

so the second determinantal divisor has valuation one.  Finally `(5)` gives

```text
v_p(det T_(p-1))=2                                           (19)
```

under `(16)`.  The invariant-factor valuations are therefore `(0,1,1)`,
which is `(17)`.

On the discriminant wall `p|(1-4dv)`, `(11)` still has rank one but the
determinant has valuation at least three.  For a simple wall the Smith type
becomes `(1,p,p^2)`.  Thus the quadratic discriminant is not cosmetic: it is
exactly where one transverse band thickens.

## 4. Exterior-square first layer

The filtered Pluecker form is even more explicit.  In the ordered basis

```text
(M wedge X, M wedge D, X wedge D)                           (20)
```

of the exterior square, put `C=Lambda^2 T_(p-1)`.  Direct `2x2` minors give

```text
C=
[p^2(4pv+4sv-1)  2p(p+s)^2             -p(p+s)             ]
[0                 -p(p+s)                2pv(p+s)            ]
[0                 -p(p+s)(4p+2s+1)      p(p+s)(4pv+2v+1)].  (21)
```

After dividing by `p` and reducing,

```text
p^(-1)C ==
[0  2s^2       -s       ]
[0  -s          2sv     ]
[0  -s(2s+1)    s(2v+1)]                         (mod p).    (22)
```

Under `(16)`, `(22)` has rank two.  Its right kernel is generated by
`M wedge X`, its left kernel is generated by `(1,-1,1)`, and each of its
three nonzero `2x2` minors is

```text
-s^2(1-4sv).                                                (23)
```

The missing input direction is not lost; it appears one layer later:

```text
p^(-2) C(M wedge X)
 ==-(1-4sv)(M wedge X)                           (mod p).   (24)
```

Thus `(17)` induces exterior-square Smith type

```text
(p,p,p^2).                                                  (25)
```

Equations `(22)--(24)` identify the precise first-order Pluecker plane and
the single second-order direction.  They are stronger than the scalar
determinant valuation: the left covector `(1,-1,1)` records the unique
first-layer relation among the three projected wedge coordinates.  Every
Smith type, kernel, and transverse direction in Sections 3--4 is a property
of the specified integral lattice `z=(M,X,D)`, not an invariant under
rational conjugacy of the underlying Gauss--Manin system.

## 5. Filtered-holotopy meaning and boundary

The reset separates the transport into

```text
one mod-p survivor  +  two height-one transverse directions. (26)
```

The fixed tail `T_p,...,T_(p+s-2)` consists of matrices affine in `v`.
Consequently every fixed-offset Newton calculation is the projection of a
two-column transverse lattice through a bounded affine tail.  When that
projection is compatible with the `x`-weighted lattice, `(17)` is a candidate
mechanism for the two positive valuation bands and for pivot mutation of the
active Euclidean row.  The Smith data alone do not supply that compatibility.

The words *transverse lattice* are load-bearing.  For odd `p`, assume `v`
and `Delta=1-4dv` are units and instead use the scalar companion state

```text
Y_n=(M_n,M_(n-1),D_n)^t.                                   (27)
```

Equation `(9)` gives its transfer

```text
Y_(n+1)=
[a_n b_n c_n; 1 0 0; 0 0 d]Y_n,

a_n=2(n+1)(2n+1)v,
b_n=n(n+1)Delta,
c_n=d-n-1.                                                  (28)
```

At `n=p-1` this companion matrix has Smith type `(1,1,p)`, not
`(1,p,p)`: its reduction has rank two and its determinant has valuation one.
The exact change to the `x`-weighted state is

```text
2vX_n=[1+2(2n+1)v]M_n+n Delta M_(n-1)-D_n.                 (29)
```

Its determinant is `n Delta/(2v)`.  It is a unit at `n=p-1` but gains one
factor of `p` at the output index `n=p`.  Thus the extra transverse factor in
`(17)` is a genuine layer of the chosen weighted-response lattice, not a
coordinate-free Smith invariant of the rational Gauss--Manin system.  An
index-dependent non-unimodular gauge can move that layer between the reset
and the output lattice; in the scalar companion framing one band is absorbed
by this output index.

What `(17)` does **not** prove is equally important.  Projection through the
tail can cancel coordinates, and successive pseudo-quotients depend on
continuant minors of the whole tail.  The observed statement

```text
first separating Euclidean depth=floor(s/2)                 (30)
```

away from arithmetic pivot walls remains open.  A proof of `(30)` requires a
closed transfer-minor/continuant formula plus a wall atlas.  No arbitrary
fixed-offset theorem, `NC(2)`, `GMC(2)`, or `LRC(14)` consequence is claimed.

## 6. Exact evidence

Run

```text
python 04-computation/factorial_gauss_manin_rank_one_reset_thm3182.py
python -O 04-computation/factorial_gauss_manin_rank_one_reset_thm3182.py
```

and compare LF-normalized bytes with the declared output.  The companion is
pure integer arithmetic.  It reconstructs `M_n,X_n` directly from their
multinomial formulas, checks `(3)` and `(9)`, verifies `(5)` on a finite exact
bank, tests the complete polynomial congruence `(14)`, and checks both the
generic `(1,p,p)` determinantal divisors and the simple discriminant-wall
thickening.  It reconstructs the complete exterior-square matrix, its
rank-two first layer, both kernels, all three nonzero minors, and the
returning second-order direction.  It also freezes a `p|s` hostile.  There
is an independent exact check of the scalar-companion gauge identity, its
`(1,1,p)` reset, and the one-factor output index.  There is no floating point,
random sampling, imported executable, or assertion-sensitive test.

**QED.**
