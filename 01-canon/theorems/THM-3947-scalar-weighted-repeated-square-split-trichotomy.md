---
id: THM-3947
title: "Scalar-weighted repeated-square splits have a three-parabola trichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the complete
  scalar weighting of the THM-3944 internal split p1-p0=alpha G^2, the common
  discriminant is controlled by one cubic F(z)=(Az+B)^2-4z^3.  Away from two
  scalar seams it is the reduced union of three distinct smooth one-place
  parabolas.  At r^2=alpha it is a doubled p0-line plus one parabola; at
  r^2=-omega alpha it is a doubled p1-parabola plus one other parabola.  A
  triple component is impossible for nonzero parameters.  The two collision
  seams are index-swapped copies of THM-3944's conductor geometry.  This is a
  factorization classification, not a computation of the generic normal
  quadratic order or its Cardano class lattice.
source: root / arbitrary-scalar completion of THM-3944, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (synthesize_thm3933, 2026-08-24).  The audit
  rederived the two scalar moduli, cubic discriminant, both seam identities,
  the three distinct nonzero generic roots, both collision factorizations,
  absence of a hidden G-factor or triple component, the smooth one-place
  conic ledger, and the index-swapped THM-3944 geometry.  The assertion-free
  22-gate companion byte-matches in normal and optimized modes after LF
  normalization; frozen output and all hashes agree.
depends_on: []
related:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3946-affine-internal-factor-split-two-end-conductor-collision-dichotomy
script: 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
output: 05-knowledge/results/jc2_scalar_weighted_repeated_square_split_thm3947.out
script_sha256: d6f96faa5c62acc961a43f7154eeac1bcabf065986acf67e2c3ddd04b4212449
output_sha256: fc19766822e0dda8f28c05793aa2ddc38fbdaa00624e44ca57bd6088c6a2f8cf
semantic_sha256: 7b1ec493d249a8cf6bf757089e97bc96288431378e0fd7192b2eeed73bc8e106
hash_basis: raw LF bytes
---

# THM-3947 -- every scalar weight produces three parabolas or a conductor seam

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Fix a primitive cube
root `omega`, put `delta=omega-omega^2`, and let `P,G` be affine coordinates.
For nonzero `alpha,r in k`, set

```text
p_0=P,                         p_1=P+alpha G^2,
L_i=p_1-omega^i p_0,
D=(q_1-q_0)/2=r G L_1,
S=(q_1+q_0)/2=(alpha/r)G L_2.                         (1)
```

Then `D S=p_1^3-p_0^3`, so `(1)` determines polynomials `q_0=S-D` and
`q_1=S+D` with a common double-torus discriminant

```text
H=q_0^2-4p_0^3=q_1^2-4p_1^3.                          (2)
```

Exactly one of the following alternatives holds.

```text
(G) r^2 != alpha and r^2 != -omega alpha:
    H is the reduced union of three distinct smooth parabolas;

(C0) r^2=alpha:
    H=-P^2(4P+3alpha G^2);

(C1) r^2=-omega alpha:
    H=-4(P+alpha G^2)^2(P+alpha G^2/4).                (3)
```

Every reduced component in `(3)` is isomorphic to `A1` and has one smooth
place at infinity.  In `(G)` the **full** branch is nevertheless reducible.
In `(C0)` and `(C1)` it is nonreduced; the doubled component is respectively
`p_0=0` and `p_1=0`.  No parameter gives a triple component.

This exhausts only the scalar-weighted allocation `(1)`.  Unequal coprime
factors `p_1-p_0=FG`, simultaneous internal splits in more than one `L_i`,
the generic quadratic normalization and Cardano lattice, source attachment,
and JC(2) remain **OPEN**.

## 1. The internal split and its two scalar moduli

Since `L_0=alpha G^2`,

```text
L_0L_1L_2=p_1^3-p_0^3=D S.                            (4)
```

Writing `q_0=S-D` gives

```text
q_0=G(AP+BG^2),
A=(alpha/r)(1-omega^2)-r(1-omega),
B=alpha(alpha/r-r).                                    (5)
```

Therefore

```text
H=G^2(AP+BG^2)^2-4P^3.                                (6)
```

On the chart suggested by the weight assignment `wt(G)=1, wt(P)=2`, define

```text
F(z)=(Az+B)^2-4z^3.                                    (7)
```

Then exactly

```text
H=G^6 F(P/G^2).                                        (8)
```

The ordinary cubic discriminant is

```text
disc_z(F)=-16 B^3(A^3+27B).                            (9)
```

The two factors in `(9)` retain the scalar allocation:

```text
rB=alpha(alpha-r^2),
A^3+27B=3delta(r^2+omega alpha)^3/r^3.                (10)
```

Since `alpha,r` are nonzero and `1+omega!=0`, the two zero loci in `(10)` are
disjoint.  This already proves that `(3)` is exhaustive.

## 2. The generic row is three distinct one-place parabolas

In case `(G)`, equations `(9)--(10)` show that `F` has three distinct roots.
Its constant term is `B^2!=0`, so none is zero.  Write them as
`z_1,z_2,z_3`.  Because the leading coefficient of `F` is `-4`, equation
`(8)` becomes

```text
H=-4 product_(j=1)^3(P-z_j G^2).                       (11)
```

The three factors are pairwise nonassociate irreducibles in `k[P,G]`; hence
the full divisor is reduced and reducible.  Each component is the graph

```text
P=z_j G^2,                                              (12)
```

so its affine normalization is `A1_G`.  Its projective closure

```text
P Z-z_j G^2=0                                          (13)
```

is a smooth conic.  At `Z=0`, equation `(13)` gives the unique point
`[P:G:Z]=[1:0:0]`, and `partial(13)/partial Z=P` is nonzero there.  Thus each
component has exactly one smooth normalization place at infinity.

This is a genuine improvement in componentwise place geometry, but not an
irreducible full discriminant and not yet a normal-order or Keller result.

## 3. The two collision seams

If `r^2=alpha`, then `B=0` and `A^2=-3alpha`.  Substitution in `(6)` gives

```text
H=-P^2(4P+3alpha G^2),                                 (14)
```

the scalar form of THM-3944: a doubled `p_0` line and one smooth parabola.

Suppose instead that `r^2=-omega alpha`.  Equivalently
`alpha=-omega^2r^2`.  Direct substitution gives

```text
H=-4(P+alpha G^2)^2(P+alpha G^2/4).                    (15)
```

Thus the doubled component is `p_1=0`, and the other component is a distinct
parabola.  After exchanging the two torus presentations and using
`P'=p_1`, equation `(15)` reads

```text
H=-P'^2(4P'-3alpha G^2),                               (16)
```

which is the same divisor and conductor geometry as `(14)` after a nonzero
scalar rescaling.

For clarity, on the second seam `B=-A^3/27`, and

```text
F(z)=-4(z-A^2/9)^2(z-A^2/36).                          (17)
```

The two roots in `(17)` are distinct.  On the first seam, `B=0` but
`A^2=-3alpha!=0`, so `F=z^2(A^2-4z)` also has a distinct simple root.
Consequently neither seam can produce a triple component.

THM-3944 supplies the detailed normalization, conductor, and character audit
for the first geometry.  Its index-swapped version applies to the second.
Those conclusions are contextual here; the present theorem proves the full
scalar factorization trichotomy itself.

## 4. Exact reproducibility

Run

```bash
python3 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
python3 -O 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
```

The assertion-free 22-gate companion verifies `(4)--(10)`, both explicit
collision factorizations `(14)--(17)`, the absence of a triple component, and
the smooth one-place projective ledger `(12)--(13)`.  Both modes reproduce
the frozen LF transcript.
