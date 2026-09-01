---
id: THM-4298
title: "Weighted-face source-normal unimodular visibility transform"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Under the reduced (2,3)
  source-normal substitution p=t(1+x^2t), y=xtp, every labelled residual
  weight-M face for wt(p)=2, wt(y)=3 maps losslessly to a minimal consecutive
  leading coefficient flag on 2n-d=M. The transform is unit lower triangular
  over Z and has an explicit alternating-binomial inverse; the minimal flag
  occupies rows ceil(M/2) through floor(2M/3), while later diagonal rows are
  redundant consequences. For M=12 it sends (U,W,Z) to
  (h0,h1,h2)=(U,6U+W,15U+5W+Z) and gives exact source-normal equations for all
  four coefficient walls. This is a visibility/inversion theorem, not a
  Darboux lift, seam-entry theorem, or JC(2) result.
source: root/jc2-planar-20260831
related:
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4005-reduced-two-three-live-seam-invariant-support-atlas
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
primary_script: 04-computation/jc2_weighted_face_source_normal_visibility_thm4298.py
primary_output: 05-knowledge/results/jc2_weighted_face_source_normal_visibility_thm4298.out
primary_script_sha256: 13398ceaadfff89b43233c8c136c0dc233b0e82113ed4b599f268c8ed08b58d7
primary_output_sha256: e0922e0419e6df3c1c459f9c7d90e1297005726076538402844c329e64f1bd17
independent_audit_script: 04-computation/jc2_weighted_face_source_normal_visibility_independent_audit_thm4298.py
independent_audit_output: 05-knowledge/results/jc2_weighted_face_source_normal_visibility_independent_audit_thm4298.out
independent_audit_script_sha256: 3bb989a46980b8e3af4c297971bfe9688958b017644b5adb0fdc7e0b9e8c3991
independent_audit_output_sha256: 347cb8e514b29c9227dc81a1aacff2fe55ac24991c908bedc7a50224ab5fe5e3
hash_basis: raw LF bytes
audit: >
  PASS. The primary standard-library path checks 68,921 substituted monomial
  terms, every nonempty face through M=200, 77,486 transform entries and
  closed inverses, a mixed polynomial through M=100, and 4,913 M=12 wall
  triples. An independent SymPy path starts from the generating-function
  substitution, checks both matrix products, forward/reverse identities, and
  the symbolic wall dictionary. LF-normalized normal and optimized streams
  byte-match the frozen outputs.
---

# THM-4298 -- Weighted-face source-normal unimodular visibility transform

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS IS A LOSSLESS
COEFFICIENT TRANSFORM AND ROW-DEPTH INVOICE. IT DOES NOT PROVE SEAM ENTRY OR
`JC(2)`.**

## 1. Statement

Let `K` be any commutative ring and give `K[p,y]` the weights

```text
wt(p)=2,                         wt(y)=3.                     (1)
```

Use the reduced `(2,3)` source-normal substitution

```text
p=t(1+x^2t),                    y=xtp.                        (2)
```

For `M>=0`, list the monomials of weight `M` in increasing `y`-exponent:

```text
j_r=j_0+2r,                    i_r=(M-3j_r)/2,               (3)
```

where `j_0` is the parity of `M`. The list is empty only for `M=1`. Write

```text
R_M=sum_r c_r p^(i_r)y^(j_r),
h_s=[x^(j_s)t^((M+j_s)/2)] R_M(p(x,t),y(x,t)).             (4)
```

> **Theorem.** The coefficient vectors are related by
>
> ```text
> h_s=sum_(r<=s) binom((M-j_r)/2,s-r)c_r,                  (5)
> c_s=sum_(r<=s)(-1)^(s-r)binom((M-j_r)/2,s-r)h_r.         (6)
> ```
>
> Thus `(5)` is unit lower triangular over `Z`, is invertible over every
> coefficient ring, and loses no labelled information from the weight face.
> This minimal leading flag occupies the consecutive source-normal rows
>
> ```text
> ceil(M/2), ceil(M/2)+1, ..., floor(2M/3).                (7)
> ```

For `M=12`, this gives an exact three-channel coordinate system for the
coefficient walls in the current planar-Jacobian frontier.

The inheritance pass is:

- closest mechanism: THM-3997's all-row Hasse transform;
- canonical hostile: two distinct weight-twelve monomials have different
  three-channel packets, so the scalar weight is not a face coordinate;
- corrected near miss: finite source observers do not themselves provide a
  Darboux lift or upper-bound residual weight;
- least-used sidecar: the invariant diagonal `2n-d` of a source-normal
  coefficient `[x^d t^n]`.

The concept board was

```text
{weighted face, source-normal diagonal, binomial flag,
 integral inverse, M12 walls, seam entry}.                         (8)
```

## 2. Exact monomial expansion and diagonal isolation

For every `i,j>=0`, substitution `(2)` gives

```text
p^i y^j
 =x^j t^(i+2j)(1+x^2t)^(i+j)
 =sum_q binom(i+j,q)x^(j+2q)t^(i+2j+q).                   (9)
```

Every term `[x^d t^n]` on the right satisfies

```text
2n-d=2(i+2j+q)-(j+2q)=2i+3j.                            (10)
```

Consequently the diagonal `2n-d=M` receives exactly the weight-`M` face and
no other residual weight. This proves the face-isolation predicate before any
matrix calculation.

For the ordered monomial `(i_r,j_r)`, its contribution reaches the channel
indexed by `s` when `q=s-r`. Also

```text
i_r+j_r=(M-j_r)/2.                                       (11)
```

Equations `(9)--(11)` give the forward formula `(5)`. Its diagonal entries
are one, proving unimodularity. This is stronger than nonsingularity over the
characteristic-zero seam field: the inverse is defined integrally over `Z`.

## 3. Generating-function proof of the inverse

Put

```text
z=x^2t,                    L=(M-j_0)/2,
F(w)=sum_r c_r w^r.                                         (12)
```

Extend the diagonal coefficients beyond the minimal flag by writing
`hhat_s` for the coefficient at row `ceil(M/2)+s`, `0<=s<=L`, and put
`Hhat(z)=sum_(s=0)^L hhat_s z^s`. The first values `hhat_s` are the `h_s`
in `(4)`. Since `(M-j_r)/2=L-r`, the complete face contribution is

```text
Hhat(z)=sum_r c_r z^r(1+z)^(L-r)
       =(1+z)^L F(z/(1+z)).                               (13)
```

The inverse fractional-linear substitution is

```text
F(w)=(1-w)^L Hhat(w/(1-w)).                               (14)
```

Extracting coefficients in `(14)` gives `(6)`. Thus the forward and inverse
maps are exact mutually inverse natural transformations of coefficient
modules, rather than a finite-field rank observation.

## 4. The minimal lossless row interval

The first admissible `j` is `j_0 in {0,1}`, so the first row in `(4)` is

```text
(M+j_0)/2=ceil(M/2).                                      (15)
```

Successive `j_r` differ by two, hence successive rows differ by one. The
largest admissible `j` has `i>=0` and gives

```text
(M+j_max)/2=floor(2M/3).                                  (16)
```

This proves `(7)` for the minimal independent flag. The full transformed face
can continue on the same diagonal through row `M`, but `(13)--(14)` make
every later coefficient a redundant integral combination of
`(h_0,...,h_m)`. The numerical weight `M` is only the address of a flag; the
lossless object is its ordered coefficient vector.

## 5. Exact `M=12` wall coordinates

The weight-twelve face is

```text
R_12=U p^6+W p^3y^2+Z y^4.                               (17)
```

Equations `(5)--(6)` become

```text
(h0,h1,h2)=(U,6U+W,15U+5W+Z),                            (18)
(U,W,Z)=(h0,h1-6h0,h2-5h1+15h0).                         (19)
```

The transform matrices are

```text
[1 0 0]          [ 1  0 0]
[6 1 0]    and   [-6  1 0].                              (20)
[15 5 1]         [15 -5 1]
```

Therefore the four coefficient walls are exactly

```text
U=0             iff h0=0,
Z=0             iff h2-5h1+15h0=0,
Lambda=0        iff h2-4h1+10h0=0,
D=0             iff h1^2+8h0h1-24h0^2-4h0h2=0.          (21)
```

The channels `(h0,h1,h2)` occur in `G` at rows `t^6,t^7,t^8`, or after the
standard division `Q=G/t`, at rows `t^5,t^6,t^7`. Hence an entry argument
which stops before the last of these rows cannot, in general, distinguish
`Z` or `Lambda`.

Return now to THM-4297's algebraically closed characteristic-zero exact-`M=12`
setting. There is an exact bridge to its contact geometry. The transverse
coefficient is

```text
A=2U+W=h1-4h0.                                              (22)
```

The contact algebra is `K[b_local]/(b_local^12-Lambda)`. Thus the same three
source-normal channels distinguish twelve simple contacts from the one
length-twelve collision. On `Lambda=0`, equation `(21)` reduces to

```text
D=A^2=(h1-4h0)^2,                                           (23)
```

so `h1-4h0!=0` is exactly the transverse `A_23` gate. The hostile face
`(U,W,Z)=(1,-2,1)`, whose channel vector is `(1,4,6)`, has both
`Lambda=0` and `A=D=0`; the collision channel alone is therefore
insufficient without the transverse sidecar.

## 6. Hostile controls and exact scope

The scalar label `M=12` does not determine the face. The two monomials

```text
p^6                         and                         p^3y^2
```

have respective channel vectors

```text
(1,6,15)                    and                         (0,1,5). (24)
```

This is the minimal witness against treating the natural-number rank of a
face as its coefficient state. The preserved object is the labelled flag;
the scalar weight is a lossy quotient.

The theorem supplies no:

1. actual source-normal pair `(A,C)` satisfying the Poisson bracket;
2. converse from prescribed rows to membership in the planar Keller atlas;
3. bound on the highest nonzero residual weight;
4. bridge from source-row depth to raw graph-conductor descent;
5. exact-`M=12` seam entry, another-cell exclusion, or `JC(2)` conclusion.

Its cheapest next use is precise: construct the actual bracket/Hasse rows
through `G`-row eight in the fixed THM-3992 gauge, transform them by
`(18)--(21)`, and ask whether one of the three still-open walls is forced.
The `D=0` equation is quadratic in the visible channels and should be kept as
such rather than replaced by a pairwise orientation or scalar ordinal.

## 7. Reproduction

Run

```bash
python3 -B 04-computation/jc2_weighted_face_source_normal_visibility_thm4298.py
python3 -B -O 04-computation/jc2_weighted_face_source_normal_visibility_thm4298.py
python3 -B 04-computation/jc2_weighted_face_source_normal_visibility_independent_audit_thm4298.py
python3 -B -O 04-computation/jc2_weighted_face_source_normal_visibility_independent_audit_thm4298.py
```

The primary path uses only the Python standard library. The independent path
reconstructs `(13)--(14)` symbolically in SymPy. **QED.**
