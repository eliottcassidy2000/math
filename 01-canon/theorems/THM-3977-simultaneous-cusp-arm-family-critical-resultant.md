---
id: THM-3977
title: "Simultaneous cusp-arm family critical resultant"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Over an
  algebraically closed field of characteristic zero, on the height-two
  exact-volume completion put A_(c,r)=c p+y^2+r x with c nonzero. Its
  restrictions to the added boundary and retained arm are y^2 and c p.
  Every r!=0 row has an affine critical point, detected by an exact
  resultant. The r=0 row is submersive but its rational generic fibre has
  six nonzero logarithmic residues, so it has no rational mate with bracket
  in k(A)^*. Hence no member of this lowest-seam family has a polynomial
  constant-Jacobian mate. No unrestricted boundary family, finite map, or
  counterexample to JC(2) is claimed.
source: jc-degree6-one-place + all-frontiers / post-THM-3973 simultaneous boundary-jet hostile, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (audit-3977-candidate, 2026-08-24). The
  audit independently rederived the generalized-c critical resultant, its
  nonzero finite critical root and stable t-leading rows, the zero-seam
  submersion, the rational generic-fibre parametrization, the scaled
  discriminant, and every residue. It enforced the two scope guards:
  squarefreeness is over k(P), not every scalar fibre, and nonzero-seam
  criticality excludes regular mates only. Normal and optimized 25-gate
  runs match the frozen output after LF normalization; all hashes pass.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3974-height-tower-few-weight-darboux-support-obstruction
  - THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate
script: 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
output: 05-knowledge/results/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.out
script_sha256: 8fc3507c40e31a54117ba96c7f56261f49386a0a61c5769800d3ae5ef56b9451
output_sha256: 1534633377e5722a13a933db89da66a22cb480e22f0f042acd753642c0b78c72
semantic_sha256: 8aadeddb4635f490ea0acde0694b539d41f6145ac66e93df36af5543c8e6ae47
hash_basis: raw LF bytes
---

# THM-3977 -- the lowest simultaneous cusp-arm seam has no Darboux mate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. On the height-two
completion of THM-3973 put

```text
z=1+x^2t,                 p=zt,                 y=xzt^2,
B_2=k[x,z,p,y] subset k[x,t].                              (1)
```

Its added boundary and retained arm are

```text
D=V(x,z,p) isomorphic to A1_y,
L_1=V(x,z-1,y) isomorphic to A1_p.                         (2)
```

For `c in k^*` and `r in k`, define

```text
A_(c,r)=c p+y^2+r x.                                      (3)
```

Then no `A_(c,r)` has a polynomial constant-Jacobian mate in `B_2`.
More precisely:

1. if `r!=0`, `A_(c,r)` has a genuine critical point on the source plane;
2. if `r=0`, `A_(c,0)` is submersive, but there are no
   `C in k(x,t)` and `h in k(A_(c,0))^*` with
   `J(A_(c,0),C)=h(A_(c,0))`.

The two branches use different mechanisms: exact affine elimination for the
nonzero seams and logarithmic residues on a rational generic fibre at the
submersive endpoint.

## 1. The simultaneous boundary data select the seam

The restrictions of `(3)` are exactly

```text
A_(c,r)|D=y^2,                 A_(c,r)|L_1=c p.           (4)
```

Thus the same function has the first cuspidal coordinate on `D` and a
nonconstant affine-linear coordinate on `L_1`. Among generator-linear
corrections, `x` is the only seam direction vanishing on both: a linear
combination of `x,z,p,y,1` that vanishes on `D` first loses its `y` and
constant terms, and vanishing on `L_1` then loses its `z` and `p` terms.
This explains why `(3)` is the minimal simultaneous family rather than an
arbitrary bounded search.

## 2. Every nonzero seam has an affine critical point

Direct differentiation gives

```text
A_x=r+2cxt^2+2xt^4+8x^3t^5+6x^5t^6,
A_t=c+2cx^2t+4x^2t^3+10x^4t^4+6x^6t^5.                 (5)
```

Their exact resultant in `t` is

```text
Res_t(A_x,A_t)=96x^24 H_(c,r)(x),                        (6)

H_(c,r)=486r^5x^12+216cr^4x^9+12c^2r^3x^6-64r^4x^5
        +2c^5rx^4+8c^3r^2x^3+c^6x-2c^4r.               (7)
```

Suppose `r!=0`. The polynomial `H_(c,r)` is nonconstant and

```text
H_(c,r)(0)=-2c^4r!=0.                                    (8)
```

Since `k` is algebraically closed, it has a root `alpha`, and `(8)` makes
`alpha!=0`. At that root the `t`-leading rows remain

```text
deg_t(A_x)=6,       lc_t(A_x)=6alpha^5,
deg_t(A_t)=5,       lc_t(A_t)=6alpha^6.                  (9)
```

Thus the resultant zero represents a finite common `t`-root, not a root at
infinity caused by degree drop. The resulting point lies in `D(x)` on the
source plane and satisfies `A_x=A_t=0`. Any regular `C` has
`J(A_(c,r),C)=0` there, so it cannot have constant bracket one. This proves
the `r!=0` branch. It does not exclude a rational function allowed to have a
pole at that critical point.

## 3. The zero seam is submersive

Set

```text
A=A_(c,0)=c p+y^2.                                       (10)
```

One has

```text
Res_t(A_x,A_t)=96c^6x^25,
J(A,y)=-czt^2.                                           (11)
```

If `A_x=A_t=0`, the second identity forces `t=0` or `z=0`. At `t=0`,
`A_t=c!=0`. At `z=0`, one has `x!=0`, `t=-x^-2`, and direct substitution
gives

```text
A_t=-c,                   A_x=2c/x^3.                    (12)
```

Hence `(10)` has no source-plane critical point. The remaining obstruction
is global exactness on its generic fibre.

## 4. The rational generic fibre has six logarithmic poles

Let `P` be transcendental and work over `K=k(P)` on the generic fibre
`P=A`. Put

```text
q=P-y^2=c p,
E=q^3-c^3y^2=c^3(p^3-y^2).                              (13)
```

The determinantal identity

```text
py+xy^2=xp^3                                             (14)
```

gives the complete rational parametrization

```text
x=c^2yq/E,             z=q^3/E,             t=E/(cq^2). (15)
```

These formulas reconstruct `(1)` and show that the generic function field is
exactly `K(y)`. Under `(15)`, equation `(11)` becomes

```text
J(A,y)=-E/(cq).                                          (16)
```

If some `C in k(x,t)` and `h in k(A)^*` satisfied `J(A,C)=h(A)`, restriction
to the generic fibre and the coordinate rule in `(16)` would give

```text
dC=-h(P)cq/E dy.                                         (17)
```

The pole polynomial has exact discriminant

```text
Disc_y(E)=64P^3c^12(27P^2+4c^3)^2.                      (18)
```

This is nonzero in `K`, so `E` is squarefree in `K[y]`. Moreover
`gcd(E,q)=1`: if `q=0` at an `E`-root, then `y=0` and hence `P=0`, impossible
for the transcendental generic value. Consequently `(17)` has a nonzero
residue

```text
-h(P)cq(a)/E'(a)                                         (19)
```

at each of its six simple roots `a`. The differential of a rational function
has zero residue at every place, contradicting `(17)`. This proves the
stronger rational no-mate statement at `r=0`.

The scalar fibres `P=0` and `27P^2+4c^3=0` are nonsquarefree; they are not
being called squarefree and play no role in the generic-fibre contradiction.

## 5. Preservation, loss, and exact scope

The connection to the adjacent completion results is typed:

```text
THM-3973 preserves  B_2, D, L_1, and the exact boundary addresses;
THM-3974 preserves  the weight grading but gives only few-support nonentry;
THM-3975 preserves  the generic-fibre exact-differential method;
this theorem uses   a rational fibre with six logarithmic poles, not the
                    hyperelliptic equation of THM-3975.                 (20)
```

The result closes only the lowest generator-linear seam family `(3)`. It
does not exclude nonlinear simultaneous corrections, prove finiteness for a
target map, construct a Darboux pair, or settle `JC(2)`.

## 6. Exact companion

The companion checks 25 exact identities and hostile gates: all three
determinantal rows, both boundary restrictions, derivatives, the full
resultant, stable leading coefficients, endpoint submersivity, the generic
parametrization, discriminant, and both generic gcds.

Reproduce with

```bash
python3 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
python3 -O 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
```

Both runs print `CHECKS=25` and, after LF normalization, byte-match the frozen
output. **QED.**
