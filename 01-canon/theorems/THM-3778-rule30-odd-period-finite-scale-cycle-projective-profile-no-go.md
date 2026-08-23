---
id: THM-3778
title: "Rule 30 odd-period finite-scale-cycle projective-profile no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every odd
  spatial period n and finite scale-cycle length r, the fully saturated
  algebraic projective-profile cycle locus is a finite union of open subsets
  of projective eigenspaces.  It has no elliptic component.  For the physical
  Rule 30 profiles of THM-3512 over Q_2, the determinant-valuation invoice
  forces the constant mode and consecutive gap one, so no such exact cycle
  exists.  This closes only the stated periodic-profile ansatz.  It proves no
  signalizer, bounded-gap, center, density, balance, or other Rule 30 prize
  statement; every Rule 30 prize remains open.
source: root + lrc14_cover_defect_bridge / 2026-08-23
audit: >
  PASS.  An independent hostile audit repaired the period-three status,
  homogeneous saturation boundary, and platform replay claim, then rederived
  the general odd-period product law, typed twisted composition, projective
  eigenspace bijection, Q_2 holonomy collapse, determinant-unit filter, and
  physical gap invoice.  Even period two, primitive holonomy, finite-field
  component censuses, composition orientation, and coordinate-dependent
  normalization were explicit hostiles.  Normal and optimized runs agree,
  their line-normalized content equals the frozen output, and the companion
  has no Python assert gate.
depends_on:
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
related:
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
script: 04-computation/rule30_odd_period_finite_scale_cycle_no_go_thm3778.py
output: 05-knowledge/results/rule30_odd_period_finite_scale_cycle_no_go_thm3778.out
script_sha256: a2142833f11cdaef46aec8f72ef7657e815742680976b8aaf21c27113758fc4d
output_sha256: 20a4849e69c00f323c0737dd45cfd4f453f3b9526ed8c251f931060cb4a3e9dd
semantic_sha256: 71adc0d815681eddae05c5bc7c74087820b7fbe83e6c324d0d4fbc7c7a391022
hash_basis: raw LF bytes
---

# THM-3778 -- odd-period finite scale-cycles linearize and fail physically

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
closes exact finite scale-cycles inside the odd-period cyclic
projective-profile ansatz.  Its primary proved dependency is THM-3512.
THM-3511 and THM-3516 type the discarded owner/section/carry boundary but are
not used in the algebraic classification.  Every Rule 30 prize remains
**OPEN**.

## 1. Algebraic statement

Let `K` be an algebraically closed field of characteristic zero.  Fix an odd
integer `n>=1`.  For a labeled `n`-periodic profile `g=(g_j)`, define

```text
R(g)_j=-g_(2j)g_(2j+1)(1-g_(2j+2))/(1-g_(2j)),       (1)
```

with indices modulo `n`.  Fix `r>=1`, put `g^(k)=R^k(g)`, and call the orbit
**fully r-saturated** when

```text
g^(k)_j != 0,1       for every 0<=k<r and every j.    (2)
```

Thus every one of the first `r` applications of `R` is defined, every ratio
amplitude is nonzero, and no sibling sum used below vanishes.

For `q in K^*`, let `V_q` be the space of bi-infinite sequences satisfying
`A_(j+n)=qA_j`, represented by `(A_0,...,A_(n-1))`.  Define

```text
T_q:V_q -> V_(q^2),
(T_q A)_j=A_(2j)+A_(2j+1),                            (3)
```

where the right side uses the extension `A_(kn+a)=q^kA_a`.  Operators in a
composition act from right to left.  Put

```text
C_(0,h)=I,
C_(k,h)=T_(h^(2^(k-1))) ... T_(h^2)T_h,              (4)
```

so `C_(k,h):V_h -> V_(h^(2^k))`.

### Theorem A (finite scale-cycle eigenspace classification)

Let `N=2^r-1`.  The reduced fully saturated locus of labeled solutions of
`R^r(g)=g` is the disjoint union, over `h^N=1`, of the images of

```text
P ker(C_(r,h)-lambda I)
```

as `lambda` ranges over the eigenvalues of `C_(r,h)`, after deleting every
hyperplane

```text
(C_(k,h)A)_j=0,       0<=k<=r, 0<=j<n.               (5)
```

The map from an allowed projective amplitude to its profile is

```text
g_j=-A_(j+1)/A_j.                                    (6)
```

It is one-to-one for fixed `h`; its inverse is the recursive amplitude lift,
unique up to common scalar.  Consequently the reduced saturated cycle locus
is a finite union of Zariski-open subsets of projective linear spaces.  In
particular it has no elliptic irreducible component.  Any one-dimensional
irreducible component is an open subset of `P^1`.

This is a statement about labeled profiles.  No quotient by phase rotation,
owner, or Mealy equivalence is taken.

## 2. Proof of Theorem A

Let

```text
p_k=product_j g^(k)_j.
```

Since doubling permutes the indices modulo odd `n`, multiplying (1) cancels
the `1-g` factors and gives

```text
p_(k+1)=-p_k^2,       p_k=-p_0^(2^k) for k>=1.        (7)
```

Full saturation makes `p_0` nonzero.  If `g^(r)=g`, then

```text
p_0^(2^r-1)=-1.                                      (8)
```

Choose `A_0!=0` and recursively impose (6).  Because `n` is odd,

```text
A_n=(-1)^n p_0 A_0=-p_0A_0=hA_0,       h=-p_0.       (9)
```

Equations (8)--(9) give `h^N=1`.

Let `B=T_hA`.  From (6),

```text
B_j=A_(2j)(1-g_(2j)),
B_(j+1)=A_(2j+2)(1-g_(2j+2)),
A_(2j+2)/A_(2j)=g_(2j)g_(2j+1).
```

Therefore

```text
-B_(j+1)/B_j=R(g)_j.                                 (10)
```

Also `B_(j+n)=h^2B_j`, proving the target typing in (3).  Induction gives

```text
g^(k)_j=-(C_(k,h)A)_(j+1)/(C_(k,h)A)_j.              (11)
```

Because doubling is a permutation modulo `n`, the conditions
`g^(k)_j!=0,1` are equivalent to nonvanishing of all coordinates of
`C_(k,h)A` and `C_(k+1,h)A`.  Thus (5) is exactly the homogeneous saturation
boundary; it is not a product after rational cancellations.

Since `h^N=1`, the final twist is

```text
h^(2^r)=h.
```

If `g^(r)=g`, equality of the adjacent ratios in (11) forces the coordinate
quotients `(C_(r,h)A)_j/A_j` to be one common nonzero scalar `lambda`.  Hence
`C_(r,h)A=lambda A`.  Conversely, an eigenline satisfying (5) makes every
intermediate ratio defined and (11) returns to (6).  This proves the stated
bijection.

For fixed `r`, the equation `h^N=1` has finitely many roots.  For each one,
the eigenvector set of the fixed linear operator `C_(r,h)` is the finite
union of its ordinary eigenspaces.  Removing (5) proves the geometric
conclusion.  No diagonalizability assumption is used: generalized
eigenvectors do not solve the required equation unless they are ordinary
eigenvectors.

## 3. Physical statement

Retain THM-3512's exact Rule 30 quantities

```text
q_m=2^m,
G_m(s)=-U_m(s+q_m)/U_m(s),
U_m(s)+U_m(s+q_m)=2^(d_m)U_(m+1)(s),                 (12)
```

where every `U_m(s)` is an odd `Q_2` unit, every `d_m` is a positive integer,
and `d_m` is independent of phase `s`.  Fix a phase `t` and define the full
translation-ray profile

```text
g^(0)_j=G_m(t+jq_m),          j in Z.                 (13)
```

Assume (13) is periodic with period dividing an odd `n`, and assume its exact
projective evolution closes after `r>=1` scales:

```text
R^r(g^(0))=g^(0).                                    (14)
```

Physical profiles automatically satisfy full saturation, because `G_m` is a
nonzero odd unit and `1-G_m=2^(d_m)` times a nonzero unit.

### Theorem B (physical odd-period finite-cycle no-go)

No profile satisfying (12)--(14) exists, for any odd `n>=1` and finite
`r>=1`.

This theorem excludes eventual scale-periodicity as well: if an exact cycle
began at a later scale, apply the statement at that scale.  It does not
exclude spatially nonperiodic profiles or nonclosing scale evolution.

## 4. Proof of Theorem B

Use the actual normalized amplitudes

```text
A^(k)_j=U_(m+k)(t+jq_(m+k)).                          (15)
```

Equation (12), with phase `t+jq_(m+k+1)`, gives the exactly oriented lift

```text
T_(h^(2^k)) A^(k)=2^(d_(m+k))A^(k+1).                (16)
```

The initial holonomy `h=A^(0)_(j+n)/A^(0)_j` lies in `Q_2^*`; Theorem A gives
`h^N=1` for the odd integer `N=2^r-1`.  The only odd-order root of unity in
`Q_2` is `1`.  Indeed, a nontrivial such unit can be written
`h=1+2^a u`, with `a>=1` and `u` odd.  For odd `N`, the first binomial term in
`h^N-1` has valuation `a` and every later term has larger valuation, so the
sum cannot vanish.  Hence

```text
h=1,
C_(r,h)=T_n^r,                                        (17)
```

where

```text
(T_nA)_j=A_(2j)+A_(2j+1)       (indices modulo n).   (18)
```

Let `S` be cyclic shift, `(SA)_j=A_(j+1)`, and let `P_2` be the permutation
matrix `(P_2A)_j=A_(2j)`.  Oddness of `n` makes `P_2` a permutation, and

```text
T_n=P_2(I+S),
|det T_n|=det(I+S)=2.                                 (19)
```

For the last equality, if `zeta` ranges over the `n`th roots of unity, then
`product_zeta(1+zeta)=2` for odd `n`, obtained by evaluating `x^n-1` at
`x=-1`.

The constant vector is an eigenvector of `T_n` with eigenvalue `2`.  Split
the monic integral characteristic polynomial in `Qbar_2`.  Every eigenvalue
is a 2-adic algebraic integer and has nonnegative valuation.  Their product
has valuation `nu_2(det T_n)=1`; the displayed eigenvalue `2` already
contributes one.  It follows simultaneously that:

```text
2 is a simple eigenvalue of T_n,
every other eigenvalue of T_n is a 2-adic unit.        (20)
```

Thus the only eigenvalue of `T_n^r` having positive valuation is `2^r`, and
its eigenspace is the constant line.

Iterating (16) gives

```text
T_n^r A^(0)=2^D A^(r),
D=sum_(k=0)^(r-1)d_(m+k).                             (21)
```

The projective return (14) makes `A^(r)=uA^(0)` for a common scalar `u`.
Both vectors have odd-unit coordinates, so `u` is an odd unit.  Hence the
eigenvalue in (21) satisfies the load-bearing physical invoice

```text
nu_2(lambda)=D>=r.                                   (22)
```

Equations (20)--(22) force `lambda=2^r`, `D=r`, and the constant amplitude
line.  Therefore every gap in the cycle equals one and the profile is
`g_j=-1`.  If `r>=2`, the cycle itself contains consecutive unit gaps.  If
`r=1`, `R(g)=g` repeats the unit gap at the next scale.  Both contradict
THM-3512's proved no-consecutive-gap-one law.  This proves Theorem B.

## 5. Period-three specialization and inherited fixed-profile lane

For `n=3,r=2`, Theorem A gives `h^3=1` and

```text
C_(2,h)=T_(h^2)T_h.
```

At `h=1`, the eigenvalues are `4,1,1`: the `lambda=1` plane gives the open
rational involution `t -> -t`, while the `lambda=4` line is constant.  At
primitive cubic `h`, all three eigenlines meet the homogeneous saturation
boundary.  Equation (22) excludes the rational plane; `lambda=4` forces the
gap pair `(1,1)`.  This recovers the independently audited period-three
classification and no-go.

For `r=1`, Theorem A forces `h=1` and uses `T_n`.  At `n=5,7`, its spectra
are exactly the fixed-profile spectra recorded in the 2026-08-21 reflection.
Theorems A--B show that increasing `n` or replacing fixedness by any exact
finite scale-cycle does not create an elliptic or physical lane inside this
same cyclic-profile ansatz.

## 6. Boundaries and hostiles

- **Even spatial period:** excluded.  Doubling is not a permutation modulo
  `n`; the product cancellation, `P_2`, and determinant argument change type.
- **Unsaturated boundary:** excluded.  A zero amplitude or sibling sum makes
  the ratio lift or an application of `R` undefined.  Simplifying rational
  products can hide a zero times a pole; the homogeneous hyperplanes (5) are
  the authoritative boundary.
- **Characteristic two or positive characteristic:** excluded from Theorem
  A as stated.  Eigenvalue collisions and the root-of-unity count can change.
- **Algebraic versus physical:** primitive holonomies exist over algebraic
  extensions and are part of Theorem A.  The conclusion `h=1`, the valuation
  invoice, and the determinant/unit filter belong only to Theorem B over
  `Q_2` with actual normalized Rule 30 amplitudes.
- **Phase-dependent normalization:** excluded.  Equation (22) uses the
  inherited fact that one positive `d_m` normalizes every sibling pair at a
  fixed scale.  An ambient coordinate-dependent odd-unit system is not a
  physical control.
- **Approximate, nonclosing, or infinite scale evolution:** untouched.
- **Owners and carries:** the projective profile discards THM-3511's phase
  owner and THM-3516's signed gauge/carry.  Since Theorem B is a necessary-
  condition no-go, it need not reconstruct them and makes no sufficiency
  claim for any surviving abstract eigenline.
- **Minimal periods:** not assumed.  “Period `n`” means period dividing `n`;
  the constant spatial-period-one boundary is intentionally retained.

## 7. Verification and strict prize scope

The assertion-free companion
`04-computation/rule30_odd_period_finite_scale_cycle_no_go_thm3778.py` checks:

- the product law at odd periods `1,3,5`;
- the `n=2` hostile, where the product law changes and `det T_2=0`;
- the rightmost-first `T_h` orientation, including an independent period-five
  nonclosing holonomy control at `h=13`;
- the complete period-three ordinary and primitive eigenspaces and exact
  homogeneous boundary;
- exhaustive saturated finite-field hostiles at `p=5,7,11,13`;
- the inherited fixed spectra at periods five and seven;
- `|det T_n|=2` and the one-dimensional constant eigenspace for every odd
  `1<=n<=15`; and
- an ambient odd-unit vector whose sibling valuations are `(1,2,2)`, showing
  why THM-3512's phase-independent physical normalizer is essential.

These are exact controls for the displayed universal proofs, not substitutes
for them.  Normal and optimized runs are identical, and their line-normalized
content equals the frozen output.  Final raw-byte hashes are

```text
a2142833f11cdaef46aec8f72ef7657e815742680976b8aaf21c27113758fc4d
  rule30_odd_period_finite_scale_cycle_no_go_thm3778.py
20a4849e69c00f323c0737dd45cfd4f453f3b9526ed8c251f931060cb4a3e9dd
  rule30_odd_period_finite_scale_cycle_no_go_thm3778.out
```

This draft proves no finite or bounded orbit-signalizer graph, no bound on
the actual innovation gaps beyond excluding this exact periodic ansatz, no
center nonperiodicity, no limiting density, no balance, no prediction lower
bound, and no other Rule 30 prize statement.  Every Rule 30 prize remains
**OPEN**.
