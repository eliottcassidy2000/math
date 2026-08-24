---
id: THM-3889
title: "Maximally confluent quadratic binary cubics still have two infinity places"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The normalized
  split quadratic binary-cubic
  row whose leading discriminant is the single eighth power C^8 remains a
  two-place trap after every homogeneous-linear perturbation.  All four
  index coefficients vanish at the origin, so the resulting cubic order is
  globally nonmonogenic; nevertheless, whenever its discriminant is
  irreducible, the unique projective infinity point has at least two
  normalization branches.  The obstruction is an exact three-case Newton
  boundary, not a bounded coefficient search.
source: jc_sparse_direct_search / post-THM-3801 nonmonogenic cubic-order search, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-23).  The audit
  rederived every Newton support halfspace and the complete
  `delta/gamma/eta` seam split.  It checked nonzero edge roots including all
  double-root parameter seams, used Newton--Puiseux only over the declared
  algebraically closed characteristic-zero field, and verified that the
  different `x`-orders give distinct completed branches and normalization
  places.  It also checked the reducible last seam, index-form unit gate, and
  conditional normal `S3` conclusion.  Normal and optimized runs byte-match
  the frozen 65,675-check transcript.  The two censuses remain FINITE-EXACT
  side evidence and are not used in the proof.
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
script: 04-computation/jc2_maximally_confluent_quadratic_binary_cubic_thm3889.py
output: 05-knowledge/results/jc2_maximally_confluent_quadratic_binary_cubic_thm3889.out
script_sha256: 077167e4ce237d77524a9f8611b6ae61db38af0e18a7333dc76185bf9be33ed0
output_sha256: 24161fe5b2667f8e4e2eed69df363a3007d64f27842ea51e8a3ee444e9900f19
semantic_sha256: 8cfd47a570c914882e86a317a26e77e211f1d7cf9ec4899a34aaaf488de1597a
hash_basis: raw LF bytes
---

# THM-3889 -- the first confluent quadratic row still pays two places

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of
characteristic zero and put `R=k[A,C]`.  Let

```text
a=A^2+AC+alpha A+alpha_1 C,
b=2AC+C^2+beta A+beta_1 C,
c=C^2+gamma A+gamma_1 C,
d=delta A+eta C,                                                (1)
```

where the eight Greek letters are arbitrary scalars in `k`.  Give

```text
S=R 1 direct_sum R omega direct_sum R theta                    (2)
```

the Delone--Faddeev multiplication law belonging to the binary cubic

```text
f(U,V)=aU^3+bU^2V+cUV^2+dV^3.                                  (3)
```

Its discriminant is

```text
Delta=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.                     (4)
```

Then:

1. the binary index form of `S/R` is `-f(X,Y)`, and it represents no unit
   of `R`; hence `S/R` is globally nonmonogenic;
2. `Delta` has degree eight and its degree-eight part is exactly `C^8`;
   consequently its projective closure meets the line at infinity at the
   single point

   ```text
   P_infinity=[1:0:0],                                         (5)
   ```

   with total intersection multiplicity eight; but
3. if `Delta` is irreducible, the normalization of `V(Delta)` has at least
   two distinct points above `P_infinity`.

Thus the most economical quadratic repair of the homogeneous Veronese row
passes the common-zero nonmonogenicity gate and concentrates all projective
intersection at one point, but it still fails the **one normalization
place** gate.  A single projective support point is not a single place.

If, in addition, the generic Delone--Faddeev algebra is a field, irreducible
`Delta` makes `(2)` a normal finite-flat nonmonogenic `S3` cubic order: the
order is free and hence `S2`, while the squarefree height-one discriminant
has valuation one and forces maximality there.  The theorem says that even
this algebraically positive boundary cannot have the one-place branch
passport required by the current Keller construction lane.

## 1. Why this is the maximally confluent split leading row

The quadratic part of `(3)` factors as

```text
f_2(U,V)
 =U(AU+CV)((A+C)U+CV).                                        (6)
```

Its three roots are split.  Along `C=0` all three coalesce; two root
differences have order one and the third has order two.  Directly,

```text
Disc(f_2)=C^8.                                                (7)
```

This is the first nontrivial way a split quadratic coefficient row can
concentrate all eight discriminant intersections in one target direction.
It is the coefficient-space analogue of the two-place tradeoff in THM-3879:
the support has coalesced, but the normalization addresses have not.

Since every coefficient in `(1)` lies in the maximal ideal `(A,C)`, every
value of `f(X,Y)` also lies in `(A,C)`.  The universal index determinant is

```text
det(1,T,T^2)=-f(X,Y),             T=X omega+Y theta.           (8)
```

It can therefore never be a scalar unit.  This proves assertion 1.  Notice
also the unavoidable singularity tariff

```text
Delta in (A,C)^4.                                             (9)
```

Any reduced discriminant in this grammar has multiplicity at least four at
the coefficient-zero point.  For a rational plane discriminant this already
forces degree at least five, because a reduced multiplicity-four singularity
has delta invariant at least six.

## 2. The exact local polynomial at infinity

Use the chart `A!=0` at `(5)` and write

```text
x=C/A,                         z=Z/A.                         (10)
```

Multiplying the four coefficients by `z^2` gives the regular local row

```text
a_tilde=1+x+z(alpha+alpha_1 x),
b_tilde=2x+x^2+z(beta+beta_1 x),
c_tilde=x^2+z(gamma+gamma_1 x),
d_tilde=z(delta+eta x).                                     (11)
```

Let

```text
H(x,z)=Disc(a_tilde,b_tilde,c_tilde,d_tilde)
      =z^8 Delta(1/z,x/z).                                  (12)
```

Equation `(7)` says

```text
H(x,0)=x^8.                                                  (13)
```

The rest of the proof reads the complete lower Newton boundary of `(12)`.
Only `delta`, `gamma`, `eta`, `beta`, and `gamma_1` occur on that boundary;
the remaining parameters affect strictly higher terms.

## 3. First case: delta is nonzero

If `delta!=0`, every exponent `(i,j)` of a nonzero monomial `x^i z^j` in
`H` satisfies

```text
i+3j>=6,                         i+5j>=8.                    (14)
```

The lower boundary has the two primitive edges

```text
(0,2) -- (3,1) -- (8,0),                                    (15)
```

with universal vertex coefficients

```text
[z^2]H=-27delta^2,       [x^3z]H=4delta,       [x^8]H=1.     (16)
```

On the first edge, putting `z=w x^3` gives

```text
delta w(4-27delta w).                                      (17)
```

Its nonzero simple root gives a normalization branch

```text
z=(4/(27delta))x^3+higher terms.                            (18)
```

On the second edge, putting `z=w x^5` gives

```text
1+4delta w,                                                (19)
```

and therefore the distinct branch

```text
z=-(1/(4delta))x^5+higher terms.                            (20)
```

Their different first exponents make the two normalization places
unmistakably distinct.

## 4. Second case: delta vanishes and gamma is nonzero

Now suppose

```text
delta=0,                         gamma!=0.                    (21)
```

The four persistent boundary coefficients are

```text
[z^3]H=-4gamma^3,
[x^2z^2]H=-27eta^2+36eta gamma-8gamma^2,
[x^4z]H=4(eta-gamma),
[x^8]H=1.                                                     (22)
```

### 4.1. The unequal seam `eta!=gamma`

Every monomial satisfies

```text
i+2j>=6,                         i+4j>=8.                    (23)
```

The first edge `(0,3)--(4,1)` has edge polynomial, after
`z=w x^2`,

```text
w[-4gamma^3w^2
  +(-27eta^2+36eta gamma-8gamma^2)w
  +4(eta-gamma)].                                             (24)
```

The bracket is a quadratic with nonzero leading and constant terms, so it
has a nonzero root over `k` and contributes a branch of order two.  The
separate primitive edge `(4,1)--(8,0)` has polynomial

```text
1+4(eta-gamma)w                                               (25)
```

under `z=w x^4`, and contributes a distinct branch of order four.

### 4.2. The equality seam `eta=gamma`

Here `[x^4z]H=0`, but the apparent collision does not join the places.  The
lower boundary becomes

```text
(0,3) -- (2,2) -- (8,0),                                    (26)
```

and

```text
[x^2z^2]H=gamma^2.                                          (27)
```

The first edge gives

```text
gamma^2w^2(1-4gamma w)                                      (28)
```

under `z=w x^2`, hence one branch with leading coefficient
`1/(4gamma)`.  The second edge, under `z=w x^3`, gives

```text
1+2(2beta+gamma-2gamma_1)w+gamma^2w^2.                      (29)
```

Its constant and quadratic coefficients are nonzero, so it has a nonzero
root and supplies a second, order-three branch.  This closes every seam in
`(21)`.

## 5. The last coefficient seam is reducible

It remains to set

```text
delta=gamma=0.                                               (30)
```

Then

```text
c=C(C+gamma_1),                         d=eta C.             (31)
```

Every term in `(4)` is divisible by `C`, so

```text
C divides Delta.                                             (32)
```

Because the degree-eight term is `C^8`, the quotient is nonconstant.
Therefore `Delta` is reducible.  In particular an irreducible discriminant
must lie in Section 3 or Section 4, where two distinct infinity places have
already been constructed.  This proves assertion 3.

## 6. Finite exact side evidence and scope

The companion contains two bounded hostile checks which are deliberately
separate from the proof.

First, each of `a,b,c,d` is independently chosen from the sixteen forms

```text
{A,C,A+C,A-C}+{0,A^2,AC,C^2}.                                (33)
```

Among all `16^4=65,536` labelled rows, exactly `5,890` have a leading
discriminant which is a pure coordinate power.  Exactly `92` of those have a
raw lower Newton chain consisting of a single lattice-primitive edge at the
unique infinity point, but in every
one of the `92` rows the full discriminant is divisible by the transverse
coordinate.  Thus all `92` are reducible and there are zero
irreducible-eligible primitive-edge survivors.  This is a **FINITE-EXACT
census**, not a theorem about arbitrary quadratic rows.

Second, in the all-parameter family `(1)`, restrict each linear coefficient
to `{A,C,A+C,A-C}`.  Of the `4^4=256` rows, `192` lie in `delta!=0`, `48`
lie in `delta=0,gamma!=0`, and the remaining `16` have reducible
discriminant by `(32)`.  Thus this bounded universe has `240` certified
two-place rows and no one-place survivor.

THM-3889 does not classify every degree-two binary-cubic row.  It closes the
first maximally confluent **split leading** grammar, where the leading
discriminant already pays only one projective point.  Nonsplit quadratic
leading rows, higher coefficient degree, discriminants with a different
Puiseux collision, and every subsequent affine-atlas/Jacobian condition
remain open.  The result supplies no planar Jacobian counterexample and no
general Jacobian obstruction.
