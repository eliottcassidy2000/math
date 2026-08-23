---
id: THM-3880
title: "Opposite marked-root signs on one normalization fibre obstruct descent"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  On the normalization of any irreducible
  marked component, a regular square root z of 1+(2/3)AC fixes one global
  sign u=epsilon*z^3 for every polynomial carrier.  Two addresses over the
  same (A,C) point with z-values z0 and -z0, z0!=0, therefore make descent
  impossible.  The undivided u argument includes A=0; z=0 is a genuine
  silent boundary, and same-sign conductor-jet descent remains open.
source: jc_sparse_direct_search / post-THM-3876 intrinsic sign-monodromy abstraction, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The assertion-free symbolic companion
  checks the cusp identity, both global-sign carrier formulas, their exact
  opposite-sign jumps, the THM-3876 primitive-root specialization, the
  undivided A=0 seam, and positive z=0/A=0 controls.  Normal and optimized
  runs byte-match the frozen 14-gate transcript.  Independent hostile audit
  remains required.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3876-monomial-marked-root-normalization-nondescent
script: 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
output: 05-knowledge/results/jc2_marked_root_opposite_sign_nondescent_thm3880.out
script_sha256: 93c911547f9de56f32058b6fd947345d49831b7adbdf6fe4e0e4174879f87566
output_sha256: 6cd51c5840517f469c71525cd3a463b0c5f02618db90f7da9442c02c3c49d4fc
semantic_sha256: e6c43877736a4ee1c1451f8670aabc51ebdde5397278fc8de1ab8aaaa3c10ea2
hash_basis: raw LF bytes
---

# THM-3880 -- opposite marked-root signs cannot come from one polynomial carrier

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  For `b in k[A,C]`, put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b,                       (1)

P=1+(2/3)AC,                         u=1+AC+A^2b.               (2)
```

Let `Gamma subset A2_(A,C)` be an irreducible reduced affine curve and

```text
nu: X -> Gamma                                                    (3)
```

its normalization.  We use the same letters `A,C` for their pullbacks to
the integral normal curve `X`.  Assume that the marked square has a regular
root

```text
z in O(X),                         z^2=P.                       (4)
```

Suppose there are two distinct geometric points `p,q in X(k)` with

```text
nu(p)=nu(q),                 z(q)=-z(p),                 z(p)!=0. (5)
```

Then there is **no** polynomial `b(A,C)` for which `Delta_b` vanishes on
`Gamma`.

The boundary statement is sharp and typed:

- the obstruction remains valid when the common point has `A=0`; one must
  use `u`, not divide by `A^2`;
- when `z(p)=0`, the two signs coincide and the argument is silent; the
  hyperbola `2AC+3=0` supplies an actual polynomial-carrier control;
- when two addresses have the same `z`-value, equality of point values is
  automatic, but descent can still fail at higher conductor/contact jets.
  This theorem asserts no sufficiency in that same-sign case.

Thus THM-3880 closes every **opposite nonzero sign** normalization fibre,
not every nonnormal marked component.

## 1. One global sign on an integral normalization

The universal cusp identity is

```text
A^2 Delta_b=27(P^3-u^2).                                      (6)
```

Assume for contradiction that `Delta_b|_Gamma=0`, and write `B=nu^*b`.
Substitution of `(4)` into `(6)` gives the identity in the domain `O(X)`

```text
u^2=z^6.                                                       (7)
```

Therefore

```text
(u-z^3)(u+z^3)=0,
```

and integrality of `X` forces one sign globally:

```text
u=epsilon z^3,                         epsilon in {+1,-1}.      (8)
```

This is a global sign, not a choice made separately at each normalization
address.  No UFD or factoriality of `O(X)` is needed.

There is also an intrinsic way to recover the marked root from a putative
carrier.  If `P` is not the zero function, then `z_b=u/P` in `k(X)` obeys
`z_b^2=P`.  At a point where `P` has order `d`, equation `(7)` gives
`2 ord(u)=3d`; hence `d` is even and

```text
ord(z_b)=ord(u)-ord(P)=d/2>=0.                                 (9)
```

Thus `z_b` is regular on the normalization.  Any prescribed regular root
`z` in `(4)` differs from it by one global sign.  This explains why `(8)`
is intrinsic rather than an artifact of selecting a square root.

## 2. The fibre contradiction without division

Because `p` and `q` map to the same point of `Gamma`, every function pulled
back from `k[A,C]` has equal values there.  In particular

```text
A(p)=A(q),              C(p)=C(q),              B(p)=B(q),
u(p)=u(q).                                                       (10)
```

But `(5)` and `(8)` give

```text
u(q)=epsilon z(q)^3=-epsilon z(p)^3=-u(p).                     (11)
```

Equations `(10)--(11)` imply `2u(p)=0`.  Characteristic zero and `(8)` then
give `z(p)=0`, contradicting `(5)`.  This proves the theorem.

Notice that this proof never divides by `A`.  If the common point lies on
`A=0`, then `P=1` and `z(p),z(q)` are `+1,-1`; equation `(11)` still gives
opposite nonzero `u`-values.  This is precisely the seam that a formula
with denominator `A^2` would obscure.

If `Gamma` itself is the vertical line `A=0`, then `z^2=1` in the domain
`O(X)`, so `z` is globally `+1` or globally `-1`; the hostile configuration
`(5)` cannot occur.  The theorem therefore has no hidden vertical-line
exception.

## 3. Exact carrier formulas away from A=0

On the open set `A!=0`, equation `(8)` determines the carrier value.  Since

```text
AC=(3/2)(z^2-1),                                               (12)
```

the two global signs give

```text
A^2B_+(z)= z^3-1-AC
          =(z-1)^2(2z+1)/2,                                   (13)

A^2B_-(z)=-z^3-1-AC
          =-(z+1)^2(2z-1)/2.                                  (14)
```

At two opposite addresses, with the order chosen as `z -> -z`, these obey

```text
B_epsilon(-z)-B_epsilon(z)
 =-2 epsilon z^3/A^2.                                         (15)
```

Thus `(15)` is nonzero whenever `A z!=0`.  Equation `(11)` is the regular,
denominator-free extension of the same obstruction across `A=0`.

## 4. The z=0 and same-sign boundaries

At `z=0`, one has

```text
P=0,                         2AC+3=0,                         (16)
```

and both signs in `(8)` give `u=0`.  The jump `(15)` vanishes.  This is a
genuine positive boundary rather than a missing proof.  On the irreducible
hyperbola `2AC+3=0`, take

```text
b_0=2C^2/9.                                                     (17)
```

Direct calculation gives

```text
Delta_(b_0)=-(C^2/3)(2AC+3)^2,                                (18)
```

so the zero-root component really is carried by a polynomial profile.

Likewise `A=0` is not forbidden in the absence of a sign flip.  The profile

```text
b_A=C^2/6
```

satisfies

```text
Delta_(b_A)=-(A C^3/4)(3AC+4),                                (19)
```

and therefore contains the vertical component `A=0` with constant marked
root `z=1`.  What Section 2 excludes is specifically two **opposite** roots
over one point of that component.

Finally, if `z(p)=z(q)`, formulas `(8)` and `(13)--(14)` assign equal point
values.  This is only a zeroth-order condition.  Membership in the
nonnormal branch ring can also impose agreement of conductor jets, so
same-sign fibres and injective cuspidal parametrizations remain open to a
higher-contact analysis.

## 5. THM-3876 is exactly the sign-flip specialization

In THM-3876, after gcd reduction the normalization data are

```text
A=r^M,
C=6r^N(1+r^(M+N)),
z=1+2r^(M+N).                                                  (20)
```

For `M>=3`, choose a primitive `M`th root `zeta`, put

```text
eta=zeta^N,
U=r^(M+N)=-1/(eta+1).                                         (21)
```

At the two colliding addresses `r` and `zeta r`, the marked roots are

```text
z(r)=1+2U=(eta-1)/(eta+1),
z(zeta r)=1+2eta U=-z(r).                                     (22)
```

The root in `(22)` is nonzero because `eta` has order `M>=3`.  Therefore
THM-3880 immediately gives the THM-3876 non-descent.  Moreover, using
`U^2=A^2 r^(2N)`, the plus-sign jump `(15)` becomes

```text
B(zeta r)-B(r)
 =-2r^(2N)(eta-1)^3/(eta+1),                                 (23)
```

exactly the primitive-root difference computed there.  The monomial proof
is thus not a separate coincidence: it is the first explicit realization
of the intrinsic sign-monodromy gate.

## 6. Exact replay and scope

Run

```text
python3 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
python3 -O 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
```

and compare both streams byte-for-byte with the frozen output.  The exact
companion checks `(6)`, both carrier formulas and jumps, the specialization
`(21)--(23)`, and the positive boundary factorizations `(18)--(19)`.

The theorem is a necessary non-descent criterion.  It does not assert that
same-sign point values suffice for branch-ring membership, does not close
conductor jets, and does not classify marked curves for which `P` has no
global square root before normalization.  Those are the precise remaining
directions after the opposite-sign fibres are removed.
