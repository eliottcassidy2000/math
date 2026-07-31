---
id: THM-2950
title: "Three-conjugate-pair V-four torsor and quartic resolvent frame"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Three
  unordered conjugate pairs have eight oriented
  half-systems and four classes modulo global conjugation.  The latter
  form a V_4 torsor; pair permutations act as GL(2,2)=S_3, and the
  resulting affine action is AGL(2,2)=V_4 semidirect S_3=S_4.  The
  four half-norm traces therefore form a canonical quartic resolvent
  frame.  Its three pairing coordinates carry the S_3 quotient, with
  kernel V_4, and the quartic and cubic-resolvent discriminants agree
  identically.  The full nonnegative norm is common to all four
  half-systems and cannot recover their phase.  No canonical C_3
  ordering, modular-group action, SFC(4) closure, or Jacobian
  conclusion is claimed.
source: codex-gmc-modular-two-three-synthesis-2026-07-29
audit: >
  An independent hostile audit rederived the intrinsic torsor rather
  than based-vector-space typing, the signed-pair permutation quotient,
  AGL(2,2)=S4 action, trace invariance, V4 kernel and S3 image on
  labeled pairings, all three resolvent-difference signs, and the
  discriminant identity.  It verified the norm/phase hostiles,
  recomputed both LF hashes, and replayed normal, optimized, and stored
  output byte-for-byte.  The final wording explicitly allows extra
  information loss when numerical resolvent roots collide.
depends_on:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
related:
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
script: 04-computation/gmc_three_pair_v4_quartic_resolvent_thm2950.py
output: 05-knowledge/results/gmc_three_pair_v4_quartic_resolvent_thm2950.out
script_sha256: 651319fc197e9c3c1c20599fa9f99f5916dc386e5c61a37f255ece3d1c64f89c
output_sha256: 443790e8f19fdbcf6bf9e80730c1f4e2d69a08a64269693ddd354af0cae5bf50
hash_basis: LF-normalized bytes
---

# THM-2950 -- three-conjugate-pair V-four torsor and quartic resolvent frame

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Three pairs produce a four-point affine plane

Let

```text
Z={P_1,Pbar_1,P_2,Pbar_2,P_3,Pbar_3}                 (1)
```

be a reduced finite real scheme consisting of three nonreal conjugate
pairs.  An oriented half-system chooses one member of each pair, so
after choosing temporary representatives `P_i` the eight choices are

```text
epsilon in F_2^3.                                    (2)
```

Global conjugation sends

```text
epsilon |-> epsilon+(1,1,1).                         (3)
```

Therefore the set `Omega` of half-systems modulo global conjugation has
four elements.  Intrinsically `Omega` is a torsor for

```text
V=F_2^3/<(1,1,1)> ~= F_2^2 ~= V_4.                  (4)
```

It is a torsor, not a canonically based vector space: choosing one
half-system as origin supplies the identification in `(4)`.

Permuting the three conjugate pairs acts on `V`, fixes zero after a
base choice, and permutes the three nonzero vectors faithfully.  Hence

```text
S_3 ~= GL(2,2).                                      (5)
```

The complete affine action is consequently

```text
V_4 semidirect S_3
  ~= AGL(2,2)
  ~= S_4,                                            (6)
```

acting on the four points of `Omega`.  Equivalently, the signed
pair-permutation group `(F_2^3) semidirect S_3`, modulo the central
global conjugation in `(3)`, is `S_4`.

This is the exact `2`-by-`3` structure.  The involution is genuine
complex conjugation; the three objects are its pair orbits.  Equation
`(6)` does **not** choose a cyclic order on those three orbits.

## 2. Four half-norm traces

Let `f` be any real function on `Z` with

```text
z_i=f(P_i),                    f(Pbar_i)=zbar_i.       (7)
```

For a half-system `epsilon`, put

```text
h_epsilon
 =product_i (z_i if epsilon_i=0 else zbar_i).         (8)
```

Complementary choices give conjugate products.  Thus each class
`omega=[epsilon] in Omega` has a well-defined real trace

```text
T_omega=h_epsilon+hbar_epsilon=2 Re(h_epsilon).       (9)
```

The four traces are permuted by the `S_4` frame `(6)`.  Hence

```text
Q_f(X)=product_(omega in Omega)(X-T_omega)            (10)
```

is independent of every temporary labelling and orientation used to
write `(7)--(9)`.

By contrast, every half-system has the same product norm

```text
h_epsilon hbar_epsilon
 =product_i |z_i|^2
 =N(f).                                               (11)
```

The nonnegative complete-intersection norm in THM-2843 is exactly of
the form `(11)`.  It is invariant under the whole `V_4` translation
torsor and therefore cannot select a half-system or recover the four
traces.

This failure is sharp even with nonzero values.  The two triples

```text
(z_1,z_2,z_3)=(1,1,1),
(z_1,z_2,z_3)=(i,1,1)                                (12)
```

both have `N=1`, while their four trace vectors are respectively

```text
(2,2,2,2),                    (0,0,0,0).              (13)
```

For the exact Gaussian-integer control

```text
(z_1,z_2,z_3)=(1+i,2+i,3+2i),                        (14)
```

the common norm is `130` and the four traces are

```text
(-6,18,14,22),                                       (15)
```

so the torsor data can be fully separated while the norm remains one
scalar.

## 3. The cubic resolvent is the S_3 quotient

Write the four roots in `(10)` as `r_0,r_1,r_2,r_3`.  The three
unordered pairings of four objects give the standard resolvent roots

```text
u_1=r_0 r_1+r_2 r_3,
u_2=r_0 r_2+r_1 r_3,
u_3=r_0 r_3+r_1 r_2.                                 (16)
```

The action of `S_4` on these three pairings has

```text
kernel V_4,                         image S_3.         (17)
```

Thus the action on the three **labeled** resolvent roots factors
exactly through the `S_3` quotient of `(6)` and forgets the affine
`V_4` position.  If two numerical roots in `(16)` collide, the
specialized cubic can of course forget additional information.

The discriminant identity is immediate from

```text
u_1-u_2=(r_0-r_3)(r_1-r_2),
u_1-u_3=(r_0-r_2)(r_1-r_3),
u_2-u_3=(r_0-r_1)(r_2-r_3).                          (18)
```

Squaring and multiplying gives

```text
disc(product_i(X-r_i))
 =disc(product_j(X-u_j)).                            (19)
```

This is the structural reason a quartic and its cubic resolvent have
the same discriminant: the cubic sees the quotient
`S_4/V_4 ~= S_3`.

## 4. Consequences and stopping boundary

On the reduced branch of THM-2843, `Z=V(Q,C)` has precisely the form
`(1)`.  THM-2947 uses only the conjugate-pair involution to prove
even corank; the present theorem identifies the finer four-valued
phase object which the full norm discards.

For a hypothetical quartic cover, a cubic-resolvent theorem can
constrain:

```text
* the three pairings;
* the S_3 quotient monodromy; and
* the common discriminant.                            (20)
```

It cannot by itself reconstruct:

```text
* a point of the V_4 half-system torsor;
* a chosen half-norm phase; or
* the full quartic cover from its cubic quotient.     (21)
```

The missing sidecar is therefore exact: an oriented half-system, or an
equivalent nonzero `V_4`-covariant observable.  A nonnegative norm is
orthogonal to that task.

Nor does `S_3` in `(5)` give a canonical `C_3` generator.  Selecting a
three-cycle requires an additional orientation.  Consequently the
identity

```text
PSL_2(Z)=C_2 * C_3                                   (22)
```

is an illuminating grammar comparison here, not an action supplied by
the complete intersection.  No modular-group, SFC(4), planar
Jacobian, or Dixmier consequence follows from `(1)--(21)`.

The half-system description also requires the reduced six-point
branch.  THM-2947's nonreduced parity theorem survives when individual
points and half-systems do not.

## 5. Exact companion

The companion

```text
python 04-computation/gmc_three_pair_v4_quartic_resolvent_thm2950.py
python -O 04-computation/gmc_three_pair_v4_quartic_resolvent_thm2950.py
```

checks with integer arithmetic:

1. eight half-systems and four complement classes;
2. the fixed-point-free `V_4` translations;
3. all `24` affine permutations and equality with `S_4`;
4. the `V_4` kernel and `S_3` image on the three pairings;
5. four exact quartic/resolvent discriminant controls; and
6. the phase-loss and phase-separation examples `(12)--(15)`.

Both modes reproduce the stored transcript byte-for-byte.
