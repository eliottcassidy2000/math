---
id: THM-3067
title: "Tetrahedral modular two-three flag quotient and origin loss"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  After choosing
  an affine origin and a cyclic orientation of the
  three nonzero directions of F_2^2, its twelve directed point-direction
  flags form a regular A_4 bitorsor.  Right edge reversal S and right
  direction rotation R satisfy S^2=R^3=(SR)^3=1 and give the tetrahedral
  quotient of C_2*C_3.  Forgetting the point has exact kernel V_4 and image
  C_3.  Forgetting edge orientation gives six K_4 edges, but R does not
  descend.  An odd affine reflection fixes S and conjugates R to R^-1; it is
  an additional sidecar, not the modular involution already used as a V_4
  translation.  No physical LRC/JC intertwiner or closure is claimed.
source: root-2026-08-01-modular-two-three-tetrahedral-flags
audit: >
  An independent read-only hostile audit rederived the right-action
  convention and tetrahedral relation, the regular A4 identification, the
  C3 direction image and V4 kernel, the six-edge no-descent hostile, the odd
  reflection laws, and the complete fixed-R involution census.  It replayed
  normal and optimized modes against the stored eleven-line transcript and
  independently matched both LF hashes.  The final wording distinguishes
  the induced-action kernel from the set projection and separates this A4
  translation quotient from THM-2595's linear-reflection S3/S4 lifts.
depends_on:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
  - THM-3049-k4-matching-monomial-tropical-root-extraction-clutch
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
external:
  - "Mario Krenn, Xuemei Gu, and Daniel Soltesz, Questions on the Structure of Perfect Matchings Inspired by Quantum Physics, arXiv:1902.06023v2."
script: 04-computation/modular_tetrahedral_flag_bitorsor_thm3067.py
output: 05-knowledge/results/modular_tetrahedral_flag_bitorsor_thm3067.out
script_sha256: faa8205afbab9aebf48c2be9aadfa59ef193997380d138588dad7d35f6af40b0
output_sha256: 7d90364b2116766a9e248376c1057a98a0e96f422cb2bd58a005f79276222cc1
hash_basis: LF-normalized bytes
---

# THM-3067 -- tetrahedral modular two-three flag quotient and origin loss

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The twelve flags

Let

```text
V=F_2^2,                    D=V\{0},
Omega=V x D.                                                (1)
```

The flag `(x,d)` is equivalently the directed edge

```text
x -> x+d                                                    (2)
```

of the tetrahedron on the four points of `V`.  Thus `Omega` has twelve
elements.  Choose an order-three linear map

```text
rho in GL(V)=S_3                                             (3)
```

which cycles the three elements of `D`.  Define two right moves on `Omega`:

```text
S(x,d)=(x+d,d),                 R(x,d)=(x,rho d).            (4)
```

Here `S` reverses one directed edge and `R` cycles the three outgoing
directions while keeping the tail fixed.  These are exactly the binary and
ternary local moves; the choice of `rho` is orientation data.

## 2. The tetrahedral modular quotient

Put

```text
A=V semidirect <rho>.                                       (5)
```

This is the orientation-preserving affine group of `V`, hence

```text
A=AGL^+(2,2)=V_4 semidirect C_3=A_4.                        (6)
```

Fix `d_0 in D`.  The map

```text
(a,rho^k) |-> (a,rho^k d_0)                                 (7)
```

is a bijection `A -> Omega`.  Right multiplication in `A` by

```text
s=(d_0,1),                    r=(0,rho)                      (8)
```

becomes precisely `S,R` in `(4)`.  The conjugates of `s` by powers of `r`
are the three nonzero translations, so `s,r` generate all of `A`.
Consequently the right action is regular.

Directly,

```text
S^2=1,                         R^3=1,                        (9)
```

and

```text
(SR)^3:
  x |-> x+d+rho d+rho^2 d=x,
  d |-> rho^3d=d.                                            (10)
```

Thus the modular presentation

```text
PSL_2(Z)=C_2*C_3                                             (11)
```

has this finite tetrahedral quotient; the extra spherical relation is
`(SR)^3=1`.  Equivalently, the chosen moves realize the standard triangle
group `(2,3,3)=A_4`.  This is a quotient of the modular grammar, not a
faithful action of the free product.

There is simultaneously a left geometric action

```text
(a,rho^k).(x,d)=(a+rho^k x,rho^k d).                        (12)
```

It is also regular and commutes with both right moves.  Hence `Omega` is an
`A_4` bitorsor: left multiplication is geometric tetrahedral relabelling,
while right multiplication is the binary/ternary move grammar.

## 3. Exact point-origin loss

Projection to the direction coordinate gives

```text
pi:Omega -> D,                 pi(x,d)=d.                    (13)
```

On the right action,

```text
pi S=pi,                       pi R=rho pi.                  (14)
```

The induced right action on directions has image `C_3`, and its kernel is
exactly the four translations `V_4`.
Thus

```text
A_4/V_4=C_3.                                                (15)
```

Each direction has four possible point origins.  In particular the modular
order-two generator `S` lies in the lost kernel: after forgetting the point,
the binary free factor becomes invisible rather than turning into the
reflection needed for a full `S_3` direction action.

This is the exact finite information ledger behind the recurring `2/3`
picture.  The ternary quotient is real, but it does not remember which of the
four affine point sheets carried it.

## 4. Six edges are not a descended tournament

Quotienting only by edge reversal gives

```text
Omega/<S> = {unordered edges of K_4},       |Omega/<S>|=6.  (16)
```

The ternary move does **not** descend to this quotient.  With a basis chosen
so `rho(1)=2`, the two representatives of one edge give

```text
(0,1) ~ (1,1),
R(0,1)=(0,2)       represents {0,2},
R(1,1)=(1,2)       represents {1,3}.                   (17)
```

The two images are different edges.  Forgetting orientation before applying
the ternary move is therefore illegal.

Nor is there an `A_4`-invariant choice of orientation for the six edges.  For
the edge `{x,x+d}`, left translation by `d` stabilizes that unordered edge and
swaps its endpoints.  Any invariant orientation would have to equal its own
reverse.  The six-edge carrier is consequently not an intrinsic tournament:
it needs an endpoint/origin sidecar before a binary orientation can be used.

## 5. The second binary role is genuinely different

Let `j in GL(2,2)` be a reflection of the three directions.  On flags put

```text
J(x,d)=(jx,jd).                                            (18)
```

Then

```text
J^2=1,                         JSJ=S,
JRJ=R^-1.                                                   (19)
```

Adding `J` extends the left affine group from `A_4` to

```text
AGL(2,2)=V_4 semidirect S_3=S_4.                            (20)
```

Thus the direction reflection is an additional odd involution.  It is not
the same operation as `S`, which has already been used as a double
transposition in the translation kernel.

This distinction is sharp.  Fix a three-cycle `R` on four points.  The ten
solutions of `S^2=1` in `S_4`--the identity and nine nontrivial
involutions--classified by their position relative to the fixed point of `R`,
give:

```text
involution type                         number   |<S,R>|
identity                                  1          3
transposition inside the R-orbit          3          6
transposition meeting the fixed point     3         24
double transposition                      3         12.     (21)
```

The last row is the tetrahedral flag quotient above.  A cross transposition
instead generates the full `S_4`, while an inner transposition leaves one
point fixed.  Merely naming an order-two generator does not determine which
information ledger is present.

## 6. Relation to the current frontiers

This theorem makes the following existing boundaries compatible without
erasing their scope.

1. The Farey tessellation underlying
   `THM-2056-kelvin-polar-farey-defect-certificate` admits the literal local
   `C_2/C_3` grammar.  Mapping those local moves to `(4)` still requires a
   chosen affine origin, a cyclic direction orientation, and transport of the
   polar owner inequality.  None is supplied here.

2. `THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go`
   sends the modular involution to a **linear direction reflection** and
   classifies its affine `S_3/S_4` lifts.  Here the modular involution is
   instead a **translation** in `V_4`, so the image is `A_4`.  These are
   different homomorphisms, separated exactly by the census `(21)`.

3. `THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin`
   classifies the six unordered quartic edges and their association-scheme
   and tournament boundaries.  The new datum here is the twelve-flag cover,
   its regular right action, and the explicit proof `(17)` that the ternary
   move does not survive the six-edge quotient.

4. `THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame`
   supplies the four-point `V_4` torsor and the three-direction `S_3`
   quotient.  After choosing its missing `C_3` orientation, `(1)--(20)` give
   the exact tetrahedral modular quotient.  Without that choice, odd
   reflections act semilinearly as in `(19)`.

5. `THM-3049-k4-matching-monomial-tropical-root-extraction-clutch` identifies
   the three disjoint perfect matchings of `K_4` as the positive
   Krenn--Gu--Soltesz `n=4,d=3` carrier.  They are exactly the three
   direction fibres in `(13)`.  THM-3067 lifts a matching direction to its
   four endpoint flags; it does not recover a matching-fibre amplitude after
   contraction.

6. `THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary`
   proves that no signed-pair-equivariant linear descent recovers the balanced
   phase sector from the fifth compound.  The finite flag bitorsor does not
   manufacture the missing cross-pair contraction.

7. `THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor`
   and `THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary`
   show that a fixed sheet plus a true cofactor contrast is load-bearing for
   Keller decoding.  A combinatorial point in `V_4` is not yet that physical
   sheet or cofactor.

The result therefore identifies the exact finite object on which the binary
and ternary grammars coexist, and the first two quotients that destroy it.  It
does not construct a carrier-preserving map from the Farey fan, a physical
LRC target action, a quartic Keller exclusion, an SFC closure, or any ledger
decrement.

## 7. Exact companion

Run

```text
python 04-computation/modular_tetrahedral_flag_bitorsor_thm3067.py
python -O 04-computation/modular_tetrahedral_flag_bitorsor_thm3067.py
```

The exact companion checks:

1. all twelve flags and the regular left/right `A_4` actions;
2. `S^2=R^3=(SR)^3=1`;
3. the three-element direction image and four-element origin kernel;
4. all six edge-reversal orbits and the explicit no-descent hostile `(17)`;
5. the absence of an invariant edge orientation;
6. the odd-reflection conjugacy laws and the full `S_4` extension; and
7. all ten solutions of `S^2=1` in the census `(21)`.

It uses explicit `require` gates and exact finite arithmetic.  Normal and
optimized executions reproduce the stored eleven-line transcript
byte-for-byte after LF normalization.
