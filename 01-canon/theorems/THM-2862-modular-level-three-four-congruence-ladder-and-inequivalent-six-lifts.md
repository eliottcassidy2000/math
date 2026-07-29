---
id: THM-2862
title: "Modular level-three/four congruence ladder and inequivalent six-point lifts"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  The level-three quotient
  SL2(F3)/{+-I} is A4 on four projective points, and those points are
  exactly its four C3 complements to the normal V4.  The first 2-adic
  thickening SL2(Z/4)/{+-I} is S4 on its four S3 complements, with
  reduction kernel V4 and marked triangle orders (2,3,4).  Its natural
  six-point projective line and the six edges of the complement torsor
  have the same three-block S3 quotient and the same V4 block kernel but
  are nonisomorphic: their point stabilizers are C4 and a nonnormal V4.
  The projective action is conjugation on the six four-cycles and retains
  quartic parity, whereas the edge action has identically even ambient
  sign.  Over one matching the common D8 stabilizer has three index-two
  subgroups, and the orientation, edge, and discriminant characters obey
  one exact Klein-four product law.  This is finite congruence and
  resolvent anatomy, not a graph-quartic realization or a Keller/LRC
  exclusion.
source: root/modular-level-three-four-congruence-2026-07-28
depends_on: []
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2746-c3-quotient-affine-v4-lifts-and-a4-projector-defect
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
script: 04-computation/modular_level_three_four_six_lifts_thm2862.py
output: 05-knowledge/results/modular_level_three_four_six_lifts_thm2862.out
script_sha256: ccd087583cdeaa230e1e0bf3ab36ee3e5de140e04834cda3952dbe9cb56fd0c0
output_sha256: 0a3a6bcd739a3c743aa15e1505fea00f1f5e273153dbeae9f7e879428da3a083
hash_basis: LF-normalized bytes
---

# THM-2862 -- levels three and four give two different six-point lifts

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

THM-2596 separates the actual Bass--Serre/Farey grammar from the slogan
that binary and ternary trees are literally the two finite factors of
`PSL2(Z)`.  THM-2595, THM-2743, and THM-2746 then identify two affine
quartic rows:

```text
A4 -> C3,                    S4 -> S3.                    (1)
```

The present theorem finds both rows inside honest modular congruence
quotients and then proves a sharp carrier warning.  At level four there
are two natural six-point lifts of the same three-point quotient.  One is
the quartic edge carrier of THM-2753.  The other is the projective line
over `Z/4`, equivalently the conjugacy class of four-cycles.  They retain
different information.

## 1. Marked modular quotients

Use

```text
S=[0 -1],                  C=[0 -1],
  [1  0]                     [1  1],

T=S^(-1)C=[1 1].
            [0 1]                                      (2)
```

To avoid ambiguity about `PSL2` over a ring, define

```text
G3=SL2(F3)/{+-I},             G4=SL2(Z/4)/{+-I}.         (3)
```

The images of `S,C` have orders two and three in both groups.  The
parabolic `T` has order three in `G3` and order four in `G4`.

### 1.1 Level three is the marked A4/C3 row

There are eight choices for the first nonzero column of a matrix in
`SL2(F3)`.  For each, the determinant-one second columns form an affine
line of size three.  Therefore

```text
|SL2(F3)|=24,                    |G3|=12.                (4)
```

The faithful action on the four points of `P^1(F3)` has census

```text
class size                 1       3       8
cycle type on P^1(F3)     1^4     2^2     3,1.           (5)
```

Hence `G3` is `A4`.  Its identity and three involutions form the normal
Klein four `K3`.  Each point stabilizer is a `C3`, every `C3` complement
to `K3` is a point stabilizer, and there are exactly four.  Thus

```text
P^1(F3)  ~=  {C3 complements to K3}                      (6)
```

equivariantly.  In the marking `(2)`, `S` is a double transposition, `C`
is a three-cycle with one fixed point, and `T` has order three.  After
quotienting by `K3`, the binary generator disappears and the ternary
generator remains.  This is exactly the finite `A4 -> C3` mechanism of
THM-2746, now realized by reduction modulo three.

### 1.2 Level four is the marked S4/S3 row

Reduction modulo two is onto `SL2(F2)`, because the reductions of `S,C`
generate it.  An element in its determinant-one kernel upstairs has the
form

```text
I+2A,                       tr(A)=0 in F2.                (7)
```

There are eight such matrices.  Quotienting by `+-I` leaves the four
classes

```text
I, I+2E12, I+2E21, I+2(E12+E21),                         (8)
```

so

```text
1 -> K4 ~= V4 -> G4 -> SL2(F2) ~= S3 -> 1,               (9)
|SL2(Z/4)|=48,                    |G4|=24.               (10)
```

The quotient acts faithfully on the three nonidentity elements of `K4`.
The extension splits: for example

```text
A=[1 1]
  [1 2]
```

has projective order three, and `<S,A>` maps isomorphically onto `S3`.
The standard `S3` action has no fixed nonzero vector in `K4`; equivalently,
the exact complement census gives four `K4`-conjugate `S3` complements.
Conjugation on those four complements is faithful.  Consequently

```text
G4 ~= V4 semidirect S3 ~= S4.                            (11)
```

On the complement torsor the marked generators have types

```text
S: 2,1,1,                   C: 3,1,                   T: 4. (12)
```

In particular the fixed pair of `S` misses the fixed point of `C`; this is
the nonzero affine compatibility class in THM-2595/2743.  The complete
four-point census is

```text
1*1^4,  6*(2,1,1),  3*2^2,  8*(3,1),  6*4.             (13)
```

Thus the two live marked quartic rows are not merely analogous to the
modular free factors:

```text
level 3:  (2,3,3) -> A4 -> C3,
level 4:  (2,3,4) -> S4 -> S3.                           (14)
```

The `4` in the second row is not another free factor.  It is the order
acquired by the parabolic at the first `2`-adic thickening.

## 2. Two six-point lifts of the same ternary quotient

Let `X` be the four `S3` complements in `(11)`, and let

```text
E=binom(X,2),                         |E|=6,
M={three perfect matchings of X},     |M|=3.             (15)
```

The map sending an edge to its opposite-edge pair gives the `S3` action
on `M` with kernel `K4`, as in THM-2753.

Now put

```text
P=P^1(Z/4),                           |P|=6.              (16)
```

A primitive vector modulo four is taken modulo multiplication by the two
units `+-1`.  Reduction modulo two gives

```text
P -> P^1(F2),                          6 -> 3,             (17)
```

with three fibres of size two.  The induced action on `(17)` is again the
`S3` quotient in `(9)`, and its kernel is again `K4`.

Both full six-point actions are faithful.  They are nevertheless not
isomorphic.  Their exact cycle tables, indexed by the common four-point
class in `(13)`, are

| four-point class | size | `P^1(Z/4)` | six edges `E` |
|---|---:|---|---|
| `1^4` | 1 | `1^6` | `1^6` |
| `2,1,1` | 6 | `2^3` | `2^2,1^2` |
| `2^2` | 3 | `2^2,1^2` | `2^2,1^2` |
| `3,1` | 8 | `3^2` | `3^2` |
| `4` | 6 | `4,1^2` | `4,2` |                      (18)

The smallest hostile is already the marked involution `S`:

```text
S on P: 2^3,                    S on E: 2^2,1^2.          (19)
```

Equivalently, a point stabilizer in `P` is cyclic of order four, with
element-order multiset

```text
1,2,4,4,                                                 (20)
```

whereas an edge stabilizer is a nonnormal Klein four, with multiset

```text
1,2,2,2.                                                 (21)
```

The normal `K4` in `(9)` must not be confused with the nonnormal edge
stabilizer in `(21)`.  It is the kernel of the common **three-block**
action, not the kernel of either faithful six-point action.

This proves the precise six-object no-go:

```text
same cardinality + same 3x2 blocks + same S3 block action
+ same V4 block kernel
does not determine the six-point carrier.                          (22)
```

## 3. The projective six-set is the four-cycle orientation lift

Let `p_infinity=[1:0]` in `P`.  Its stabilizer is `<T>`, because `T`
fixes it, `<T>` has order four, and transitivity in `(16)` makes every
stabilizer have order four.  Under `(11)--(12)`, `T` is a four-cycle on
`X`, whose centralizer in `S4` is exactly `<T>`.  Therefore

```text
g p_infinity |-> g T g^(-1)                                (23)
```

is a well-defined equivariant bijection from `P` to the six four-cycles
of `S4`.  Replacing `T` by `T^(-1)` gives the only other equivariant
bijection; fibrewise inversion exchanges the two projective lifts over
each mod-two point.  Thus `(23)` includes an explicit orientation marking
rather than pretending that the inverse choice is canonical.

The two inverse four-cycles in a cyclic subgroup have the same square.
Their square is one of the three nonidentity elements of the normal
`K4`, acting on `X` as a perfect matching.  Hence the common quotient in
`(17)` can be written

```text
four-cycle gamma |-> gamma^2 in K4\{1}.                   (24)
```

The corresponding quotient of `E` sends an edge to the nonzero `K4`
translation which exchanges it with its opposite edge.  Thus the two
six-sets are two different double covers of the same ternary matching
set:

```text
P: inverse orientations over a matching,
E: two constituent opposite edges over a matching.       (25)
```

This is the rigorous `3 times 2` co-occurrence object.  It is a block
cover, not a six-vertex tournament.  Any tournament orientation is extra
gauge data; the involutions in `(19)` already forbid treating it as an
intrinsic invariant relation.

## 4. Parity is retained by one lift and erased by the other

Read the ambient permutation signs from `(18)`.  On `P`, the class
`2,1,1` acts as three transpositions and the class `4` acts as one
four-cycle.  On all other classes both signs are positive.  Therefore

```text
sgn(Sym(P) action of g)=sgn(Sym(X) action of g)           (26)
```

for every `g in G4`.  The projective/four-cycle lift remembers the
quartic sign character.

On `E`, every row of `(18)` is even:

```text
sgn(Sym(E) action of g)=+1                                (27)
```

for every `g`.  Equation `(27)` is THM-2753's edge-parity erasure.
Equation `(26)` gives a new positive sidecar: parity can be recovered on
six objects, but only after selecting the correct six-object action.
Cardinality six and the common `S3` quotient do not choose it.

## 5. The D8 character triangle over one matching

Fix one matching `m in M`.  Its stabilizer

```text
D=Stab_G4(m)
```

has order eight and is dihedral.  It contains exactly three index-two
subgroups:

```text
K4                 normal Klein four, kernel of the matching action;
C4=<gamma>         stabilizer of one oriented four-cycle above m;
V4_edge            stabilizer of one constituent edge above m.  (28)
```

Let the corresponding quotient characters be

```text
chi_disc,             chi_or,             chi_edge:D -> {+-1}. (29)
```

Here `chi_disc` is the restriction of the four-point sign character,
because `D intersect A4=K4`.  The three nonzero characters of
`D_ab ~= V4` satisfy

```text
chi_disc=chi_or chi_edge.                                  (30)
```

Thus the missing bit in `(22)` is not decorative.  The orientation and
edge double covers are the two independent binary faces over the ternary
matching base; their product is the discriminant face.

There is an exact field-theoretic form.  Let `L/F` be an `S4`-Galois
closure and use the subgroup marking `(28)`.  Then `L^D/F` is the degree-
three matching-resolvent field, and

```text
L^K4/L^D,             L^C4/L^D,             L^V4_edge/L^D (31)
```

are its three quadratic subextensions inside `L`, corresponding to
discriminant, four-cycle orientation, and edge choice.  In characteristic
different from two, their square classes obey the product law `(30)`.
This does not say that a graph quartic supplies any preferred generator or
polynomial model for the two sextic fields in `(31)`.

## 6. What this changes -- and what it does not

The finite modular picture is now exact:

```text
C2*C3 --mod 3--> A4 -> C3,
C2*C3 --mod 4--> S4 -> S3.                                (32)
```

There is no coefficient-ring map which canonically sends the first
congruence quotient to the second.  After a marking, `A4` is of course the
even subgroup of `S4`, but that embedding is not induced by reduction
between levels three and four.

The useful consequence for the quartic Keller frontier is a stopping rule
and a concrete next object.

1. The abstract free-factor relations cannot exclude either live quartic
   row: both occur as honest congruence images.
2. The cubic matching resolvent retains only the common three-point base.
   A proof which needs parity, sheet orientation, or edge incidence must
   specify which quadratic lift in `(31)` it uses.
3. The cheapest new algebraic test is to construct the `C4`-stabilizer and
   nonnormal-`V4`-stabilizer sextic resolvents of a graph quartic and compare
   their divisor, different, and Jelonek sidecars over the same cubic
   resolvent.  No such comparison is claimed here.

Likewise no LRC carrier, owner current, or endpoint phase is produced.
THM-2596's Gram-owner sidecar remains essential before any modular/Farey
motion can preserve the Euclidean defect.  The theorem gives no graceful-
tree labelling, partial-cube embedding, Feuerbach incidence, or canonical
tournament on six vertices.

## 7. Exact companion and boundary ledger

Run

```text
python 04-computation/modular_level_three_four_six_lifts_thm2862.py
python -O 04-computation/modular_level_three_four_six_lifts_thm2862.py
```

Both executions byte-match

```text
05-knowledge/results/modular_level_three_four_six_lifts_thm2862.out.
```

The standard-library companion uses exact modular arithmetic and explicit
exceptions.  It enumerates both special linear groups and their projective
quotients, proves the generator orders and group closures, finds every
complement, checks every action law and cycle census, constructs `(23)`,
identifies both block quotients, and verifies every stabilizer, sign, `D8`
subgroup, and character-product claim.  There is no Python `assert`,
floating-point decision, optional CAS, or scratch dependency.

```text
PROVED IN THE CANDIDATE: level-3 A4/complement-torsor realization;
                         level-4 V4 extension and S4 complement torsor;
                         marked (2,3,3)/(2,3,4) rows;
                         two exact six-point actions and their common quotient;
                         C4 versus V4 stabilizer nonidentification;
                         projective=four-cycle conjugacy action;
                         parity retention versus edge erasure;
                         D8 three-character and fixed-field triangle.

NOT PROVED:              a graph-quartic or Keller realization;
                         exclusion of A4 or S4 monodromy;
                         a canonical map between levels three and four;
                         sextic divisor/different/Jelonek control;
                         a physical LRC carrier or endpoint current;
                         JC(2), DC(2), G1, or LRC(14).                  (33)
```

**Candidate proof complete; independent hostile audit pending.**
