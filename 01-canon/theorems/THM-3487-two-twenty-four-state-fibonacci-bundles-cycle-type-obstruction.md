---
id: THM-3487
title: "Two twenty-four-state Fibonacci bundles and the cycle-type obstruction"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  AUDIT PENDING.  Synchronizing the proved six-state Fibonacci channel-order
  cycle with the four-state P1(F3) cycle gives two 12-cycles, whereas its
  lawful affine-V4 owner bundle has four 6-cycles.  Their mismatch is the
  zero versus nonzero class in H^1(C6;V4): a nonzero seam class repairs the
  permutation type but forbids every closed owner section.  The three repair
  classes are the three K4 matching/XOR directions.
source: codex-2026-08-16-two-fibonacci-24-state-bundles
audit: >
  self-contained cycle-holonomy and graph-cohomology proof candidate; exact
  companion pins THM-3339 and the four-state frame-line probe, reconstructs
  both affine gauges, checks all cycle types, all three nonzero seam classes,
  conjugacy counts, matching-gauge fibres, normal/optimized replay, and AST
  security; independent proof audit pending
depends_on:
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
related:
  - THM-3486-critical-harmonic-transform-of-periodic-polynomial-words
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
script: 04-computation/fibonacci_two_24_state_bundle_obstruction_thm3487.py
output: 05-knowledge/results/fibonacci_two_24_state_bundle_obstruction_thm3487.out
script_sha256: 13ba2822f9f020623f59f361b8a2aac37da73eeb23c55b00157ab167adb677f9
output_sha256: 6c43abb71d127d4efe80c250f0651eeff1683c9a8df2781ddc81883f1b488a98
semantic_sha256: 9c338ec5eebe1f93326ae12c00f37da58382fa2084339adbd1de422e7d05d70e
hash_basis: LF-normalized bytes
---

# THM-3487 -- the missing transplant is a nonzero seam class

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
AUDIT PENDING.**

Equal cardinality does not identify two finite dynamical bundles.  In the
present `6*4=24` collision, the complete obstruction is sixth-step holonomy.
It also says exactly what extra datum repairs the abstract permutation and
exactly what that repair destroys.

## 1. Typed setup

Write

```text
V4=F_2^2={0,p,q,r},                    r=p+q.             (1)
```

THM-3339 proves that the consecutive Fibonacci/Farey flanks have the six
channel orders

```text
(q,r,p), (q,p,r), (p,q,r), (p,r,q), (r,p,q), (r,q,p),   (2)
```

and, on the same residue hexagon, the closed affine-owner section

```text
o=(0,p,p,r,q,q),
Delta o=(p,0,q,p,0,q).                                  (3)
```

The pinned four-state probe gives the Fibonacci action on `P1(F_3)` as the
4-cycle

```text
G=(0 1 2 3).                                             (4)
```

Synchronizing both packets by the same Fibonacci index defines the
frame-line state set and shift

```text
X=Z/6 x Z/4,
T_X(i,j)=(i+1,j+1).                                      (5)
```

This is a statement about the indexed Fibonacci path.  It is not a free
identification of the projective line with the affine owner set and not an
action of the full ternary Berggren branch monoid.

For the owner packet, let `L_i in GL_2(F_2)` be the unique linear map carrying
the three entries of order `i` in (2), position by position, to order `i+1`.
Put

```text
t_i=o_(i+1)+L_i o_i,
A_i(u)=L_i u+t_i,                                        (6)
Y=Z/6 x V4,
T_Y(i,u)=(i+1,A_i(u)).                                   (7)
```

Thus `A_i(o_i)=o_(i+1)` by construction.  The exact companion reconstructs
all six `L_i,t_i`; their ordered affine product is

```text
A_5 A_4 ... A_0=id_V4.                                  (8)
```

Equation (8) also follows conceptually: the ordered linear product is the
identity because the displayed frame closes, and the affine product fixes
`o_0` because the section (3) closes.

## 2. Cycle-holonomy lemma

Consider any skew product over an `m`-cycle,

```text
T(i,x)=(i+1,A_i x),                                     (9)
```

with bijective fibre maps.  Its return map over base point zero is

```text
H=A_(m-1)...A_0.                                        (10)
```

Each length-`s` orbit of `H` produces one length-`ms` orbit of `T`, and every
orbit of `T` arises this way.  Indeed a state can return to its base fibre
only after a multiple of `m` steps, and the induced motion there is a power
of `H`.  This elementary lemma is the whole cycle classifier needed below.

For (5), the sixth return is

```text
G^6=G^2=(0 2)(1 3).                                    (11)
```

It has two 2-cycles.  For (7), equation (8) gives four fixed points.  Hence

```text
cycles(T_X)=12^2,
cycles(T_Y)=6^4.                                        (12)
```

Equivalently,

```text
Fix(T_X^6)=empty,
Fix(T_Y^6)=Y.                                           (13)
```

Cycle type is a complete invariant of conjugacy for finite permutations, so
there is no bijection `Phi:X->Y` with

```text
Phi T_X=T_Y Phi.                                        (14)
```

This is stronger and better typed than saying only that the two fibres have
different interpretations: their one-step dynamical systems are not even
abstractly conjugate.

## 3. The obstruction is not a frame artefact

There is also a fixed-channel-frame owner transport

```text
B_i(u)=u+o_i+o_(i+1).                                   (15)
```

Its six translations telescope to zero.  Both the moving transport (6) and
the fixed transport (15) therefore have trivial return holonomy.  The fibre
gauges defined recursively by

```text
phi_0=id,
phi_(i+1)=B_i phi_i A_i^(-1)                            (16)
```

close at `i=6` and give the base-preserving conjugacy

```text
phi_(i+1) A_i=B_i phi_i.                                (17)
```

The companion constructs all six `phi_i` and checks (17) on all 24 states.
Thus the `12^2` versus `6^4` obstruction survives lawful changes between the
moving and fixed owner frames.  It is projective return holonomy versus
owner return holonomy, not a coordinate accident.

## 4. Exact H1 classification

After the linear trivialization (16), any affine owner transport around the
hexagon is represented by six edge translations

```text
c=(c_0,...,c_5) in V4^6.                               (18)
```

A vertex gauge `g=(g_0,...,g_5)` changes them by

```text
c_i -> c_i+g_i+g_(i+1).                                (19)
```

Consequently the only gauge-invariant datum is the seam sum

```text
h=sum_i c_i in V4,                                     (20)
```

and every value occurs.  This is the explicit cellular identification

```text
H^1(C_6;V4)=V4.                                        (21)
```

The lawful owner packet (3), (6) is class `h=0`.  Its closed section is the
corresponding flat section.  Changing the final seam by a nonzero

```text
h in {p,q,r}                                           (22)
```

makes the return map translation by `h`.  Every nonzero translation of V4
has cycle type `2^2`; the cycle-holonomy lemma therefore gives

```text
cycles(T_(Y,h))=12^2.                                  (23)
```

So each nonzero cohomology class repairs the abstract cycle type of `X`.
Conversely (23) forces `h!=0`; the classification is exact.

The tariff is equally exact.  A closed section would give a fixed point of
the return translation, but `u -> u+h` has no fixed point for `h!=0`.
Therefore:

```text
h=0       iff a closed affine-owner section can exist;
h!=0      iff the 12^2 projective cycle type is matched. (24)
```

No bundle in this family has both properties.  The repair does not merely
move the particular section (3); it forbids every closed section.

## 5. Why the three nonzero classes are K4/XOR matchings

On any four-point set, the identity together with the three double
transpositions is the canonical normal Klein subgroup of `S4`.  Its three
nonidentity elements are the three perfect matchings of `K4`.  Equation
(11) says that the projective return `G^2` selects one such matching.

On the affine owner set, the three nonidentity translations are

```text
u -> u+p,       u -> u+q,       u -> u+r.               (25)
```

They are the same three matching directions in XOR coordinates.  A point
identification `P1(F_3)->V4` transports `G^2` to exactly one map in (25).
Among the `4!=24` point identifications, exactly eight choose each direction.
Thus the projective packet determines a matching but the existing data do
not intrinsically identify that matching with `p`, `q`, or `r`: the missing
sidecar is precisely the channel-to-matching gauge already isolated in
THM-3339.

This explains the apparent three-way ambiguity without forcing a tournament.
It is an affine/projective torsor ambiguity.

## 6. Strongest survivor after forgetting the base map

Squaring the frame-line shift gives

```text
cycles(T_X^2)=6^4.                                     (26)
```

Hence `T_X^2` and `T_Y` are abstractly conjugate.  The number of conjugating
bijections is the centralizer size of a permutation of type `6^4`, namely

```text
6^4 4! = 31104.                                        (27)
```

This survivor does not give a bundle transplant over the displayed base.
The base component of `T_X^2` is `i->i+2`, while that of `T_Y` is
`i->i+1`.  A map preserving the base projection would have to equate these
two maps and cannot exist.

For each nonzero repaired class, the common type is `12^2`, so the number of
abstract conjugacies with `T_X` is

```text
12^2 2! = 288.                                         (28)
```

Again these bijections need not preserve the channel-order base, owner
section, branch letters, or any current.

## 7. Periodic harmonic address grades

There is one safe connection to subsets of the harmonic series.  Declare a
uniform 24-address enumeration, repeat it periodically, and mark a subset of
states that is invariant under one of the displayed permutations.  Such a
subset is a union of cycles.  Its reciprocal-address sum has logarithmic
coefficient equal to its cardinality divided by 24.  Therefore

```text
T_X or T_(Y,h!=0):  {0,1/2,1},
T_Y or T_X^2:        {0,1/4,1/2,3/4,1}.                 (29)
```

This is the degree-zero periodic mechanism developed more generally in the
THM-3486 candidate.  Equation (29) is **not** a time average along one orbit:
an invariant union has constant indicator on each orbit.  Nor does it say
anything about arbitrary subsets of the natural numbers.

## 8. Tournament and branch boundaries

The six base states in (2) are the six total orders of three matching
channels.  They may be drawn as transitive tournaments on three channel
labels, but they are not six vertices of a tournament.  The four projective
lines and four affine owners are different four-element torsors; neither is
an intrinsic tournament.  The double transposition (11) selects an unoriented
perfect matching, not six pairwise arc orientations.  Missing edges,
bidirectional edges, and XOR are therefore not silently collapsed into a
cosmetic tournament.

Likewise, (21) is an explicit finite path-cycle model for how a word/owner
transport can become an `H^1` class.  It suggests a useful discipline for the
harder LRC word-current versus Jacobian-flux map: identify both holonomies,
their coefficient local systems, and the lost gauge before comparing them.
It does **not** write that D5 map.

This candidate proves no full Berggren-tree conjugacy, no ternary branch
action, no physical current, no Jacobian flux, no LRC bispectrum
nonvanishing, and no case of LRC(14).

## 9. Exact companion

Run

```bash
python -B 04-computation/fibonacci_two_24_state_bundle_obstruction_thm3487.py
python -B -O 04-computation/fibonacci_two_24_state_bundle_obstruction_thm3487.py
```

The companion pins the source computations, reconstructs the unique linear
and affine edge maps, checks both owner gauges and their trivial holonomy,
computes (12)--(13), (26)--(28), checks all three nonzero seam holonomies,
and verifies the `8+8+8` matching-gauge census.  These finite checks support
but do not replace the proof above.  The candidate remains outside the
proved dependency graph until independent audit.
