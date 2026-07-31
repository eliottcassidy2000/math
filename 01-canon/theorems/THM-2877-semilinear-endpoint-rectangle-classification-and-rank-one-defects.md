---
id: THM-2877
title: "Semilinear endpoint-rectangle classification and rank-one defects"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  No affine, prime-field semilinear, or coordinatewise projective
  permutation gives a full q0/q3/q11 endpoint-rectangle bridge.  Exact
  nearest maps leave one rank-one factor defect: nine points for q0 to
  q3, ten lost plus ten gained for the q3 edge, and nine lost plus nine
  gained for the q11 edge.  The positive q0-to-q3 injection preserves
  the character-three line but moves the named origins and carries no
  physical interval or ancestry action.  No row exclusion or LRC(14)
  proof follows.
source: root/lrc-semilinear-rectangle-bridge-2026-07-28
depends_on:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2859-horn-collar-q0-hinge-minimal-v4-globalization-and-witt-endpoint-obstruction
related:
  - THM-2863-endpoint-prony-splitter-and-carry-character-three-intertwiner
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression
  - THM-2876-q3-affine-defect-source-phase-annihilation-and-localized-augmentation-shadow
script: 04-computation/lrc14_semilinear_endpoint_rectangle_classification_thm2877.py
output: 05-knowledge/results/lrc14_semilinear_endpoint_rectangle_classification_thm2877.out
script_sha256: 10a4c965d02f9fab60f135d0bf10184096eeb70d30d0d9ef7ad4fcf5fc1aa447
output_sha256: 9ab7526b1e268176a98f0cf21238bc4cbabad0ce30c28b6565d0ca7b9d9a751e
hash_basis: LF-normalized bytes
---

# THM-2877 -- Semilinear endpoint-rectangle classification and rank-one defects

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem is not a row exclusion or an LRC(14) conclusion.  It pins and
reconstructs the endpoint masks from THM-2859/2847, exhausts the natural
affine groups, and is interpreted below against the promoted
THM-2868/2874 coefficient atlas and Kummer--Galois clutch.

## Verdict

There is no affine, monomial-substitution, or prime-field semilinear
permutation of `F_13^2` carrying the full q0 endpoint rectangle to a full q3
or q11 endpoint rectangle.  There is also no such map carrying either the
q3 or q11 step-2 mask to its step-68 mask.

The obstructions occur at different levels.

1. q3 has mass `90`, while q0 has mass `81`.  No bijection of the endpoint
   plane can identify them.
2. q11 has mass `81`, but every affine equivalence between the relevant
   proper rectangles is forced to be monomial by Cauchy--Davenport.  The
   required one-dimensional affine factor equivalences do not exist.
3. `F_13` has no nontrivial Frobenius automorphism.  Thus
   `A Gamma L(2,13)=AGL(2,13)` here.
4. Diagonal coefficient phases do not change Boolean support.  Even
   allowing a separate `PGL_2(F_13)` Möbius transformation on each
   coordinate gives no full q0/q11, q3-edge, or q11-edge rectangle bridge
   in either the direct or swapped orientation.  Individual non-affine
   factor automorphisms and equivalences do exist, so the full-rectangle
   qualifier is essential.

The useful positive is a low-rank failure anatomy rather than an action.
The q0 start mask embeds affinely and character-compatibly in the q3 start
mask with one extra nine-point horizontal fibre.  The closest q3 edge
differs by one ten-point vertical fibre after a vertical dilation; the
closest axis-preserving q11 edge differs by one nine-point vertical fibre.
Those two dilations do not preserve the labelled character-three line, and
the q11 map also moves the named origins.

## 1. Exact rectangle bank

The script reconstructs the full source endpoint masks at the physical
step-2 and step-68 intervals.  Source and target endpoint frames agree at
these states by THM-2859.  In point-support convention:

```text
q0, step 2:
  A0={0,1,2,3,4,5,6,7,12}
  B0={0,1,2,3,4,5,8,9,10}

q0, step 68:
  A0
  B0'={0,1,2,5,6,7,8,9,10}

q3, step 2:
  A3={0,1,2,3,4,5,6,7,8,9}
  B3={0,3,4,5,8,9,10,11,12}

q3, step 68:
  A3
  B3'={0,5,6,7,8,9,10,11,12}

q11, step 2:
  A11={0,1,2,3,4,5,8,9,12}
  B11={0,1,2,3,4,5,8,11,12}

q11, step 68:
  A11
  B11'={0,1,2,5,6,7,8,11,12}.
```

The q0 physical label `Z^8` acts on mask functions by pullback.  Therefore
its action on support points is `b -> b-8=b+5`, explaining the apparent
sign difference between THM-2859's operator convention and the displayed
`B0 -> B0'`.

## 2. Why full affine mixing cannot help

Let `S=A x B`, with `|A|=|B|=9`, and suppose

```text
g(x)=Mx+t
```

maps it to another `9 x 9` coordinate rectangle.  A row `(r,s)` of `M`
has projection `rA+sB`.  If `r,s` are both nonzero, Cauchy--Davenport gives

```text
|rA+sB| >= min(13,9+9-1)=13.
```

The target projection has size nine, so each matrix row has at most one
nonzero entry.  Invertibility makes `M` monomial.  The audit also verifies
directly that every mixed q0 projection is all of `F_13`.

The same argument applies to the q3 `10 x 9` edge: a mixed row again has
full projection.  Hence the exhaustive `AGL_2` search is not concealing a
non-axis-aligned equivalence.

For q0-to-q11, the one-dimensional search already stops both monomial
orientations.  A particularly transparent invariant is

```text
F_13\A0={8,9,10,11},
```

a four-term affine arithmetic progression, whereas the complements of
`A11` and `B11` are `{6,7,10,11}` and `{6,7,9,10}`.  Neither is an affine
four-term progression.  The complete factor table checks all step-2 and
step-68 direct and swapped possibilities; none completes both axes.

The computation exhausts

```text
|GL_2(F_13)|=26,208,
|AGL_2(F_13)|=4,429,152,
```

and finds no q0-to-q11 map, no q3 edge, no q11 edge, and no single affine
map carrying the whole q0 step-2/step-68 pair to the q3 or q11 pair.

## 3. The q0 germ is not uniquely affine

The q0 edge itself has exactly eight affine realizations.  They are the
product of the two automorphisms

```text
a -> a,          a -> -a+6
```

of `A0` with the four maps

```text
b -> b+5,
b -> 5b+8,
b -> 8b+7,
b -> -b+10
```

from `B0` to `B0'`.  The first vertical map is the support-point form of
the physical `Z^8` pullback.

This ambiguity disappears only after retaining the carry character.  For
an affine state map `y=Mx+t`, pushforward sends a Fourier covector `k` to
`kM^{-1}`; translation changes only its scalar normalization.  Requiring
the labelled vertical line `(0,3)` to remain `(0,3)` leaves exactly the two
maps with vertical rule `b -> b+5`.

Neither of these maps preserves the two THM-2863 named origins

```text
O={(0,0),(12,0)}
```

even setwise.  Their images are respectively

```text
{(0,5),(12,5)},    {(6,5),(7,5)}.
```

In fact no one of the eight exact q0 edge maps fixes `O`.  Thus the facts
that both named q0 origin occupancies remain nonzero before and after the
physical step do not define an equivariant two-origin column.  This is a
support-level boundary complementary to promoted THM-2868, not an
obstruction to that theorem: THM-2868 uses a signed origin difference,
26 actual multiplier samples, local Prony transport, and a projective
coefficient ratio; it never asserts that one positive q0 full-mask map
fixes the named origins.  Support transport, character transport, and
named-column coefficient transport remain three distinct requirements.

## 4. Nearest maps and the surviving low-rank defects

### 4.1 q0 to the q3 start mask

The affine translation

```text
g(a,b)=(a+1,b+8)
```

is an exact injection of the q0 step-2 mask into the q3 step-2 mask:

```text
E_(q3,2) = g(E_(q0,2)) disjoint-union ({9} x B3).        (1)
```

The defect is exactly nine points, and this is globally optimal over all
`4,429,152` affine maps.  The vertical slope is one, so `(1)` preserves the
labelled character-three line up to a scalar phase.  It moves the named
origins, however, and it is only an endpoint-address injection: it does not
map the physical interval, the empty q-shifted `U` ancestry, or the current
provenance.

Equation `(1)` is the smallest natural support enlargement found in this
audit.  It is more precise than the bare mass mismatch: the q3 excess is one
horizontal fibre, not nine unrelated points.

### 4.2 The q3 edge

The globally closest affine map is

```text
g3(a,b)=(a,2b).
```

It preserves the horizontal factor and fixes the named origin coordinates,
but

```text
E_(q3,68)
 = (g3(E_(q3,2)) \ (A3 x {3})) union (A3 x {12}).        (2)
```

Thus the minimum Hamming defect is `20`, one lost and one gained ten-point
vertical fibre.  The map sends the labelled vertical character `3` to
character `8`, so the best character-three-preserving affine map has the
larger Hamming defect `40`.

A coefficient automorphism sending `omega` to `omega^2` would formally send
that character `8` back to `3`.  This is **not** THM-2857's relative carry
Galois action: the relative group fixes `omega`, while such a base-field
automorphism also moves the coefficient normalization.  More sharply, among
the `78` automorphisms of `Q(zeta_2366)` inducing
`omega -> omega^2`, none fixes THM-2863's order-182 left Prony node; a
13-element coset fixes the order-14 right node.  Because the two node orders
differ, no such automorphism preserves the unordered node pair.  Testing a
**transformed** four-sample Prony bank is therefore a legitimate next
experiment, but the existing normalized bank cannot be imported unchanged.

### 4.3 The q11 edge

An axis-preserving closest map can be taken as

```text
g11(a,b)=(a,6b+1).
```

It has the exact defect

```text
E_(q11,68)
 = (g11(E_(q11,2)) \ (A11 x {10}))
     union (A11 x {11}).                                (3)
```

The minimum Hamming defect is `18`.  This map moves the named origins and
sends character `3` to character `7`.  Imposing either exact
character-three preservation or pointwise origin preservation raises the
minimum Hamming defect to `36`; imposing both still gives `36`.
Correcting the label by a coefficient automorphism
`omega -> omega^6` has the same Prony boundary: all `78` possible lifts move
the left node, while exactly `13` fix the right node.

For comparison, the constrained direct-edge minima are:

```text
             unrestricted   chi3-preserving   origins fixed   both
q0                 0               0                36          36
q3                20              40                20          40
q11               18              36                36          36.
```

No entry except the first two q0 entries is zero.

## 5. Horizontal, carry, origin, and E3 typing ledger

```text
map / object:
  affine endpoint-address permutation or partial injection;

preserved by the strongest positive (1):
  Boolean support inclusion, coordinate axes, q3 horizontal factor after
  one-fibre completion, and the labelled character-three line up to phase;

destroyed or moved:
  the named THM-2863 origins, absolute endpoint phase, physical interval,
  owner/current provenance, and q-shifted U ancestry;

E3:
  q0, q3, and q11 all lie in the E3 block on the 20-cell horn, so the
  present support comparisons do not cross macro truth.  Nothing here
  transports q11 to the q7 complement block or contracts THM-2847's
  rank-one E3 mapping cone;

projective boundary:
  the full coordinatewise Möbius table has non-affine individual factor
  maps (for example four non-affine maps among the six maps
  A11 -> A11 and among the six maps A11 <-> B11), but none supplies the
  companion factor needed for a full q0/q11 or direct q3/q11 edge
  rectangle bridge.  Promoted THM-2868 proves exactly the surviving
  projective route: it acts on signed coefficient ratios/projectors, not
  as a positive coordinatewise support bijection.  THM-2874 then gives
  its explicit Q(zeta_91)-rational clutch to the centered endpoint Galois
  orbit and locates the Bockstein seam obstruction.
```

## 6. Cheapest next decisive experiment

THM-2876 executes the q3 one-fibre signed experiment.  Its endpoint
numerators retain the complete Hermitian atlas, but the inverse source
phase annihilates the physical boundary termwise.  That route is closed.

The remaining positive control is the q0-to-q3 injection `(1)`.  Its
missing `{9} x B3` fibre has rank one as a rectangle factor and its
vertical translation already preserves the character-three line.  The
next decisive test is to adjoin precisely this fibre, then replay the named
origins, common source coefficient, 26 Prony samples, macro-E3 truth,
physical interval, and `QA/QAB` ancestry.  Any failure must identify the
first lost coordinate; any success must produce a nonzero full current,
not another endpoint-only coefficient shadow.

## Reproduction

```bash
python3 04-computation/lrc14_semilinear_endpoint_rectangle_classification_thm2877.py
python3 -O 04-computation/lrc14_semilinear_endpoint_rectangle_classification_thm2877.py
```

The companion contains no executable Python `assert` statements.
Normal and optimized outputs are byte-identical.  SHA-256:

```text
script  10a4c965d02f9fab60f135d0bf10184096eeb70d30d0d9ef7ad4fcf5fc1aa447
output  9ab7526b1e268176a98f0cf21238bc4cbabad0ce30c28b6565d0ca7b9d9a751e
```

The independent audit rebuilt all six masks and source/target frames,
repeated one-pass exhaustion of all `4,429,152` affine maps, and reproduced
every exact count, overlap minimum, constrained minimum, q0 map census,
one-row defect, covector law, named-origin image, and cyclotomic rechart.
It separately enumerated coordinatewise `PGL_2(F_13)`: non-affine
individual factor maps exist, but none supplies the companion factor for a
full rectangle bridge.  Normal and optimized audit modes byte-match.
