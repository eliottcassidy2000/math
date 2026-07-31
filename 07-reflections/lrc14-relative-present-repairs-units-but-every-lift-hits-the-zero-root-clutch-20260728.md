> **CORRECTION — MISTAKE-310 (2026-07-28).**  Section 4 correctly computes
> the edge-preserving root label `0`, but its conclusion that the physical
> private-root lift relation is empty is false.  The right-root-`0` and
> left-root-`1` half-tooth charts overlap on the open interval `(1,13)/182`;
> an exact audit finds a restricted coefficient-preserving rechart there.
> The census, unit tables, scalar-lift counts, and full-target coexistence
> below survive.  The no-go survives only for the unrestricted
> edge-preserving labelled row.  THM-2744 is reserved for the corrected
> partial-clutch statement.

# Relative present repairs the unit table, but every scalar lift hits the zero-root clutch

**Status: FINITE-EXACT + VERIFIED in byte-matching normal and optimized
modes.**  This is a canonical-fibre support and obstruction result.  It
does not produce a private-root physical edge, endpoint amplitude, row
exclusion, or proof of LRC(14).

## 1. The three nearby objects are not interchangeable

The new full two-target calculation
`lrc14-the-present-wall-is-a-one-target-projection-and-the-full-target-sheet-repairs-it-20260728.md`
proves that the lawful sheet `U_(s,t)` repairs the scalar wall: for every
`ell!=0`, every `s`, and every `t!=0`, the canonical row has positive

```text
E3 intersection F_(ell,s,t) intersection D^(-6)Q_(3,{1,2}),
```

and every nonzero `t`-character survives.  That result has no rail, carry,
private root, primitive unit, or physical lift.

THM-2712 instead gives 304 following-atom endpoints and a physical inner
support SCC, but no semantic THM-2305 source.  The present computation first
checks the intersection directly:

```text
full THM-2707 packet bank (3346):  (source-free,E2)=(3081,265),
THM-2712 atom bank (304):          (source-free,E2)=(280,24).
```

Every `E2` point has empty terminal word.  The complete semantic record is
constant on the whole inherited open cylinder `I`; the minimum stability
radius is `337` times its transported radius.  Thus neither existing packet
bank contains even one semantic handoff vertex.

## 2. The exact relative-present coefficient

Broaden only from THM-2707's fixed rail to the fourteen source-one rails,
keeping its exact address orbit, delayed sector zero, `h=6`, `kappa=1`, and
clocks `1,...,6`.  Before present, the semantic census contains exactly

```text
12848 E3 -> Q_(3,{1,2}) fork points,
address residue 7: 6404,
address residue 8: 6444.
```

Every point lies on exactly one source-one rail, has its forced strict private
half-tooth, and lies outside **all thirteen** one-target present labels.

It would be ill-typed to attach the old primitive-unit flag after deleting
the factor which defined it.  Recompute the seven-coordinate vector three
ways, coefficient by coefficient:

```text
V_old  = coefficient on rail intersection Present,
V_free = coefficient on the rail with Present omitted,
V_rel  = V_free - V_old
       = coefficient on rail intersection Present^c.          (1)
```

Exact delayed-carry linearity proves `(1)` before taking determinants.  All
three restricted banks have gcd `2122120=26*81620`, and division by the
canonical content `26` gives the unit rails

```text
                 residue 7                  residue 8
V_old            all 14                     none
V_free           all except rail 13          all except rail 3
V_rel             all 14                     all except rail 3. (2)
```

Thus the relative complement repairs the missing residue-8 determinant on
thirteen of fourteen rails.  It leaves

```text
6404 residue-7 endpoints,
6178 residue-8 endpoints,
79127824 directed scalar cross-lifts.                         (3)
```

All `12582` endpoints retain the whole transported open cylinder `I`; the
binding radius equals the open-cylinder radius, so only excluded endpoints
touch a wall.  There are `1083056` reverse-clock unordered pairs.

## 3. The full target sheet survives on these endpoints

At each endpoint, retain the actual THM-2365/2407 labels

```text
s: q1 safe at -s/13 and c2 safe at +s/13,
t: q2 safe at -t/13 and c3 safe at +t/13.                 (4)
```

Every relative endpoint has between nine and eleven allowable **nonzero**
`t` values; none admits `t=0`.  Hence the full `U_(s,t)` repair and the
relative-complement coefficient really coexist pointwise.  For example,

```text
n=618: rail 2, clock (1,1), edge R, root 12,
n=619: rail 2, clock (1,1), edge L, root 1,
physical numerators (7,13^6-7),
common s=(0,2,3,4,7,8,9,10,11,12),
common t=(2,3,6,7,8,9,10,11).                            (5)
```

So the missing target coordinate is no longer the obstruction.

## 4. Uniform zero-root clutching no-go

The apparent bipartite carrier in `(3)` is not a private-root lift graph.
The two residue classes force opposite half-tooth edges:

```text
residue 7: carry 12, right edge, root 12,
residue 8: carry  6, left  edge, root  1.                 (6)
```

A physical numerator from residue 7 to residue 8 has `2k=1 mod13`.
Physical translation preserves the half-tooth edge, so on the right edge it
sends root `12` to root `0`.  The reverse numerator has `2k=-1` and sends
left root `1` to left root `0`.  The lawful target in each direction exists
only after switching to the opposite edge.  Consequently **every** scalar
cross-lift in `(3)` factors through the forbidden zero-root seam, and the
induced private-root lift relation is empty.

This is the first failed implication:

```text
semantic E3 fork + nonzero t + relative unit + scalar physical lift
  -/-> private-root physical lift.                         (7)
```

The missing coordinate is now an edge-clutch transition at root zero, not a
present label, target colour, rail, clock, or determinant.

## 5. Why this is not THM-2716's `C4` arm

There is a visual `C2` resemblance, but no lawful identification.  Here the
bit is the **local left/right half-tooth chart** and the needed switch occurs
at one unchanged circle point after an edge-preserving translation reaches
root zero.  In THM-2716 the objects `E,G` are **endpoint/ghost macro types**
of a fixed abstract `C4/H` transporter, and its odd arms have relative degree
one between those types.

No map from `(edge,root)` to `E,G` has been constructed; it would not preserve
the defining observable (local comb membership versus macro endpoint/ghost
type).  Moreover THM-2716 proves that the odd-order THM-2707 translation
group has no nontrivial `C4` or deck involution.  The required recharting is
therefore not the `C4` quarter-turn arm.  At most, the root-zero overlap is a
new candidate **source** for a clutch functor into that transporter.  Such a
functor would need an explicit predicate-preserving map and amplitude; none
is inferred here.

## 6. Next decisive test and scope

The highest-leverage remaining object is a root-zero overlap map between the
two private half-tooth charts, retaining `(ell,s,t)` and the relative vector
before taking its determinant.  The cheapest hostile asks whether this chart
change intertwines the two delayed-carry vectors.  If it does not, `(7)` is a
terminal no-go for the source-one fibre; if it does, it supplies exactly the
missing physical clutch sidecar for the already positive full target sheet.

Reproduce with

```bash
python3 04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py
python3 -O 04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py
```

against

```text
05-knowledge/results/lrc14_relative_present_semantic_lift_probe_20260728.out
```

The normal transcript and the optimized transcript byte-match the stored
output.  SHA-256:

```text
script  f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f
output  861e920b945aafde3bc24c6089ba630763035f0919d738b77a0525b91ba6fa74
```

No target-action intertwiner, private-root lift, spectral endpoint pair, row
exclusion, or LRC(14) conclusion follows.
