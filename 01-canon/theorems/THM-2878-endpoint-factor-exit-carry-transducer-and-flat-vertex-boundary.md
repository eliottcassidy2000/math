---
id: THM-2878
title: "Endpoint-factor exit carry transducer and flat vertex boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Literal guard/q5 truth is the chain SS<DS<DD over every lawful endpoint
  representative, so the proposed fourth corner does not exist.  A
  faithful guard/u1/q5 rechart gives a marked square at the distinguished
  zero address, but every vertex encoding remains flat.  The directed u3
  D-to-S exit at q12-to-q0 is unique among the eighteen oriented factor
  toggles at zero address and, more sharply, is the unique carry marker
  present uniformly at all 169 addresses.  Its positive-path count equals
  floor((q+h)/13) on all 169 edges.  It derives the nonsplit carry law on
  all 2,197 compositions and the seam phase omega^3.  The event supplies
  the carry transition, not an initial ancestry state or common physical
  QA-to-QAB current; no row exclusion or LRC(14) proof follows.
source: root/lrc-missing-corner-constructor-2026-07-28
depends_on:
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
related:
  - THM-2859-horn-collar-q0-hinge-minimal-v4-globalization-and-witt-endpoint-obstruction
  - THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression
script: 04-computation/lrc14_endpoint_factor_exit_carry_transducer_thm2878.py
output: 05-knowledge/results/lrc14_endpoint_factor_exit_carry_transducer_thm2878.out
script_sha256: 8e75f9c9caa9482db17b40d220f5e8ef68d00d54d4ed481869f9b31ec146f1ca
output_sha256: fadccf2c9afa8cb386066b400ea51a0b90210ec26c0cd70d77b8354e279c0c38
hash_basis: LF-normalized bytes
---

# THM-2878 -- Endpoint-factor exit carry transducer and flat vertex boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem proves no row exclusion and no LRC(14) conclusion.

**Correction lineage.**  The first promoted wording confused a
distinguished zero-address census with a global uniqueness statement.
MISTAKE-320 records the repair.  The global survivor is stronger in the
coordinate that matters: `u3:D->S` is not pointwise unique at every address,
but it is the unique marker that works at **every** address.

## Verdict

There are two sharply different results.

First, the literal missing guard/q5 corner is globally impossible on the
selected source/target atom orbit.  The lawful endpoint characters are

```text
ell(a,b)+sW,                 s in F13,
W mod13=(1,1,1,1,1,1,0,0,0),
```

so every representative has `ell0=ell5=s`.  Exhaustion of both atoms, all
`169` addresses, all `13` representative gauges, and all `13` residues
gives

```text
(guard,q5) in {SS,DS,DD},       never SD.
```

Second, the complete nine-bit word contains a more useful object than the
missing vertex: one oriented factor-transition event realizes exactly the
THM-2851 carry cocycle.  Factor `u3` (index `3`, speed `40`) is dangerous
only at residue `12`.  With the positive unit orientation fixed, define

```text
kappa(q,h)
 = number of u3 transitions D->S
   on q,q+1,...,q+h,             0<=h<=12.
```

Then all `169` rows satisfy

```text
kappa(q,h)=floor((q+h)/13).                            (1)
```

Thus the endpoint factor path groupoid already realizes the carry
transition function.  It does not supply the initial ancestry coordinate;
the smallest faithful extra fibre remains `C13`, producing the joint
`C169` object.

## 1. Exact representative and corner no-go

The relation packet has rank `6`, its annihilator has `13^3=2197`
elements, and the old character quotient has `169` classes.  The `2197`
values

```text
ell(a,b)+sW
```

are pairwise distinct and exhaust the annihilator.  Canonical
representatives have `(ell0,ell5)=(0,0)`, while the full lawful
representative orbit projects to

```text
(ell0,ell5)=(s,s)
```

with multiplicity `169` for each `s`.

On either the source atom or its target translate, the zero-chart literal
word orbit is

```text
q0  SSSSSSSSD       q7  DSSSSDSSD
q1  SSDSSSSSD       q8  DSSSSSSSD
q2  SSDSSSSSD       q9  SSSSDSSSD
q3  SSSSSSSSD       q10 SSSSDSSSD
q4  SDSSSSSSD       q11 SSSSSSSSD
q5  DDSSSSSSD       q12 SSSDSSSSD
q6  DSSSSDSSD.
```

In particular

```text
guard danger = {5,6,7,8},
q5 danger    = {6,7},
```

so q5 danger is strictly nested inside guard danger.  The complete
`57,122`-row census is

```text
SS 39,546;  DS 8,788;  DD 8,788;  SD 0.
```

The same census holds before and after lawful marked representative
transport

```text
(ell,atom_q) -> (ell+sW,atom_(q-s)).
```

All nine factor bits, not only guard/q5, are preserved.  Source/target
swap, address change, residue reversal, or affine residue offset only
permutes this already exhausted universe.

## 2. Zero-address 36-pair perspectives and the all-address audit

At the distinguished zero address, for factor order

```text
(guard,u1,u2,u3,u4,u5,target_a,target_b,deep_C3),
```

the danger arcs are

```text
({5,6,7,8},{4,5},{1,2},{12},{9,10},{6,7},empty,empty,F13).
```

Among all `36` unordered pairs, the image-size histogram is

```text
one state: 3;  two states: 18;  three states: 14;  four states: 1.
```

the unique honest full square is `(guard,u1)`.  Its two danger arcs cross at
`q5`; neither contains the other.  The only three-state pair missing `SD`
is `(guard,u5)`, the original guard/q5 horn.  Thirteen pairs have disjoint
danger arcs and miss `DD`; the remaining pairs have a constant coordinate.

This uniqueness is deliberately local to that section.  Across all `169`
canonical addresses, the number of full pairs per address has histogram

```text
0:54, 1:64, 2:38, 3:10, 4:3,
```

with `38` distinct full-pair sets.  Each of

```text
(guard,u1), (guard,u2), (u1,u2), (u1,u4),
(u1,u5), (u2,u4), (u2,u5)
```

occurs at exactly `26` addresses.  The minimal hostile witness is address
`(1,0)`, whose unique full pair is `(u1,u5)`, not `(guard,u1)`.

At zero address, clean guard/u1 squares with the other seven bits fixed are

```text
(qSS,qSD,qDD,qDS)=(0,4,5,8),(3,4,5,8),(11,4,5,8).
```

Each directed boundary winds once around `C13`.  These are endpoint-factor
squares, not physical macro-current squares.

## 3. Faithful zero-address q7 parity rechart

Let

```text
j = q5 XOR u1.                                         (2)
```

Unlike an arbitrary binary relabeling, this loses no endpoint truth:
retaining `u1` recovers

```text
q5 = j XOR u1.                                         (3)
```

The pair `(guard,j)` has all four values.  At canonical q7,

```text
(guard,q5,u1,j)=(D,D,S,D),
```

so `j` agrees with the original q5 bit.  The q-labelled square

```text
(qSS,qSD,qDD,qDS)=(3,4,7,8)
```

retains q7 and has positive edge lengths

```text
(1,3,1,8),          sum=13.
```

At fixed physical q7, the same projected square appears on the
representative-fibre pullback at

```text
(sSS,sSD,sDD,sDS)=(9,10,0,1),
```

with words

```text
SSSSSSSSD, SDSSSSSSD, DSSSSDSSD, DSSSSSSSD.
```

Only the `s=9` vertex lies in endpoint `PAT_E3`; the other three require
vertex-dependent complement polarization.  All four restrictions are the
same physical atom and have identical endpoint values

```text
(231164267889491750,630230755085920022)
```

in the two certified fields.  Hence this coefficient boundary is flat.
Moreover, lawful co-transport sends each fixed-q7 chart to `q7-s` and
collapses all four words to canonical `DSSSSDSSD`.  The square is therefore
a useful marked representative-fibre locator, not descent on the old
`Ghat`, not a physical transporter, and not a `V4` action.

## 4. The directed full-word event is exactly carry

The guard/q5 filtration

```text
SS=0 < DS=1 < DD=2
```

has path sequences

```text
q3 -> q7 direct:       0,0,1,2,2
q3 -> q11 -> q7:       0,0,1,2,2,1,0,0,0,
                       0,0,0,0,0,0,1,2,2.
```

Both have net threshold rise `2`.  The via route has extra cancelling
down/up crossings, so every vertex filtration or signed coboundary is
blind to the carry.

The full word is finer.  Its address-uniform positive basepoint transition
is

```text
q12 -> q0:
SSSDSSSSD -> SSSSSSSSD,
```

the `u3` one-bit `D->S` event.  At zero address, exhausting all nine bits
and both state orientations shows that `(u3,D->S)` is the unique one of
`18` candidates whose positive-path count equals `(1)` on all `169`
`(q,h)` edges.

The all-address quantifier is subtler and sharper.  The same `u3:D->S`
marker works at every one of the `169` canonical addresses.  At `121`
addresses it is the only candidate; at `44` addresses one shifted `u1` or
`u2` event also works; and at `4` addresses two shifted events also work.
Thus `48` addresses have extras.  For example, address `(7,0)` admits both
`u1:D->S` and `u3:D->S`.  Each orientation of `u1` and `u2` occurs at only
`13` addresses, while `u3:D->S` occurs at all `169`; hence `u3:D->S` is
the unique **address-uniform** carry marker.

The event count is additive under concatenation.  All `2197` triples
`(q,h,k)` satisfy

```text
kappa(q,h)+kappa(q+h mod13,k)
 =kappa(q,h+k mod13)+floor((h+k)/13).                  (4)
```

Consequently `(4)` derives, rather than merely resembles,

```text
L_k L_h
 =T^floor((h+k)/13) L_(h+k mod13).                     (5)
```

On the seam,

```text
kappa(3,8)=0,        kappa(11,9)=1,
kappa(3,4)=0.
```

Thus the q3-via-q11 path differs from the direct q3-to-q7 path by one
central carry, and character three evaluates it as `omega^3`.

This is marked-gauge covariant.  Under `ell->ell+sW`, the singleton marker
moves to effective residue `q+s=12`; simultaneous physical transport
`q->q-s` restores `(1)`.  The companion checks `2197` such gauge rows, and
all `169` canonical address representatives have `ell3=0`.

The positive orientation is load-bearing.  A reversed path also sees one
`D->S` exit from the singleton, but through the other neighbour.  Truth
bits alone cannot orient it.  The lawful orientation is the fixed positive
`+1` generator and least-nonnegative `h` convention already pinned in
THM-2851.  Signed enter-minus-exit counts telescope to zero; the successful
observable is the noninvertible, one-sided positive path-semigroup count.

## 5. Why neither V4 nor a tournament sees the prize

Exact normalized cochain ranks over `F13` are

```text
dim C1=3, dim C2=9, rank B2=3, dim Z2=3, H2(V4,F13)=0.
```

Every intrinsic `mu13`-valued `V4` cocycle is therefore a coboundary.
Forgetting q-edge lengths and the positive basepoint makes every vertex
square flat.

Likewise the useful q-labelled four-cycle has edge lengths `(1,3,1,8)`.
It is a directed quiver with a based winding, not a tournament.  Pairwise
shortest-arc orientation on its four residues is transitive and reverses
or forgets the essential length-eight closing edge.  Tournament reduction
therefore destroys precisely the datum detected by `kappa`.

## 6. Sharp boundary and cheapest next test

The result pays the carry-cocycle invoice but not the carrier-state invoice.
It supplies the transition function `kappa`; it does not supply an initial
ancestry value.  A faithful object still retains

```text
(initial a in C13, q in C13, positive unit-path marker),
```

equivalently the based `C169` lift.

The next physical test is now precise:

1. refine THM-2851's oriented `q11 --9--> q7` edge into its nine positive
   unit edges;
2. attach the unique `u3` event at `q12->q0` to every one of the `449`
   pinned `QA(q11,a)->QAB(q7,a+1)` ancestry sheets;
3. retain the initial ancestry coordinate and the positive orientation;
4. on the same origin-resolved current, compare that event-decorated path
   with the direct `q3 --4--> q7` path, whose marker count is zero; and
5. require the resulting coefficient boundary, after E3/complement
   polarization, to be `omega^3` rather than the flat common-atom value.

This is the cheapest route because it uses a literal existing factor event,
not a fabricated fourth vertex.  What remains absent is a common physical
QA-to-QAB current/origin attachment across the q7 complement, not the carry
cocycle itself.

## 7. Reproduction

```bash
python3 04-computation/lrc14_endpoint_factor_exit_carry_transducer_thm2878.py
python3 -O 04-computation/lrc14_endpoint_factor_exit_carry_transducer_thm2878.py
```

Normal and optimized modes byte-match the stored output.  The script has no
executable Python `assert`.

```text
script 8e75f9c9caa9482db17b40d220f5e8ef68d00d54d4ed481869f9b31ec146f1ca
output fadccf2c9afa8cb386066b400ea51a0b90210ec26c0cd70d77b8354e279c0c38
```

The independent audit rederived the singleton danger support, all `169`
edge counts, all `2,197` composition identities and all `2,197`
marked-gauge identities.  It verified that the inherited positive
generator makes unit refinement canonical in the THM-2851 comparison
nerve, that reversal is a hostile rather than a second orientation-free
proof, and that the event realizes a carry cocycle but not a physical
current.
