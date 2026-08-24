# Imaginary quadratic class-rank challenge: source, data, and certificate audit

> **SOURCE / DATA / CERTIFICATE AUDIT, 2026-08-24.**  Every status in this
> note is explicit.  In particular, **no line in the Epoch class-rank
> challenge is certified here**.  The displayed form lines are exact
> below-threshold controls and earn no challenge credit.  PARI
> `quadclassunit` invariant factors are used only for GRH-assisted discovery
> unless an unconditional certification is stated.

## Executive verdict

The challenge remains open at every requested pair

```text
ell = 3:   r = 9,...,16
ell = 5:   r = 6,...,12
ell = 7:   r = 5,...,10
ell = 11:  r = 4,...,8
ell = 13:  r = 4,...,8.
```

Targeted searches of the public problem page, primary papers, the associated
public code and data, GitHub, and exact certificate-like strings found no
qualifying line.  This is not a proof that no unpublished example exists.  It
agrees with the problem setter's statement that no `D` is known for the tested
pairs.

The useful outcomes are instead:

1. a tiny unconditional verification path that never computes the full class
   group;
2. exact form bases at the known `3,5,7,11,13` near-frontiers;
3. an unconditionally certified one-digit source error that would otherwise
   manufacture fictitious `11`-torsion;
4. a complete census of the only public `131199`-field Bagshaw data file;
5. a concrete request for the much larger unpublished form buckets; and
6. two typed discovery routes: Mordell-curve descent for `ell=3` and
   dihedral-field collision buckets for general odd `ell`.

## 1. The problem specification has a discriminant/radicand inconsistency

The controlling specification is Epoch AI's institutional 2026 web problem,
[*Large ell-Rank in Class Groups of Imaginary Quadratic
Fields*](https://epoch.ai/frontiermath/open-problems/class-group-rank).  No
individual mathematical author is named there.

The introductory prose asks for a square-free positive `D` in
`Q(sqrt(-D))`.  The actual submission grammar instead requires a positive
integer `D` such that `Delta=-D` is a **fundamental discriminant**.  The latter
must control because the verifier consumes binary forms of discriminant
`-D`, and the page's own example has

```text
D = 4447704 = 2^3 * 3 * 47 * 3943,
```

which is not square-free while `-4447704` is fundamental.  Thus `D` is the
absolute field discriminant, not necessarily a square-free radicand.

The page also attributes to an unlinked `[Kim 2015]` the statement that for
every `ell` there are infinitely many imaginary quadratic fields of
`ell`-rank at least three.  A targeted search did not identify a matching
primary source.  Moreover, the 2024 Bagshaw--Jacobson--Scheidler--Rollick
introduction says that examples of odd-prime rank exceeding two were known
only for primes at most `19`.  The Epoch sentence is therefore
**UNAUDITED / LIKELY CITATION ERROR** and is not used below.

## 2. The strongest audited mechanisms

### 2.1 The rank-eight `3`-frontier is theorem-backed

Noam Elkies's [*Rank of an elliptic curve and 3-rank of a quadratic field via
the Burgess bounds*](https://doi.org/10.1007/s40993-024-00601-x), *Research in
Number Theory* **11** (2025), article 70, proves for

```text
D0 = 72513834653847828539450325493
Delta = -3 D0 = -217541503961543485618350976479
```

that

```text
Cl(Delta) = (Z/2)^2 x (Z/3)^8 x Z/77681 x Z/139939.       (CITED, PROVED)
```

In the paper, Theorem 1 treats the rank-sixteen Mordell curves,
Proposition 2 gives

```text
rank E_(16D)(Q) <= 1 + r3 Cl(Q(sqrt(D)))
                         + r3 Cl(Q(sqrt(-3D))),
```

and Proposition 3 removes the GRH assumption from the imaginary class-group
calculation.  Footnote 3 reports reduction of all `3^8-1=6560` nonzero
generator combinations and explicitly notes the `2(3^4-1)=160`
meet-in-the-middle alternative.

### 2.2 The best public general search produces order-`p` forms directly

Christian Bagshaw, Michael J. Jacobson, Renate Scheidler, and Nickolas
Rollick, [*Improved methods for finding imaginary quadratic fields with high
`n`-rank*](https://doi.org/10.1090/conm/796/15995), *Contemporary Mathematics*
**796** (2024), use

```text
4 m^p = y^2 - z^2 Delta.                                  (1.2)
```

Their Proposition 2.3(a) gives sufficient gcd and reduction conditions for
the associated ideal to have exact prime order `p`.  Theorem 2.4 constructs
the ideal basis

```text
{m, (x + sqrt(Delta))/2},
```

and hence the binary form

```text
[m, x, (x^2-Delta)/(4m)].
```

Algorithm 3.2 groups many norm-equation solutions by `Delta` and tests the
corresponding forms for independence.  Section 3.2 records the crucial
failure boundary: these discovered solution forms rarely detect rank above
two even when the full class group has larger `p`-rank.  Table 5.1 identifies
full class-group computation as the bottleneck.

The [public implementation](https://github.com/ChristianBagshaw/quadratic-fields-high-n-rank)
is narrower than the paper's raw computation:

* `src/rankcheck.sage` returns as soon as it detects a subgroup of order
  `p^2`;
* the public `data/D.txt` contains only `131199` discriminants, without their
  form buckets; and
* the legacy search code documents raw lines of the form
  `Delta,[A,B,C]`, but the billions-scale raw files are not present.

### 2.3 Appendix exponent tuples are not ranks

In the paper's Appendix, a tuple `(e1,...,ek)` denotes

```text
C(p^e1) x ... x C(p^ek).
```

The `p`-rank is the tuple **length** `k`.  For example,
`(6,1,1,1)` is `5`-rank four, not rank six, and `(9,1,1,1,1)` is
`3`-rank five, not rank nine.  Reading the leading exponent as the rank is a
seductive but invalid route to a challenge line.

## 3. Exact below-threshold form controls

The following lines were derived from public discriminants using PARI only
as a generator-discovery aid.  They then passed independent exact integer
checks of fundamentality, primitivity, discriminant, exact prime order, and
independence under Gauss composition.  Their spans have sizes
`3^8,5^4,7^4,11^3,13^3`, respectively.

```text
3|8|217541503961543485618350976479|133663175160114,-83514707956161,419929068229900|85450831602870,36334851996309,640314790756777|10768953543577,-3271456475729,5050449087478240|2688172755838,2248218685611,20231824570830525|183294999722469,-120074312826153,316374348609388|215222692946620,-204653468437719,301344322191753|98225989117273,-57322363599763,562039026826594|59581187049667,-59327762944611,927563289538350
5|4|1264381632596|188505,25402,1677710|415025,-353402,836862|38025,-6952,8313149|318873,-67684,994881
7|4|469874684955252968120|1958128126,876510532,60088375911|6723992369,-3407420908,17901760434|2775584655,-2184045520,42751780146|1369,-1304,85806187902712686
11|3|3035884424|81,-52,9370022|625,24,1214354|11025,3674,69147
13|3|38630907167|89199,41411,113078|557,469,17338926|4466,-3359,2163132
```

They are **FINITE-EXACT lower bounds**, but all lie below the challenge set
and are ignored for scoring.  They are collected separately in
[`imaginary_quadratic_class_rank_below_threshold_controls_20260824.txt`](../05-knowledge/results/imaginary_quadratic_class_rank_below_threshold_controls_20260824.txt).
The exact verifier, reproduction commands, and independent PARI comparison
are in
[`imaginary_quadratic_class_rank_certificate_tool_20260824.out`](../05-knowledge/results/imaginary_quadratic_class_rank_certificate_tool_20260824.out).

## 4. A one-digit `11`-rank source error is an exact hostile

Bagshaw et al. Table 5.2 prints

```text
-3532321517865683
```

as the smallest `11`-rank-at-least-three discriminant found in that search,
whereas Appendix Table A.4 prints

```text
-3532321517864683.
```

The distinction is load-bearing.

For the Table 5.2 value, `bnfinit` returned invariant factors
`[3555056,2,2]`, and the unconditional command `bnfcertify(bnf)` returned
`1`.  Hence

```text
Cl(-3532321517865683) = Z/3555056 x Z/2 x Z/2,
h = 14220224 = 2^6 * 83 * 2677,
```

and there is **no `11`-torsion**.  This negative statement is
**FINITE-EXACT**, not merely GRH-assisted.  The command and output are frozen
in the [unconditional hostile ledger](../05-knowledge/results/bagshaw_11rank_typo_bnfcertify_20260824.out).

For the Appendix value, exact form arithmetic verifies the positive lower
bound

```text
11|3|3532321517864683|9266797,1918973,95394449|18021997,-936603,49012309|19135921,-2565769,46233791
```

with all `11^3=1331` products distinct.  Therefore the digit `5` in Table
5.2 is wrong and the digit `4` in Appendix A.4 is right.

There is a second, milder typo in the same paper: introductory prose drops
the final `6` from the `5`-rank-four field.  Tables 5.2/A.2 and the exact
four-form control establish the fundamental discriminant
`-1264381632596`.

## 5. The cheapest reliable certificate is not a class-group computation

Given a proposed line, set `Delta=-D` and perform the following exact checks.

1. Verify that `Delta` is fundamental.  Equivalently, `Delta` is square-free
   and `Delta=1 mod 4`, or `Delta=4d` with square-free `d=2 or 3 mod 4`.
   For a very large constructed `D`, retain a factorization and primality
   certificates as a sidecar.
2. For every `(a,b,c)`, check `a>0`, `gcd(a,b,c)=1`, and
   `b^2-4ac=Delta`, then reduce the form canonically.
3. Use principal form `(1,1,(1-Delta)/4)` for odd `Delta`, and
   `(1,0,-Delta/4)` for even `Delta`.
4. Check that each class is nonprincipal and that its `ell`th power is
   principal.  Since `ell` is prime, its order is exactly `ell`.
5. Split the `r` forms into sets of sizes `u` and `v`.  Hash all elementary
   `ell`-products in each half.  Each half must have its full expected size,
   and the two subgroup hashes may intersect only in the principal class.

For the largest challenge ranks, the largest half contains only

```text
3^8 = 6561,   5^6 = 15625,   7^5 = 16807,
11^4 = 14641, 13^4 = 28561
```

reduced forms.  This proves the requested lower bound without a class number,
a full class group, a GRH assumption, or enumeration of all `ell^r`
relations.  PARI/Magma may be used to discover candidate generators, but a
final claim should depend only on this exact replay.

The committed implementation is
[`imaginary_quadratic_class_rank_certificate_tool_20260824.py`](../04-computation/imaginary_quadratic_class_rank_certificate_tool_20260824.py).

## 6. What the public data actually contains

The public Bagshaw `data/D.txt` was exhaustively rescanned.  PARI's
GRH-assisted invariant factors report

```text
5-rank 2: 129889 fields
5-rank 3:   1309 fields
5-rank 4:      1 field, D=1264381632596
5-rank 5+:     0 fields.
```

This full-rank histogram is **FINITE-COMPUTED / GRH-ASSISTED DISCOVERY**;
`bnfcertify` was not run field by field.  The four displayed generators for
the unique maximum independently prove its rank-at-least-four lower bound.
See
[`bagshaw_public_5rank_census_20260824.out`](../05-knowledge/results/bagshaw_public_5rank_census_20260824.out).

Table 5.1 records a much larger, specially biased search population:

| `p` | discriminants found | full class groups computed | uncomputed |
|---:|---:|---:|---:|
| `3` | `20609841975` | `20609841975` | `0` |
| `5` | `1331448842` | `1331448842` | `0` |
| `7` | `402708300` | `297354233` | `105354067` |
| `11` | `13236853` | `10342190` | `2894663` |
| `13` | `5013641` | `2522501` | `2491140` |

The high-value external request is therefore not merely a list of `D` values.
It is the raw sorted bucket for each `D`, including every form or
`(m,y,lambda z)` triple, bucket multiplicity, the uncomputed discriminants,
and retained factorization data.  Those files would permit immediate exact
incremental rank tests and avoid repeating thousands of CPU-days.

## 7. Why random PARI search is not competitive

The Bagshaw populations are already biased to have at least two explicit
`p`-classes.  Nevertheless, more than a billion fully computed `p=5` fields
produced rank four as the top reported tier, hundreds of millions of
computed `p=7` fields produced only `67` rank-at-least-four fields, and the
`p=11,13` runs stopped at rank three.  A random-discriminant scan discards
this useful bias and is strictly less targeted.

The Cohen--Lenstra distribution gives a heuristic rank-`r` scale of roughly
`p^(-r^2)` up to `p`-dependent factors.  This is **HEURISTIC**, but it
correctly identifies rarity, rather than the subsecond cost of many
individual `quadclassunit` calls, as the dominant barrier.

The source-error hostile makes the proof-economics contrast concrete:
exact certification of three supplied `11`-classes is tiny, whereas
unconditionally certifying the entire erroneous class group took about
`231` seconds.  Discovery and acceptance should be separate programs.

## 8. Three live frontier routes

### Anchor: a rank-eighteen clean Mordell curve for `(3,9)`

Under Elkies Proposition 2's square-free and congruence hypotheses, an
exhibited rank-at-least-eighteen curve `E_(16D)` forces the two mirror
`3`-ranks to sum to at least `17`.  Scholz reflection then forces the
imaginary mirror to have rank at least nine.  A rank-seventeen curve does not
have this consequence; Elkies Section 3.4 explicitly reports a rank-seventeen
example whose mirror fields still have ranks `8` and `7` under GRH.

Thus the clean computational target is not a random quadratic field but a
Mestre--Nagao search for a rank-eighteen Mordell curve, followed by
3-isogeny descent, GRH-assisted form extraction, and the unconditional
certificate replay above.

### Niche: retain projective form novelty in the Bagshaw buckets

The public code should be conceptually extended from the predicate
`rank>=2` to incremental meet-in-the-middle membership.  Every new reduced
order-`p` form is tested against the current span, and buckets are scored by
the number of distinct projective order-`p` lines rather than by raw
norm-equation solution count.  Different powers of one class and repeated
solutions otherwise create a large but rank-poor bucket.

Section 3.2's failure warning remains: even a rich norm-equation bucket may
miss true class-group directions.  Full `p`-Sylow computation should be
reserved for the highest novelty/multiplicity buckets, while the
`lambda`-pairs and ideal-norm supports are varied to expose missing
directions.

### Wildcard: collide dihedral fields by quadratic resolvent

Class field theory gives a typed dual search.  If a quadratic field has
`p`-rank `r`, then it has

```text
(p^r-1)/(p-1)
```

unramified degree-`p` extensions.  Equivalently, Bagshaw et al. explain that
there are that many degree-`p` fields with dihedral Galois closure sharing
the quadratic resolvent and discriminant power.  At the minimum challenge
ranks, the required collision multiplicities are already

```text
(3,9): 9841,  (5,6): 3906,  (7,5): 2801,
(11,4): 1464, (13,4): 2380.
```

The map is

```text
nonzero Hom(Cl,Cp) / Fp^x  ->  unramified Cp extension
                            ->  degree-p dihedral field.
```

It preserves the number of projective `p`-directions but destroys a chosen
generator basis and scalar labels.  Explicit ideal classes or binary forms
are therefore the necessary sidecar.  A database indexed by quadratic
resolvent offers an orthogonal collision search, not a certificate by itself.

## 9. Final claim ledger

* **CITED / PROVED:** Elkies's exact `3`-rank-eight class group and the
  Mordell/reflection mechanism.
* **CITED:** Bagshaw et al. Proposition 2.3, Theorem 2.4, Algorithm 3.2,
  computational tables, and the failure/bottleneck observations.
* **FINITE-EXACT:** the five below-threshold form bases; the corrected
  Appendix `11`-rank-three basis; and the absence of `11`-torsion at the
  erroneous Table 5.2 discriminant.
* **FINITE-COMPUTED / GRH-ASSISTED:** the full invariant-factor census of the
  public `131199`-field data file and any uncertified negative rank statement
  from `quadclassunit`.
* **HEURISTIC:** Cohen--Lenstra rarity estimates and bucket-priority advice.
* **OPEN:** every requested challenge pair.  No qualifying output line is
  claimed.
