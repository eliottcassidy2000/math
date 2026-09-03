---
id: THM-4386
title: "LRC14 raw-carrier lattice sum, zero-defect obstruction, and complete short shell"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373/4384 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. For primitive w, the cross product
  induces Z^3/Zw isomorphic to ker(C dot w), the mod-three owner gate and
  component length give an exact rank-two lattice sum, and two independent
  zero-defect relations force owner collision. The primitive full-support
  ternary-unit l1<=14 shell has exactly thirteen patterns; the three omitted
  by THM-4384 have sharp maxima 564/20405, 12/539, and 444/18179. This closes
  only the declared short-relation shell, not arbitrary triples or LRC(14).
source: root + lrc_relation_coverage + lrc_carrier_cleanroom / continuation session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
  - THM-4384-lrc14-small-defect-short-relation-master-formula-and-sharp-sector-atlas
primary_script: 04-computation/lrc14_raw_carrier_formula_referee_thm4386.py
primary_output: 05-knowledge/results/lrc14_raw_carrier_formula_referee_thm4386.out
primary_script_sha256: b4e4bf1a06bdb75302ee6e7dc9aeab6fa320d4bc4a5f735b4dc5f6f44366b34d
primary_output_sha256: 513963ce734a84faa220aa6f79e20cb0263fcbde2db31e5a03dc0046e9c9978d
incidence_script: 04-computation/lrc14_relation_incidence_referee_thm4386.py
incidence_output: 05-knowledge/results/lrc14_relation_incidence_referee_thm4386.out
incidence_script_sha256: 3b576101fb660c04ecf9fcc5e25205021d8c6523d421ab4adac92900db33163a
incidence_output_sha256: 8a057fa64997415c6f6c014afd34eae42e063759f4f717b9112b0b70ad8a54c1
shell_script: 04-computation/lrc14_missing_short_relation_shell_referee_thm4386.py
shell_output: 05-knowledge/results/lrc14_missing_short_relation_shell_referee_thm4386.out
shell_script_sha256: 243cc63c167e03673c00e4b141b165d02bc720c2511be2d9ac16862f59dab649
shell_output_sha256: 9612f9103d8acd8773738374d0681fcf700f2db57b5f7195ca9895917241d0cf
independent_referee_script: 04-computation/lrc14_raw_carrier_short_relation_shell_independent_referee_thm4386.py
independent_referee_output: 05-knowledge/results/lrc14_raw_carrier_short_relation_shell_independent_referee_thm4386.out
independent_referee_script_sha256: d52dc856301f4ef492e1a78f913b0451d227a486f669d6e7222da292bc2af75c
independent_referee_output_sha256: 0985a1ad10f9e82158fc451b19be72914678786a485b25c0f125b300f2c37cfb
hash_basis: raw LF bytes
audit: >
  PASS. The dependency-free clean-room referee constructs the inverse with
  n=C cross u for u dot w=1, compares the literal and lattice combs on 2,916
  triples and 15,366 components, proves the independent-relation obstruction,
  enumerates all thirteen short patterns, and rederives the three sharp
  missing-sector maxima and the H=199 presentation-incidence census. It runs
  99,241 live checks under normal, optimized, and hash-seeded Python. The
  positive-measure H=199 census, its unique uncovered peak, the 24 empty
  presentations, and the norm-16 boundary probe are producer-only
  VERIFIED-EXACT sidecars, not part of the independent audit.
---

# THM-4386 -- Raw carriers, empty double resonance, and the complete short shell

**PROVED ELEMENTARY RELATIVE TO THM-4373/4384 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) OPEN.** This theorem closes only the declared
primitive full-support ternary-unit relation shell `||c||_1 <= 14`. Universal
nonresonant control, seam entry, synchronization, all-tail transfer, and
LRC(14) remain **OPEN**.

The inherited result is THM-4384,
`01-canon/theorems/THM-4384-lrc14-small-defect-short-relation-master-formula-and-sharp-sector-atlas.md`,
relative to THM-4373,
`01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md`.
THM-4384 sharply closes ten coefficient-one relation sectors but explicitly
does not classify their union or arbitrary triples.

## Outcome

There are four main conclusions.

1. For primitive `w`, the raw carrier

   ```text
   C = w cross n
   ```

   is a complete address for a nearest-integer component:

   ```text
   Z^3 / Z w  -->  Lambda_w := {C in Z^3 : C dot w = 0}
   [n]        |->  w cross n
   ```

   is an isomorphism. The owner condition is exactly `3 does not divide
   C_1 C_2 C_3`, and there is an exact rank-two lattice-sum formula for the
   physical comb.

2. Two rationally independent integer relations which annihilate both `w`
   and the nearest-integer lift `n` force owner collision. Consequently every
   genuine double zero-defect overlap is pointwise empty. This explains all
   THM-4384 cross-sector and within-sector overlaps and shows why sector
   maxima cannot be added.

3. The coefficient-one atlas misses three relation directions at the same
   short-defect scale. The complete primitive, full-support, ternary-unit
   carrier shell `||c||_1 <= 14` has thirteen magnitude patterns, not ten.
   The missing patterns are

   ```text
   (2,5,5), (2,5,7), (4,5,5).
   ```

   Their exact sharp maxima are respectively

   ```text
   564/20405 at {5,53,55};
   12/539    at {7,17,77} and {7,53,77};
   444/18179 at {5,49,53}.
   ```

   Thus THM-4384 plus these three rows sharply closes the entire declared
   `l1 <= 14` full-support ternary-unit relation shell. Its overall maximum
   remains the THM-4384 value `6/77` at `{1,5,11}`.

4. In the finite hostile universe of all primitive unordered distinct
   positive odd three-unit triples with maximum speed at most 199, the ten
   THM-4384 sectors cover only 5,674 of 47,499 triples. The first uncovered
   height is 23, with four exact positive-comb hostiles. This is presentation
   incidence, not phase-space coverage, and it exposes the rank-two
   nonresonant lattice sum as the real remaining object.

## Inheritance pass and concept board

The closest proved mechanism is THM-4384's bounded integer defect followed by
one-ray quadrature. The canonical hostile is the unrestricted rank-two
nearest-integer lift lattice of a primitive triple. The corrected near miss is
to treat a coefficient presentation as an intrinsic partition: the same speed
triple can have several presentations, and a small relation containing a
coefficient divisible by three need not have defect zero. The least-used
sidecar is the un-normalized cross product `C=w cross n`; normalizing it to a
primitive direction discards the component scale.

The live board was:

```text
raw carrier C; primitive carrier direction; owner residues;
integer relation defects; coefficient-pattern incidence;
rank-two lattice geometry; first defect-three shell.
```

The source-to-target map is now explicit. A lift class `[n]` maps to `C`,
preserving the exact component length and owner gate. Normalizing `C` preserves
only the zero-defect relation direction and destroys the integer scale which
addresses the component. A coefficient-sector label preserves the existence
of a relation but destroys presentation multiplicity and does not partition
either speed triples or phase space.

## 1. Raw carrier theorem

Let `w=(w_1,w_2,w_3)` be primitive. Put

```text
Lambda_w = {C in Z^3 : C dot w = 0}.
```

### Proposition 1.1: integral isomorphism

The map

```text
Phi_w : Z^3/Zw -> Lambda_w,       [n] |-> w cross n
```

is an isomorphism.

**Proof.** The cross-product matrix is

```text
[  0   -w_3   w_2 ]
[ w_3    0   -w_1 ]
[-w_2   w_1    0  ].
```

Its rank is two. The gcd of its entries is `gcd(w_1,w_2,w_3)=1`, and the gcd
of its `2 x 2` minors is `gcd(w_1,w_2,w_3)^2=1`. Hence its Smith normal form is
`diag(1,1,0)`. Its kernel is the rational line `Qw` intersected with `Z^3`,
which equals `Zw` because `w` is primitive. Its image is a saturated rank-two
sublattice of `w^perp`; therefore it is all of `Lambda_w`. QED.

### Proposition 1.2: exact owner gate

At a physical body phase use the nearest-integer lift

```text
n_i = nint(w_i y),       e_i=w_i y-n_i,       |e_i|<r,
r=3/14,
o_i=-w_i^(-1)n_i mod 3.
```

For `{i,j,k}={1,2,3}`,

```text
C_i=w_j n_k-w_k n_j
   =w_j w_k(o_j-o_k) mod 3.
```

Since each `w_i` is a ternary unit, the owners are pairwise distinct if and
only if every coordinate of `C` is nonzero modulo three:

```text
{o_1,o_2,o_3}=F_3  <=>  3 does not divide C_1 C_2 C_3.       (1)
```

### Proposition 1.3: exact component length and comb sum

For `C in Lambda_w`, define

```text
L_w(C)=max(0,min(
  2r/w_1, 2r/w_2, 2r/w_3,
  r/w_1+r/w_2-|C_3|/(w_1 w_2),
  r/w_1+r/w_3-|C_2|/(w_1 w_3),
  r/w_2+r/w_3-|C_1|/(w_2 w_3))).                         (2)
```

Then, in the normalization of THM-4373/4384,

```text
mu(F_w) = sum_(C in Lambda_w, 3 does not divide C_1 C_2 C_3) L_w(C).   (3)
```

The sum is finite. Indeed, positivity of the three pair terms gives explicit
strict bounds on the three carrier coordinates.

**Proof.** For a lift `n`, the three eligible real intervals have centers
`n_i/w_i` and radii `r/w_i`. The length of their common intersection is the
minimum of every upper endpoint minus every lower endpoint, clipped at zero.
The three same-index differences give `2r/w_i`; combining the two directed
differences for a pair gives the corresponding absolute-value term in (2).
Proposition 1.1 says that every lift modulo physical translation has exactly
one raw carrier and every carrier occurs. Equation (1) is exactly the physical
owner gate. QED.

The primitive direction `C/gcd(C_1,C_2,C_3)` is the unique primitive
zero-defect relation attached to a physical component: any integer vector
annihilating both `w` and `n` lies on the one-dimensional common orthogonal
line. The raw scale is not redundant; it is the component address.

The independent referee compared (3) with a literal interval-comb
implementation on all 2,910 primitive triples through height 79 and six high
controls: 2,916 triples and 15,366 base components in total. Exact dictionaries
of raw carriers and lengths agreed component by component.

## 2. Short relations and the pointwise overlap obstruction

### Lemma 2.1: full-support short relations have zero defect

Let `w` be a positive odd three-unit triple, and let `c in Z^3` satisfy

```text
c dot w=0,
c_i is nonzero mod 3 for i=1,2,3,
||c||_1 <= 14.
```

At every alleged physical failure phase, with lift `n` as above,

```text
c dot n=0.                                               (4)
```

**Proof.** Put `delta=c dot n=-c dot e`. Strict eligibility gives

```text
|delta| < r ||c||_1 <= 3.                                (5)
```

Set `A_i=c_i w_i mod 3`. Every `A_i` is nonzero, while
`A_1+A_2+A_3=c dot w=0`. Three nonzero elements of `F_3` sum to zero only when
they are all equal, say to `A`. Because `n_i=-w_i o_i mod 3` and the distinct
owners are `{0,1,2}`,

```text
delta = -A(o_1+o_2+o_3)=0 mod 3.
```

The only multiple of three satisfying (5) is zero. QED.

This is the true zero-defect criterion. A coefficient equal to one is not
required.

### Lemma 2.2: independent zero-defect relations force an empty comb

Let `w in Z^3` be primitive and have ternary-unit coordinates. Suppose an
alleged physical component has nearest-integer lift `n in Z^3`. If
`c,d in Z^3` are `Q`-linearly independent and

```text
c dot w=d dot w=c dot n=d dot n=0,                       (6)
```

then that component cannot have distinct owners.

**Proof.** Both `w` and `n` lie in the one-dimensional common orthogonal of
`c,d`, so `n=lambda w` for some rational `lambda`. Primitivity gives an
integer Bezout vector `z` with `z dot w=1`; hence
`lambda=z dot n` is an integer. Therefore

```text
o_i=-w_i^(-1)n_i=-lambda mod 3
```

for every `i`, and all owners collide. QED.

The quantifiers are pointwise: whenever two relations vanish on one alleged
component, that component is impossible. If two fixed relations meet the
hypotheses of Lemma 2.1, they vanish on every alleged component, so the entire
physical comb is empty. This is not an instruction to add or subtract sector
maxima.

### Corollary 2.3: one-ray formula

If primitive `c` satisfies Lemma 2.1 and `c dot w=0`, then (4) makes
`C=w cross n` parallel to `c`. Primitivity makes the multiplier integral, and
the owner gate excludes precisely the multiples of three. Thus

```text
mu(F_w)=sum_(k in Z, 3 does not divide k) L_w(kc)
       =2 sum_(k>=1, 3 does not divide k) L_w(kc).        (7)
```

Unlike the coefficient-one endpoint formula, (7) works without choosing a
distinguished endpoint or requiring a unit coefficient.

## 3. Complete `l1 <= 14` coefficient shell

For a primitive physical carrier direction, full support and (1) require all
three coefficient magnitudes to be nonzero modulo three. Odd speeds and
`c dot w=0` imply that the coefficient sum is even modulo two, equivalently
`||c||_1` is even. Sorting magnitudes and imposing gcd one leaves exactly
thirteen patterns through norm 14:

```text
l1=4:  (1,1,2)
l1=6:  (1,1,4)
l1=8:  (1,2,5)
l1=10: (1,1,8), (1,2,7), (1,4,5)
l1=12: (1,1,10), (1,4,7), (2,5,5)
l1=14: (1,2,11), (1,5,8), (2,5,7), (4,5,5).            (8)
```

THM-4384 covers the ten rows containing a coefficient one. The first missing
coefficient complexity is therefore `(2,5,5)` at norm 12, not a norm-16
coefficient-one relation. Three minimal concrete hostiles are:

| role | speed triple | primitive carrier direction | exact measure |
|---|---|---|---:|
| first missing pattern by norm | `{1,11,25}` | `(5,-5,2)` | `38/1925` |
| remaining `(2,5,7)` pattern | `{5,17,25}` | `(7,-5,2)` | `4/425` |
| first ten-sector miss by height | `{5,19,23}` | `(4,5,-5)` | `2/665` |

The last triple also has the short relation `(3,-2,1)`, but its physical
defects for that presentation are `-1,+1`; a coefficient divisible by three
breaks the ternary zero-defect argument. Likewise `(5,17,25)` has relation
`(1,10,-7)` with physical defects `-3,+3`, while its true zero-defect carrier
is `(7,-5,2)`.

## 4. Sharp maxima for the three missing rays

Fix a primitive short relation direction `c` and set

```text
f(t)=L_w(tc),       t in R.
```

Formula (2) shows that `f` is even, compactly supported, and nonincreasing on
the nonnegative axis. Its integral is independent of the particular speed
triple on `c dot w=0`:

```text
I(c)=integral_R f(t) dt
    =area(c^perp intersect [-r,r]^3)/||c||.              (9)
```

For the magnitude pattern `(a,b,c)`, elementary three-box convolution gives

```text
I(a,b,c)=r^2/(2abc) sum_(S subset {1,2,3}) (-1)^|S|
  (a+b+c-2 sum_(i in S) a_i)_+^2.                       (10)
```

For `S_h=sum_(k in Z) f(hk)`, monotonicity gives the elementary lattice-rule
error bound

```text
|S_h-I/h| <= f(0).
```

Since the owner gate makes `mu=S_1-S_3` and
`f(0)=2r/w_max=3/(7w_max)`, one obtains

```text
mu(F_w) <= (2/3)I(c)+6/(7w_max).                         (11)
```

After a candidate `M` is found, (11) gives a rigorous finite maximum-speed
cutoff. Exhaustive exact evaluation of (7) below that cutoff, independently
compared to the literal physical comb on every retained triple, gives:

| pattern | `I(c)` | baseline `2I/3` | strict `w_max` threshold | finite triples | sharp maximum | maximizer(s) |
|---|---:|---:|---:|---:|---:|---|
| `(2,5,5)` | `81/2450` | `27/1225` | `204050/1333` | 208 | `564/20405` | `{5,53,55}` |
| `(2,5,7)` | `9/343` | `6/343` | `539/3` | 442 | `12/539` | `{7,17,77}`, `{7,53,77}` |
| `(4,5,5)` | `36/1225` | `24/1225` | `64925/366` | 252 | `444/18179` | `{5,49,53}` |

Thus the finite enumerations use respectively `w_max <=153,179,177`; every
larger triple is excluded analytically by (11). In particular the `(2,5,7)`
proof is not based on the earlier discovery scan through 399. A separate
optimization-safe exact replay through **maximum speed 1200** (strictly beyond
539 even in the literal height coordinate) found 20,784 `(2,5,7)` triples and
the same two maximizers. That replay is corroboration; inequality (11) is the
completeness proof.

The sharp referee checks 902 sector triples in the three finite proof windows
and 6,602 literal base components. Normal and optimized runs have identical
stdout.

## 5. Exact incidence versus union coverage through height 199

The hostile universe is

```text
U_H={unordered distinct positive odd triples w:
     3 does not divide w_1w_2w_3, gcd(w)=1, max(w)<=H}.
```

For this section, a triple has **presentation incidence** in sector `(h,m)` if
some signed permutation of `(1,h,m)` annihilates it. Its **ten-sector union
coverage** means incidence in at least one of the ten THM-4384 sectors. Neither
notion partitions the physical `y` circle, and neither counts components of
`F_w`.

At `H=199` the exact counts are:

```text
universe                                      47,499
positive physical comb                       47,463
ten-sector presentation union                 5,674
ten-sector uncovered                         41,825
positive and presentation-covered             5,650
positive and presentation-uncovered          41,813
presentation-covered but empty comb              24
```

The number of distinct sector labels per triple is

```text
0 labels: 41,825
1 label:   5,663
2 labels:     10
3 labels:      1.
```

Individual sector incidences, which deliberately sum to more than the union,
are

```text
(1,2):1023  (1,4):515  (1,8):257  (1,10):203
(2,5):825   (2,7):588  (2,11):370
(4,5):800   (4,7):589  (5,8):516.
```

Exact coverage by height, in columns
`(all, presentation-covered, positive, positive-and-covered)`, is

| `H` | counts |
|---:|---:|
| 11 | `(4,4,2,2)` |
| 23 | `(56,52,41,37)` |
| 43 | `(454,245,422,221)` |
| 79 | `(2910,883,2874,859)` |
| 127 | `(12231,2287,12195,2263)` |
| 199 | `(47499,5674,47463,5650)` |

The earliest uncovered height is 23, with exactly

```text
{5,19,23}:   mu=2/665
{7,19,23}:   mu=4/437
{11,19,23}:  mu=4/437
{17,19,23}:  mu=10/2261.
```

Among uncovered triples through height 199, the unique largest measure is

```text
mu(F_(1,11,175))=36/1225.                               (12)
```

This is **FINITE-EXACT**, not a universal nonresonant theorem.

## 6. Complete THM-4384 overlap atlas and its meaning

The complete cross-sector presentation-overlap list is retained here:

```text
(1,4)-(5,8):   {1,7,11}
(1,8)-(2,5):   {1,5,13}
(1,8)-(4,7):   {1,5,13}, {1,11,19}
(1,10)-(2,5):  {1,7,17}
(1,10)-(4,7):  {1,13,23}
(2,5)-(4,7):   {1,5,13}, {5,7,11}
(2,7)-(2,11):  {1,7,25}
(2,7)-(5,8):   {1,5,17}, {7,11,19}
(2,11)-(5,8):  {5,7,31}, {7,19,29}.
```

The within-sector independent duplicate-presentation rays are

```text
(2,7):  {1,5,17}
(2,11): {5,7,41}
(5,8):  {11,13,23}.
```

There are thirteen distinct double-relation rays in total:

```text
{1,5,13}, {1,5,17}, {1,7,11}, {1,7,17}, {1,7,25},
{1,11,19}, {1,13,23}, {5,7,11}, {5,7,31}, {5,7,41},
{7,11,19}, {7,19,29}, {11,13,23}.
```

Lemma 2.2 proves, without enumeration, that every one of these physical combs
is empty. The referee independently checks all thirteen directly.

This distinguishes three things that must not be conflated:

- a **presentation** is an algebraic label on a speed triple;
- the **union** is a set of speed triples having at least one label;
- `F_w` is one physical subset of the phase circle for that triple.

For a triple in several sectors, each exact sector formula computes the same
`F_w`; it does not produce a new piece to add. Even for disjoint subclasses of
speed triples, a universal pointwise bound over their union takes the maximum
of the subclass bounds, not their sum. The canonical overlap hostile
`{1,5,13}` lies in three THM-4384 sectors and has empty comb. If the coefficient
list is enlarged to include small nonzero-defect presentations, `{1,5,11}`
has four labels `(1,2),(1,6),(3,4),(4,9)` but still just one comb, of measure
`6/77`.

For completeness, the 24 ten-sector-covered empty triples through height 199
are

```text
{1,5,13}, {1,5,17}, {1,5,31}, {1,5,35}, {1,7,11},
{1,7,17}, {1,7,25}, {1,11,19}, {1,13,19}, {1,13,23},
{1,13,37}, {1,17,23}, {5,7,11}, {5,7,31}, {5,7,41},
{5,11,13}, {5,13,17}, {5,13,29}, {5,17,19}, {7,11,19},
{7,19,29}, {11,13,23}, {11,17,25}, {13,17,19}.
```

## 7. The next boundary: norm 16 needs affine defect layers

The next primitive, full-support, ternary-unit shell has exactly seven
magnitude patterns:

```text
(1,1,14), (1,2,13), (1,4,11), (1,5,10),
(1,7,8), (2,7,7), (4,5,7).                              (13)
```

Now eligibility gives only

```text
|delta| < (3/14)16=24/7,
delta=0 mod 3,
```

so `delta` may be `-3,0,+3`. Every pattern in (13) realizes nonzero defect in
the exact `H=199` probe. Earliest positive examples are:

| pattern | earliest triple | relation | observed defects | measure |
|---|---|---|---|---:|
| `(1,1,14)` | `{1,5,19}` | `(14,1,-1)` | `{-3,3}` | `4/133` |
| `(1,2,13)` | `{1,5,23}` | `(13,2,-1)` | `{-3,3}` | `2/161` |
| `(1,4,11)` | `{1,7,19}` | `(1,-11,4)` | `{-3,3}` | `8/931` |
| `(1,5,10)` | `{5,7,13}` | `(1,-10,5)` | `{-3,3}` | `4/637` |
| `(1,7,8)` | `{17,23,25}` | `(8,-7,1)` | `{-3,0,3}` | `744/68425` |
| `(2,7,7)` | `{7,11,13}` | `(2,7,-7)` | `{-3,3}` | `2/1001` |
| `(4,5,7)` | `{5,13,25}` | `(7,5,-4)` | `{-3,3}` | `2/2275` |

These are **FINITE-EXACT boundary witnesses**, not sharp global maxima. The
one-ray reduction is no longer complete; the next theorem must retain the
affine defect state `delta/3` together with the raw carrier.

There is also a useful coefficient-one side split below `h+m<=13`. Besides
the ten zero-defect patterns, exactly ten parity-compatible patterns have one
of `h,m` divisible by three:

```text
(1,6), (1,12), (2,3), (2,9), (3,4),
(3,8), (3,10), (4,9), (5,6), (6,7).
```

They force a nonzero residue class for the relation defect rather than zero.
The sole pair `(3,6)` is impossible modulo three. This side split explains why
blindly adjoining all small coefficient-one presentations does not extend the
zero-defect theorem.

## 8. Consequence for the open frontier

The raw-carrier formula changes the universal nonresonant problem into the
explicit lattice discrepancy problem

```text
sum_(C in Lambda_w, all C_i nonzero mod 3) L_w(C),        (14)
```

where `L_w` is the compact anisotropic tent (2). A short zero-defect relation
collapses the contributing part of this rank-two lattice to one ray. Outside
that case, primitive direction alone is insufficient: several directions and
raw scales can contribute, as the twelve carriers for `{1,11,175}` show.

The sharp next tasks suggested by the packet are therefore:

1. prove a geometry-of-numbers or residue-discrepancy bound for (14) under the
   absence of a ternary-unit primitive vector of norm at most 14;
2. build an affine `delta/3 in {-1,0,1}` quadrature for the seven norm-16
   patterns rather than forcing a false zero-defect reduction;
3. use Lemma 2.2 as a seam obstruction: simultaneous independent zero-defect
   certificates eliminate a phase instead of contributing two bounds;
4. test whether successive minima of the support hexagon of `L_w`, with the
   three coordinate-zero residue sublattices removed, control (14) uniformly.

The `H=199` census says the ten original sectors are far from a coverage
argument: 41,813 positive combs remain outside them. The new three sectors
complete a natural coefficient shell, not the arbitrary-triple universe.

## Reproduction and frozen artifacts

All four checked-in referees use explicit runtime exceptions and remain active
under `python -O`. Normal, optimized, and fixed-hash-seed invocations produce
the checked-in streams:

```powershell
python -B 04-computation/lrc14_raw_carrier_formula_referee_thm4386.py
python -B 04-computation/lrc14_relation_incidence_referee_thm4386.py
python -B 04-computation/lrc14_missing_short_relation_shell_referee_thm4386.py
python -B 04-computation/lrc14_raw_carrier_short_relation_shell_independent_referee_thm4386.py
```

The canonical SHA-256 values are recorded in the front matter. The independent
referee imports no repository mathematics and shares no producer code. Its
scope is exactly the one stated in the audit field; the additional
positive-measure and norm-16 observations in the incidence output remain
producer-only **VERIFIED-EXACT** evidence.
