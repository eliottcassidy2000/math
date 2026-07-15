---
id: THM-822
title: Height-two Hamming-five closure and the live-relation kernel boundary
status: PROVED FINITE-EXACT (25,344 maximin rows and exact three-face kernel joins)
source: codex-2026-07-15-S10 Hamming-five kernel continuation
depends_on:
  - THM-806  # exact half-open owner-handoff bands
  - THM-820  # Hamming-five chart and exact maximin method
related: [THM-815, THM-818, THM-821, THM-845, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_height_two_kernel_pairs_codex_S10.cpp
  - 05-knowledge/results/lrc13_hamming_five_height_two_kernel_pairs_codex_S10.out
---

# THM-822 — Height-two Hamming-five kernel boundary

Put

```text
delta=1/13,                         [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

Here “height two” means that every replacement has lift height **at most**
two, not that all five heights equal two.

## Theorem

Let `R={r_1<...<r_5}` be a five-element subset of `[12]`, put
`P=[12] minus R`, and choose

```text
h_i in {1,2},                      u_i=r_i+13h_i.
```

Then the scale-one five-replacement packet

```text
B=P union {u_1,...,u_5}                                  (1)
```

is loose:

```text
M(B)>1/13.                                                (2)
```

Across all

```text
C(12,5) 2^5=25,344                                      (3)
```

rows, the unique minimum is

```text
M(B)=1/12,
R=(1,3,5,7,9),
(h_1,...,h_5)=(1,1,1,1,1),
B=2[11] union {11}.                                     (4)
```

Its complete maximizer set is

```text
{1,5,7,17,19,23}/24.                                    (5)
```

The same census gives an exact boundary for a proposed three-face quotient.
The labelled live owner-handoff relation `H0` and the relation decorated by
all eligible integer centres `H1` have identical fibres.  Their ordered
three-face kernel join has

```text
relation rows                         111,006
off-diagonal rows                      85,662
collision fibres                        5,855
largest fibre                              20
fibres mixing exact M                    3,810
rows in exact-M-mixed fibres            15,354.           (6)
```

There are no fibres mixing tight and loose rows, but only because every row
in (3) is loose.  This is not a tightness-purity theorem for `H0` or `H1`.

If each face is instead decorated by the ordered reduced endpoints of every
strict-safe component (`H2`), the three-face join is injective on (3).  This
is a finite static reconstruction result, not an arbitrary-height Markov or
continuation theorem.

## 1. Exact maximin census

For a finite speed set, `min_w ||wt||` is piecewise linear.  A positive local
maximum occurs at either a self-cusp or an intersection of two signed linear
branches.  It is therefore enough to test exact rational candidates whose
denominator is

```text
2u,                         u+v,                         |u-v|.  (7)
```

The replay evaluates every candidate with integer residues, reduces the
resulting fractions only for storage, and retains every maximizer.  It finds

```text
rows                                      25,344
tight rows                                     0
loose rows                                 25,344
global minimum                                1/12
number of minimizers                              1.       (8)
```

The unique row and its witnesses are exactly (4)--(5).  This extends
THM-820's `792`-row height-one closure to all `32` height words over every
missing-label set.  It does not exhaust either of THM-820's arbitrary-height
finite boxes.

## 2. Three literal faces and unique gluing

Index the five replacement coordinates by `1,...,5` in increasing label
order.  Use the faces

```text
A=1234,                         B=1235,                         C=1245. (9)
```

For a fixed full missing-label set `R`, a face state consists of its four
ordered literal `(label,height)` coordinates.  The full `R` and the face role
are part of every observation key.  The corresponding partial speed set is

```text
P union {u_i:i in A},
```

and similarly for `B,C`; it has eleven speeds because the fifth base label
is still absent.

The cover (9) has an elementary effective-descent property.

> **Literal gluing lemma.** Three labelled face assignments which agree on
> every ordered overlap glue to one and only one full five-coordinate height
> assignment.

Every coordinate occurs in at least one face, while agreement on overlaps
makes repeated coordinates equal.  This proves existence and uniqueness.
Moreover every unordered pair of replacement coordinates occurs together in
at least one of `A,B,C`: pairs `34`, `35`, and `45` occur respectively in
`A`, `B`, and `C`, and all remaining pairs occur in more than one face.
Consequently the three labelled pair relations together reconstruct the full
labelled pair relation.

For an observation `H`, define its ordered kernel relation on one face by

```text
R_H={(x,y):H(x)=H(y)}.                                  (10)
```

The ordered pair `(x,y)` is never replaced by an unordered line: the two
coordinates must project in the same order to every overlap.  By the literal
gluing lemma, compatible triples of rows from the three relations (10) are in
bijection with ordered pairs of full height words having equal observations
on `A,B,C`.  Thus the join counts below can equivalently be obtained by an
explicit role-specific relation join or by grouping reconstructed full rows;
the replay uses the latter exact implementation.

## 3. `H0`: the labelled live relation

For distinct replacement labels `q,r`, write

```text
lambda=u_q/u_r,                    z=q r^(-1) in F_13^*.
```

The provider-to-owner left-handoff rule from THM-806 is

```text
-1 < z-2lambda-13m <= 1                         for some m in Z, (11)
```

or, with integer centre `k=z-13m`,

```text
(k-1)/2 <= lambda < (k+1)/2.                    (12)
```

`H0` retains the full missing-label context, face role, and labelled support
of (11).  For Tournament Analysis it also computes two deterministic
completions of silent pairs, by increasing and decreasing label, and retains
their score histograms, directed triangles, SCCs, and Hamiltonian-path counts.
Those fingerprints are functions of the labelled relation and do not refine
its fibres.

There are `792*16=12,672` unique states on each face.  Their exact kernel
relations are

| face | support | relation rows | off diagonal | collision fibres | max fibre |
|---|---:|---:|---:|---:|---:|
| `A=1234` | 5,093 | 54,614 | 41,942 | 2,886 | 16 |
| `B=1235` | 5,378 | 51,686 | 39,014 | 2,957 | 16 |
| `C=1245` | 4,776 | 57,582 | 44,910 | 2,873 | 16 |

The `A`--`B` intermediate join has

```text
relation rows                         123,242
off diagonal                           97,898
collision fibres                        5,914
exact-M-mixed fibres                    3,891
rows in those fibres                   16,441
largest fibre                              24.            (13)
```

Adding `C` gives (6).  Because every replacement pair occurs in a face, the
`ABC` partition is exactly the partition of the full five-vertex labelled
live relation; it is not a subtler invariant created by the join.

## 4. `H1`: integer centres add nothing in this bank

`H1` adjoins, in ordered provider-owner position, the centre `k` from (12) on
every live edge.  On the present bank all replacement speeds lie in
`[14,38]`.  A live edge has

```text
1/2 <= lambda <= 38/14 < 3.                              (14)
```

The lower bound below one is THM-806's subunit rule and is automatic when
`lambda>=1`.  Combining (12) and (14) gives

```text
k in {2,3,4,5,6}.                                       (15)
```

Adjacent labels are distinct, so centres congruent to zero or one cannot
occur.  Since `k=z (mod 13)` and the interval (15) contains at most one
representative of any nonzero residue, the labelled edge endpoints already
determine `k`.  Therefore `H1` and `H0` have exactly the same fibres on every
face, their `AB` joins agree, and their `ABC` joins both have `111,006` rows.

This redundancy is bounded-range mathematics.  Once ratios can reach another
centre congruent to the same `z`, the integer centre is again independent
data.

## 5. An explicit quantitative liar fibre

Take

```text
R=(1,2,3,4,5).
```

The two height words

```text
(1,1,1,1,1),             replacements (14,15,16,17,18),
(2,2,2,2,2),             replacements (27,28,29,30,31)  (16)
```

have the identical decorated live relation

```text
(2 -> 1,k=2),             (3 -> 1,k=3),             (4 -> 2,k=2), (17)
```

with every other pair silent.  Nevertheless exact maximin evaluation gives

```text
M(P union {14,15,16,17,18})=1/4,
M(P union {27,28,29,30,31})=12/37.                       (18)
```

Thus `H0=H1` is not an exact-`M` codec.  Both rows are loose, so (18) is not a
mixed tight/loose fibre and supplies no evidence that the live relation is or
is not tightness-pure outside the audited bank.

## 6. `H2`: literal residual-component stalks

For face `F`, put

```text
E_F={t:min_(q in P union {u_i:i in F})||qt||>1/13}.      (19)
```

Each one-speed strict-safe set is the exact interval union

```text
union_(0<=j<u) ((13j+1)/(13u),(13(j+1)-1)/(13u)).        (20)
```

Intersecting (20) over the eleven partial speeds gives all components of
`E_F`.  `H2` retains their ordered reduced rational endpoint word together
with the full `R` and face role.  Its face kernel census is

| face | support | relation rows | off diagonal | collision fibres | max fibre |
|---|---:|---:|---:|---:|---:|
| `A=1234` | 12,668 | 12,680 | 8 | 4 | 2 |
| `B=1235` | 12,669 | 12,678 | 6 | 3 | 2 |
| `C=1245` | 12,670 | 12,676 | 4 | 2 | 2 |

The `AB` join has only four unordered collision pairs:

```text
relation rows=25,352,                off diagonal=8.     (21)
```

The `ABC` join is exactly diagonal:

```text
support=25,344,            relation rows=25,344,
off diagonal=0.                                             (22)
```

For comparison, the endpoint word of the **final** residual set `E_R` is not
injective on speed rows.  Its full-row kernel has

```text
support                                25,330
relation rows                          25,372
off diagonal                               28
unordered collision pairs                 14
largest fibre                               2
fibres mixing exact M                       0.             (23)
```

Three partial component stalks therefore reconstruct more static information
than the final threshold residual alone.  Equation (22) is still only a
finite statement about heights `{1,2}`.  It does not prove that an `H2` state
transports faithfully under a new arbitrary-height comb, gcd descent, scale
change, or another face system.  Such an operation requires the remaining
labelled tooth bank and any arithmetic provenance it uses.

## 7. Memory preflight and computational consequence

The three `H0` face relations contain

```text
54,614+51,686+57,582=163,882 rows.                       (24)
```

Two packed 32-bit face-state identifiers use `1,311,056` bytes.  Charging the
same conservative `56` bytes per row as a packed row plus three aligned
projection indexes uses `9,177,392` bytes.  The full `H0` upper relation uses
`1,776,096` bytes at sixteen bytes per ordered pair.  The three `H2` face
relations contain `38,034` rows and cost `2,129,904` bytes under the same
conservative charge.

The join is therefore easily small enough for a canonical regression.  Its
mathematical value is more limited: `H0` is quantitatively impure, `H1` is
redundant in this band, and literal `H2` is so fine that the three-face codec
is injective.  The next useful compression experiment should lie between
them—for example an ordered component/endpoint-owner incidence word—rather
than interpreting either extreme as the missing arbitrary-height state.

## 8. Tournament Analysis and assumption challenge

For the liar fibre (16), take the five replacement labels as telemetry
vertices.  The pairwise observable is the antisymmetric difference of the two
left-handoff indicators.  Complete the seven silent pairs by increasing label
and then by decreasing label.  The gauges differ by seven edge flips and have

| gauge | score histogram | directed triangles | SCC sizes | Hamiltonian paths | tie Hamiltonian path |
|---|---|---:|---|---:|---|
| increasing | `{0:1,2:2,3:2}` | 2 | `(4,1)` | 5 | `(1,4,2,3,5)` |
| decreasing | `{0:1,1:1,2:1,3:1,4:1}` | 0 | `(1,1,1,1,1)` | 1 | `(5,4,3,2,1)` |

Both speed rows have exactly these fingerprints and the same decorated live
edges, while their values in (18) differ.  The tournament is a checksum of a
pairwise observation, not the theorem predicate.

The vertex challenge is consequently operation-specific.

- Runner vertices retain speed identity but lose simultaneous strict-safe
  component coverage.
- Residue vertices retain the labelled packet but lose height, ratio, and
  metric placement.
- Owner exits and provider teeth are the correct local vertices for collar
  reduction, but their live relation loses the residual interval union.
- Bare safe components retain current threshold geometry but lose the
  remaining labelled tooth bank and descent provenance.
- Kernel-pair rows `(x,y)`, with ordered role projections, are the correct
  vertices for testing whether a proposed quotient is legal.

For exact insertion, the predicate-preserving state is the ordered strict-safe
component union together with the labelled remaining-tooth incidence.  The
action is deterministic interval intersection.  A quotient of that state is
legal only after its continuation fibres, not merely its static fingerprints,
are proved pure.

## Reproduction and certificate digests

```bash
c++ -O3 -std=c++17 \
  04-computation/lrc13_hamming_five_height_two_kernel_pairs_codex_S10.cpp \
  -o /tmp/lrc13_h5_h2_kernel
/tmp/lrc13_h5_h2_kernel
```

The oracle digest hashes, in lexicographic `(R,h)` order, five labels, five
heights, the reduced maximum, the number of maximizers, and every reduced
maximizer as big-endian unsigned 64-bit fields.  The face digest hashes the
context, height word, live mask, ordered centre word, and every reduced
component endpoint in the same format.  They are

```text
oracle
a704bcad9cf023838f17e77d8853b03c9c98a011a92b861960f165ed08f816bd

face states
ccee941334229f071d4b118e8cf1852c39cdc5decccc8ba4c561b0924e31a297
```

For continuity with the independent exploratory implementation, the four
kernel digests concatenate its exact deterministic tuple/Fraction records:

```text
H0 ABC
06186f17284850ca02ca4f032db113991609a5095f5579e8dad88aba13332ca9

H1 ABC
606ea18f5e8768ae38140ceebd9d2ca3acb087b8dff4d1eb102ad6b334767b2d

H2 ABC
9b5b2757f76e596c6ac5be91e2cc00a4c804830cd178d643e3e9bd4e722e36ec

H2 full residual
a750c6c1d8332278f9a71adf790f7f079a0158da117bb17a52148d30b3621dce
```

An optimized build takes about `5.3` CPU seconds (`18.3` seconds wall time on
the shared busy host used for the recorded run).  The stored output SHA-256 is

```text
d08d0702d7a9b1f6c96f4889c68631de10ddfdefc0793ae32a9610ab98f2bf4a.
```

Independent optimized, unoptimized, and AddressSanitizer/UndefinedBehaviorSanitizer
builds produce byte-identical output. ∎
