---
id: THM-3387
title: "Exact cyclic sheet-cover atlas and q=2 gcd graph"
status: >
  PROVED analytic pointwise and aligned-grid criteria + PROVED q=2 gcd graph
  + FINITE-EXACT literal-body atlas + INDEPENDENTLY HOSTILE-AUDITED after the
  MISTAKE-382 endpoint repair.  Pointwise safe-image equality is equivalent
  to B_q(U) lying in the core danger union A_C.  Universally, core coverage of
  the unsupported open cells is equivalent to B_q(U) minus A_C lying on the
  removed D-grid; these criteria can differ at strict endpoint handoffs.  No
  such exception occurs for F subset {1,...,14}, where the theorem classifies
  15,393 of 23,569 body/degree rows, exactly the THM-3366/S172 body-descended
  census.  At q=2 two odd speeds form a cover edge iff
  u+v>7 gcd(u,v); the relation is an undirected graph, not a tournament.
  These are raw support-row terminals, not an additive refined-ledger
  decrement or a proof of LRC(14).
source: codex-2026-08-14-exact-sheet-cover-atlas
audit: independent pointwise proof, event-bitset and interval atlas, endpoint hostile, exact S172 key replay, q2 graph, and harmonic-lattice audit
depends_on:
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3366-all-sector-complement-clock-completion
related:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - THM-3381-reflected-residue-affine-phase-transport-and-frozen-tree-stability
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
script: 04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py
output: 05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out
script_sha256: 9b0b46874a569d674b937b37cf74a8985fca2b77e3e480a75fb4924ea602f25a
output_sha256: b4d9ce439bab4501bfd5e2cf13eb0b0e3685b7364f30e43b7d5ca9138d25cb5c
semantic_sha256: de40f9da0f3b335d0de52bad0c75f586c15508cdf4d94b082d2f69c1098e5bc8
hash_basis: LF-normalized bytes
---

# THM-3387 -- exact sheet cover is a Boolean hypergraph, not a capacity sum

**PROVED analytic pointwise and aligned-grid criteria + PROVED `q=2` gcd
graph + FINITE-EXACT literal-body atlas + INDEPENDENTLY HOSTILE-AUDITED after
the MISTAKE-382 endpoint repair.**

## 1. From a sufficient bound to the true object

Use THM-3385's notation at radius `1/14`:

```text
F={q c:c in C} disjoint union U,             pi_q(x)=qx, (1)
G_F=T minus union_(v in F)D_v.
```

For `y in T`, label the fibre by `Z/qZ` and define the transverse blocked-
sheet hyperedge

```text
E_u(y)={k in Z/qZ: (y+k)/q in D_u}.                     (2)
```

The full transverse-cover locus is

```text
B_q(U)={y: union_(u in U) E_u(y)=Z/qZ}.                 (3)
```

THM-3385 bounds `|E_u(y)|` and uses a union bound to prove `B_q(U)=empty`
under its capacity hypothesis.  The set `(3)`, rather than that sufficient
bound, is the exact object.

## 2. Exact image and iff

For every base point `y`, each core speed has sheet-independent status:

```text
qc(y+k)/q=cy+ck.                                       (4)
```

Thus a core-dangerous `y` has no safe preimage.  At a core-safe `y`, a safe
preimage exists exactly when the transverse hyperedges do not cover the
whole fibre.  Pointwise,

```text
pi_q(G_F)
 =[T minus union_(c in C)D_c] minus B_q(U).             (5)
```

Consequently

```text
pi_q(G_F)=T minus union_(c in C)D_c
 iff B_q(U) subset union_(c in C)D_c.                   (6)
```

This is the promised exact Boolean-fibre criterion.  THM-3385's capacity
sum is one sufficient way to make the left set in `(6)` empty; it is not an
equivalent reformulation.

## 3. Canonical grid form

Let

```text
L=14 lcm(F),                         D=L/q,             (7)
```

and let `S_D` be the THM-2928 projected safe-cell support.  All event
boundaries lie on the `D`-grid.  For a core clock `c`, this is THM-3385's

```text
D/(14c)=lcm(F)/(qc) in N.                              (8)
```

For a transverse speed `u`, a boundary on sheet `k` has the form

```text
y=q(n+/-1/14)/u-k,
D y=L n/u +/- L/(14u)-kD in Z.                         (9)
```

Hence every Boolean word in `(2)`--`(3)` is constant on each open `D`-cell,
and

```text
U_D(S_D)
 =[union_(c in C)D_c union B_q(U)] minus (D^(-1)Z/Z).  (10)
```

Put `A_C=union_(c in C)D_c` and `Gamma_D=D^(-1)Z/Z`.  Equation `(10)` gives
the universal, grid-correct criterion

```text
C covers U_D(S_D) iff B_q(U) minus A_C subset Gamma_D.  (10a)
```

This is weaker than pointwise image equality `(6)`: a transverse full cover
can be core-safe only at removed grid points.  MISTAKE-382 freezes the clean
six-speed witness

```text
F=(4,5,8,9,10,18), q=2, C=(2,4,5,9), U=(5,9),
(L,D)=(5040,2520),    B_2(U) minus A_C={3/14,11/14}.     (10b)
```

For the literal atlas `F subset {1,...,14}`, no endpoint-only exception can
occur.  Every core clock satisfies `c<=7`.  If an open core union covered both
punctured sides of a core-safe point, distinct clocks `c,d` would meet there
on opposite strict boundaries.  Eliminating the point from
`cy=+1/14 (mod 1)` and `dy=-1/14 (mod 1)` forces `14 | c+d`, impossible for
distinct `c,d<=7`.  Hence `(6)` and `(10a)` agree throughout that atlas.  The
grid still needs THM-3366's aligned owner clock, or its extra clock `D` when
`k=0`.

## 4. The exact `q=2` gcd graph

When `q=2`, every transverse speed is odd and blocks at most one sheet.  Let
`u,v` be distinct odd speeds.  Put the first sheet at `t`; the other is
`t+1/2`.  For `u` to block the first and `v` the second, `t` must lie in open
intervals centred at

```text
a/u                  and                  (2b+1)/(2v), (11)
```

with radii `1/(14u)` and `1/(14v)`.  If `g=gcd(u,v)`, the minimum circular
distance between the two centre grids is

```text
g/(2uv).                                                (12)
```

Indeed the numerator `2va-u(2b+1)` runs through the residue class `g mod 2g`,
and its least circular absolute value is `g`.  The intervals overlap exactly
when

```text
g/(2uv)<1/(14u)+1/(14v),
equivalently                    u+v>7gcd(u,v).           (13)
```

Define an undirected edge by `(13)`.  A set of odd transverse speeds covers
both sheets somewhere if and only if it contains an edge: any full cover
chooses one blocker on each sheet, and an edge supplies such a pair.
Therefore

```text
B_2(U)=empty iff U is independent in the gcd graph.     (14)
```

On the literal odd pool `{1,3,5,7,9,11,13}`, the nonempty independent sets
of size at most six are exactly

```text
{1},{3},{5},{7},{9},{11},{13},{1,3},{1,5},{3,9}.        (15)
```

This yields `147` one-odd/five-even and `3*binom(7,4)=105` two-odd/four-even
rows, for `252` exact degree-two rows.  The last `105` lie beyond THM-3385's
strict capacity-sum subclass.

The relation in `(13)` is intrinsic and binary, but symmetric.  Orienting it
would be a gauge choice, not a tournament observable.  For `q>=3`, the native
objects are the hyperedges `(2)`, including singleton, antipodal-pair and
parity-triple edges at `q=4,6`; missing and simultaneous edges must remain
literal.

## 5. The hidden ancestry is a `(3,5)` lattice, not a tree

There is an all-speed simplification of the nonedge relation.  Divide distinct
odd `u,v` by their gcd.  The reduced values are coprime, odd, and a nonedge
exactly when their sum is at most seven.  The only possibilities are
`{1,3}` and `{1,5}`.  Therefore the undirected **noncover graph** on all odd
positive speeds has precisely the edges

```text
{a,3a} and {a,5a},                         a odd.        (15a)
```

Strip all factors of three and five from an odd speed.  Every connected
component has a unique root `r` with `gcd(r,15)=1` and vertex set

```text
V_r={r 3^i 5^j:i,j>=0}.                                (15b)
```

It is the quadrant lattice with horizontal `x3` and vertical `x5` edges.
The square `r,3r,15r,5r` shows why the ancestry object is not a tree.

The natural speed labels give each component exact harmonic mass

```text
sum_(n in V_r) 1/n
 =(1/r)(sum_i 3^(-i))(sum_j 5^(-j))=15/(8r).           (15c)
```

More generally, for any subset `R` of odd roots coprime to fifteen, positivity
and disjointness give, with infinity allowed,

```text
sum_(n in union_(r in R)V_r) 1/n=(15/8)sum_(r in R)1/r. (15d)
```

At lattice depth `d=i+j`, the free binary ancestry has `2^d` words but only
`d+1` distinct integer values.  The value `r3^i5^(d-i)` has collision
multiplicity `binom(d,i)`.  The distinct-support and multiplicity-weighted
harmonic shell masses are respectively

```text
H_d=(1/(2r))(5/3^d-3/5^d),        15H_(d+2)=8H_(d+1)-H_d,
W_d=(1/r)(8/15)^d,                15W_(d+1)=8W_d.       (15e)
```

Thus commutation changes the recurrence as well as the count: ancestry words,
integer support, collision multiplicity, and harmonic weight are four
different profiles.  This is not the Berggren/Fibonacci tree: its branch
words retain order, whereas multiplication by three and five has already
abelianized order.  Any transport between them needs that lost word-order
sidecar.

Thus every root subset has an exact harmonic-series realization under this
ancestry decomposition.  The integer embedding is load-bearing, as in
THM-3382.  This does **not** mean that an entire component is one safe
transverse set: safe sets are cliques of the noncover graph, and the lattice
has no triangles.  Only single vertices and single edges survive globally.

## 6. Complete literal-body atlas

Enumerate every six-subset `F` of `{1,...,14}` and every `2<=q<=14` having
nonempty core and transverse parts.  This gives `23,569` body/degree rows.
The core-only cell equality

```text
U_D(S_D)=A_C minus Gamma_D
```

holds in `15,393` and fails in `8,176`.  By the preceding no-handoff lemma,
the same `15,393` rows satisfy the pointwise criterion `(6)`.

| `q` | candidate-exact rows |
|---:|---:|
| 2 | 252 |
| 3 | 588 |
| 4 | 619 |
| 5 | 1,619 |
| 6 | 1,478 |
| 7 | 2,079 |
| 8 | 1,152 |
| 9 | 1,205 |
| 10 | 1,269 |
| 11 | 1,287 |
| 12 | 1,271 |
| 13 | 1,287 |
| 14 | 1,287 |

Among the `15,393`, exactly

```text
6,420   satisfy THM-3385's strict capacity sum,
15,381  have B_q(U)=empty by exact event sweep,
12      have B_q(U) nonempty but contained in the core danger union. (16)
```

The twelve core-rescued rows occur at degrees `3,5,6`.  Their complete frozen
list is in the exact output.  A smallest one is

```text
F=(1,3,6,8,11,13),       q=6,       C=(1),
U=(1,3,8,11,13),          (L,D)=(48048,8008).            (17)
```

Here the transverse speeds do cover all six sheets on some open events, but
every such event lies in `D_1`; the core deletes the apparent obstruction.
This is the minimal reason the global-transverse count is twelve below the
body-relative atlas.

## 7. Exact identification of the S172 census

THM-3366's first body-descended scout S172 uses, for a row `(F,L,D)`,

```text
q=L/D,                    C={v/q:q divides v in F}.      (18)
```

Core danger cells are always contained in the unsupported target: if
`y in D_c`, all `q` sheets are blocked by speed `qc`.  If any subset `C'` of
`C` covers the target, then the grid-correct chain is

```text
A_C minus Gamma_D subset U_D(S_D)
 subset A_(C') minus Gamma_D subset A_C minus Gamma_D.   (19)
```

Thus core-only cell equality follows.  Conversely that equality makes `C` a
completion.  This argument needs no reconstruction of deleted endpoints.

Every atlas row has `|C|<=5`, hence injects into S172's closed row set.  The
atlas has `15,393` rows, exactly S172's independently audited total; therefore
the injection is equality.  The exact companion additionally proves every
core tuple irredundant and reproduces S172's least-completion profiles:

| completion size | rows | occurrences |
|---:|---:|---:|
| 1 | 12,658 | 3,427,329,787,389 |
| 2 | 2,029 | 182,339,009,252 |
| 3 | 409 | 35,772,924,216 |
| 4 | 150 | 1,643,011,511 |
| 5 | 147 | 12,170,729,897 |

Thus the former black-box `15,393 / 3,659,255,462,265` S172 result is exactly
the cyclic Boolean-sheet atlas, not an accidental set-cover census.

## 8. All-sector structural subcensus

Applying THM-2928's inherited support cutoffs and THM-3366's clock budgets to
the same atlas gives:

| `k` | structural rows | denominator-shape occurrences |
|---:|---:|---:|
| 0 | 15,246 | 93,732,978,513,930 |
| 1 | 15,393 | 3,659,255,462,265 |
| 2 | 15,393 | 133,947,094,813 |
| 3 | 15,393 | 4,445,930,697 |
| 4 | 2,735 | 15,952,450 |
| 5 | 1,104 | 154,344 |
| 6 | 291 | 291 |

For `k=0`, the omitted `147` rows are exactly the five-core degree-two family
which needs five pool clocks but has budget four.  The table is contained in
THM-3366's existing pool-14 terminal census and is not additive.  In
particular, only THM-3366's already intersected seven rows and `7,648`
occurrences transfer to the current refined `k=3` ledger.

## 9. Boundaries and controls

1. **Capacity is far from necessary.**  The exact atlas has `8,973` rows
   beyond the `6,420` strict sum-capacity subclass.
2. **Global transverse survival is not necessary.**  The twelve rows in
   `(16)` require the core-restriction sidecar.
3. **Deleted endpoints are a genuine sidecar.**  The witness `(10b)` has
   exact cell completion but not pointwise image equality.  It lies outside
   the literal atlas; the no-handoff argument above is load-bearing there.
4. **Strictness changes the hyperedges.**  At `q=7,u=1,y=1/2`, the strict
   comb blocks zero sheets, while replacing `<1/14` by `<=1/14` blocks the
   two boundary sheets.
5. **The gcd graph is not a tournament.**  `{1,3}` is a nonedge and survives;
   `{1,9}` is an edge and covers at `y=1/9`.  There is neither an orientation
   nor a lawful rule forcing missing/bidirectional relations into one.
6. **The grid and row key remain typed.**  The theorem uses `D=L/q` and the
   exact `(body,D)` key.  Body-only, divisor-only, and nearby quotient keys do
   not inherit `(10)`.
7. **S172 is not the universal pool census.**  S173/S174 allow clocks not
   descended from the body and close more raw rows; the present theorem
   classifies exactly the body-descended subengine.

No physical tail realization, arbitrary reflected-phase transport, new
refined-ledger subtraction, tournament reconstruction, or LRC(14) follows.

## 10. Exact verification

The standard-library companion pins THM-3385 and the independently audited
THM-3366 artifacts, then:

- checks `(13)` on all `4,950` odd pairs below `200` by an independent integer
  event sweep;
- checks `4,445` finite `(3,5)`-component edges and exact rectangular
  harmonic factorizations;
- checks all `23,569` literal body/degree rows by exact cell-interval equality
  and verifies there is no pointwise/grid discrepancy;
- independently separates `15,381` global transverse survivors from the
  twelve core rescues;
- freezes `(10b)` as a positive endpoint-only hostile outside the atlas;
- reproduces every degree, sector, occurrence, and S172 size profile using a
  downward divisor recurrence rather than Mobius summation; and
- checks strict versus closed endpoints and the capacity-equality hostiles.

It contains no floating literal or optimization-dependent `assert`.
Reproduce with

```text
python 04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py
python -O 04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py
```

Ordinary and optimized runs LF-normalized-byte-match the stored output.

**QED.**
