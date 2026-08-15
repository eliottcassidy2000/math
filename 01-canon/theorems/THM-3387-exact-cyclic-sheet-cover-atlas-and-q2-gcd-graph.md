---
id: THM-3387
title: "Exact cyclic sheet-cover atlas and q=2 gcd graph"
status: >
  PROVED analytic Boolean-fibre iff and q=2 gcd graph + FINITE-EXACT
  literal-body atlas; independent hostile audit pending.  For a degree-q
  split F=qC union U, let B_q(U) be the base points where the transverse
  speeds cover all q sheets.  The safe image is exactly
  (T minus union_c D_c) minus B_q(U), so the divided core is a complement
  tuple iff B_q(U) lies inside its danger union.  On the canonical aligned
  grid this classifies 15,393 of 23,569 body/degree rows, exactly the entire
  THM-3366/S172 body-descended census.  Of these, 15,381 have no global
  transverse full cover and twelve are genuinely rescued by the core.  At
  q=2 two odd speeds form a cover edge iff u+v>7 gcd(u,v); the relation is an
  undirected graph, not a tournament.  These are raw support-row terminals,
  not an additive refined-ledger decrement or a proof of LRC(14).
source: codex-2026-08-14-exact-sheet-cover-atlas
depends_on:
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3366-all-sector-complement-clock-completion
related:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - THM-3381-reflected-residue-affine-phase-transport-and-frozen-tree-stability
script: 04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py
output: 05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out
script_sha256: 9eb77f5293fa41b298b96fbaa29ecaf4899c9dd86054f2f804737297b806e0fb
output_sha256: 455133f5e6c7eaa480875f854142b56aa5315d1fa8055c49da0c8bc0fdc15acd
semantic_sha256: 0a38ff30719647e20b6eecd9e2ea13e700426989ee209791ef62bdeeedacd496
hash_basis: LF-normalized bytes
---

# THM-3387 -- exact sheet cover is a Boolean hypergraph, not a capacity sum

**PROVED analytic Boolean-fibre iff and `q=2` gcd graph + FINITE-EXACT
literal-body atlas; independent hostile audit pending.**

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

Both danger unions are open.  Therefore the divided core tuple `C` covers
all unsupported open cells if and only if `(6)` holds; a failure cannot hide
only on the removed finite grid.  The grid still needs THM-3366's aligned
owner clock, or its extra clock `D` when `k=0`.

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

## 5. Complete literal-body atlas

Enumerate every six-subset `F` of `{1,...,14}` and every `2<=q<=14` having
nonempty core and transverse parts.  This gives `23,569` body/degree rows.
Exact interval equality `(10)` holds in `15,393` and fails in `8,176`.

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

## 6. Exact identification of the S172 census

THM-3366's first body-descended scout S172 uses, for a row `(F,L,D)`,

```text
q=L/D,                    C={v/q:q divides v in F}.      (18)
```

Core danger is always contained in the unsupported target: if `y in D_c`,
all `q` sheets are blocked by speed `qc`.  If any subset of `C` covers the
target, then

```text
core union subset target subset subset-union subset core union, (19)
```

so equality `(10)` follows.  Conversely `(10)` makes `C` a completion.

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

## 7. All-sector structural subcensus

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

## 8. Boundaries and controls

1. **Capacity is far from necessary.**  The exact atlas has `8,973` rows
   beyond the `6,420` strict sum-capacity subclass.
2. **Global transverse survival is not necessary.**  The twelve rows in
   `(16)` require the core-restriction sidecar.
3. **Strictness changes the hyperedges.**  At `q=7,u=1,y=1/2`, the strict
   comb blocks zero sheets, while replacing `<1/14` by `<=1/14` blocks the
   two boundary sheets.
4. **The gcd graph is not a tournament.**  `{1,3}` is a nonedge and survives;
   `{1,9}` is an edge and covers at `y=1/9`.  There is neither an orientation
   nor a lawful rule forcing missing/bidirectional relations into one.
5. **The grid and row key remain typed.**  The theorem uses `D=L/q` and the
   exact `(body,D)` key.  Body-only, divisor-only, and nearby quotient keys do
   not inherit `(10)`.
6. **S172 is not the universal pool census.**  S173/S174 allow clocks not
   descended from the body and close more raw rows; the present theorem
   classifies exactly the body-descended subengine.

No physical tail realization, arbitrary reflected-phase transport, new
refined-ledger subtraction, tournament reconstruction, or LRC(14) follows.

## 9. Exact verification

The standard-library companion pins THM-3385 and the independently audited
THM-3366 artifacts, then:

- checks `(13)` on all `4,950` odd pairs below `200` by an independent integer
  event sweep;
- checks all `23,569` literal body/degree rows by exact cell-interval equality;
- independently separates `15,381` global transverse survivors from the
  twelve core rescues;
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
