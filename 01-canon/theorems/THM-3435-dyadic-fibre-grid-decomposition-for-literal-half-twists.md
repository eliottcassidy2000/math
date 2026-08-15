---
id: THM-3435
title: "Dyadic fibre-grid decomposition for literal half twists"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  AUDIT REQUIRED.  Literal half-twist masks on Q=2^a R split exactly into
  repeated dyadic grids over the odd base R.  Base projection preserves the
  support and fibre count but loses a selected coset/orientation sidecar.
  The even cap-seven census through Q=362 is FINITE-EXACT only; Q=366 and the
  all-Q even classification remain open.  No LRC(14) consequence is claimed.
source: root even-half-rank7 dyadic-fibre session, 2026-08-15
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3434-seventeen-fibre-two-sided-mass-closure
related:
  - THM-3432-order-two-fixed-half-parity-transplant
  - THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures
script: 04-computation/lrc_dyadic_fibre_grid_decomposition_thm3435.py
output: 05-knowledge/results/lrc_dyadic_fibre_grid_decomposition_thm3435.out
script_sha256: f4f613e95fea20a7cba98c5c13dd99c7eae2140256cd43e392a7e8f1d829c804
output_sha256: b713d888b3832e80b54b5b0389f53b06d41f9ed0f6c3fdcd17c54809f35fa114
semantic_sha256: c8033efed04555c46a3414c9dd99121a4ffc066b7e1c6489597c69d57922220c
finite_script: 04-computation/lrc_dyadic_fibre_grid_census_thm3435.py
finite_output: 05-knowledge/results/lrc_dyadic_fibre_grid_census_thm3435.out
finite_script_sha256: d301f0cc802b311e578109a571e85e1ecf1e81a20efdd688bc4c91b8dc6c530d
finite_output_sha256: df94ab3b8b085fb9ddaf668438d6b07af92c0095a3dfee5774c2c85e13fe5b0d
finite_semantic_sha256: 35dfbba192b4b17b72cda5a019692a9af5ee1ce406f3ca4a7d7ff551fe349096
hash_basis: LF-normalized bytes
---

# THM-3435 -- dyadic fibre-grid decomposition for literal half twists

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
AUDIT REQUIRED.**

## 1. Exact all-modulus statement

For `Q>=2` and a residue `r modulo 2Q`, write

```text
B_(Q,r)={ell in Z/QZ: ||r(2ell+1)/(2Q)||<1/14}.          (1)
```

The transverse owner universe imposes `Q does not divide r`, excluding the
universal zero residue and the empty order-two residue.  The identity below
is valid even without that exclusion.

Factor

```text
Q=2^a R,       R odd,
b=max{c<=a:2^c divides r},       r=2^b s,
d=2^(a-b).                                                (2)
```

The truncated valuation `b`, and hence `d`, is intrinsic to the residue
class modulo `2Q`.  Identify a sheet uniquely as

```text
ell=j+Rt,       j in Z/RZ,       t in Z/2^a Z.           (3)
```

On the fibre over `j`, define the dyadic grid

```text
K_(r,j)={t in Z/2^a Z:
 ||s(2j+1)/(2dR)+st/d||<1/14}.                           (4)
```

Then the literal block has the exact disjoint fibre decomposition

```text
B_(Q,r) intersect {j+Rt:t in Z/2^a Z}
 ={j+Rt:t in K_(r,j)}.                                  (5)
```

The predicate in `(4)` depends only on `t modulo d`.  Every selected
`d`-grid point is therefore repeated exactly `2^b` times in the fibre.  If
`d<7` -- equivalently `d in {1,2,4}` -- the grid arcs are disjoint and

```text
K_(r,j) is nonempty
 iff ||s(2j+1)/(2R)||<d/14;                              (6)
```

when nonempty it is one residue class modulo `d`, repeated `2^b` times.
All inequalities remain strict.

There is also an exact count at every active depth.  Suppose `b<a`, and
write the power of two

```text
d=7q+c,             c in {1,2,4}.                        (7)
```

Define the smaller-radius odd-base mask bit

```text
epsilon_j=1 iff
 ||(s+qR)(2j+1)/(2R)||<c/14.                            (7a)
```

After the affine coordinate change `n=s(t modulo d)` on `Z/dZ`, the
selected grid points form one cyclic interval, and

```text
#{n in Z/dZ:selected over j}=q+epsilon_j,
|K_(r,j)|=2^b(q+epsilon_j).                             (7b)
```

In particular, `(7b)` specializes to `(6)` when `d=2,4`; when `d>=8`,
every base fibre is nonempty and its reduced count is `floor(d/7)` or
`ceil(d/7)`.  No active branch `b<a` can meet a strict radius-`1/14`
endpoint.  Inactive pullbacks `b=a` can meet one, so their strict exclusion
must still be retained.

When `b=a`, so `d=1`, `(5)--(6)` reduce to the literal pullback

```text
B_(Q,2^a s)=pi^(-1)B_(R,s),       pi(ell)=ell mod R.     (8)
```

In this deepest branch, transversality is preserved exactly:

```text
Q does not divide 2^a s  iff  R does not divide s.       (9)
```

## 2. The two live low-depth charts

The rank-seven novelty lies only at `a=1,2`:

- On `Q=2R`, an even coefficient `r=2s` pulls an ordinary radius-`1/14`
  block on `R` back to both fibre sheets.  An odd coefficient selects one
  oriented sheet over the widened support

  ```text
  W_(R,r)={j:||r(2j+1)/(2R)||<1/7}.                     (10)
  ```

- On `Q=4R`, coefficients with truncated valuations `b=2,1,0` respectively
  select all four sheets over a radius-`1/14` support, one two-sheet coset
  over a radius-`1/7` support, or one oriented sheet over a radius-`2/7`
  support.

For `Q=2R` and odd `r` other than the empty coefficient `R`, let
`tau(ell)=ell+R`.  Equations `(4)--(6)` sharpen to

```text
tau B_(2R,r)=B_(2R,2R-r),
B_(2R,r) intersect B_(2R,2R-r)=empty,
B_(2R,r) union B_(2R,2R-r)=pi^(-1)W_(R,r).             (11)
```

Thus a complementary odd pair costs two owners and fills both sheets over
one widened block.  An unpaired odd owner still carries a genuine orientation
bit; replacing it by its projected support is not a literal-cover operation.

The same grid calculation has a sharp Boolean branch boundary.  For a
general radius `u`, write

```text
H_Q^(u)(r)={ell:||r(2ell+1)/(2Q)||<u}.
```

For `k in {1,2,3}`, assume `2^(k-1)|Q`, choose any seed `v`, and take the
complete coefficient fibre

```text
r_h=v+hQ/2^(k-1) modulo 2Q,       0<=h<2^k.            (11a)
```

For `k=1,2`, the `2^k` literal masks `B_(Q,r_h)` are pairwise disjoint and

```text
union_h B_(Q,r_h)=H_Q^(2^k/14)(2^k v).                 (11b)
```

For `k=3`, the eight masks cover every sheet with multiplicity one or two.
Two branch labels of the same parity never overlap; on a fixed sheet, a
possible double hit is an adjacent pair in a local cyclic order.  That order
is multiplied by the odd unit `2ell+1`, so it is sheet-dependent rather than
a global orientation.  The eight coefficients need not all be transverse,
and dropping one is not uniformly a cover, so this identity alone is not a
rank-seven certificate.

If `a>=3`, then `8|Q`.  The rank-four half atom on `8` from THM-3416 pulls
back to `Q`, so

```text
rho_H(Q)<=4.                                            (12)
```

Consequently `a=1,2` are the only even dyadic depths at which new cap-seven
support can occur.

## 3. Proof

Substitute `(2)--(3)` into `(1)`:

```text
r(2ell+1)/(2Q)
 =s(2j+1)/(2dR)+st/d.                                  (13)
```

This is exactly `(4)--(5)`.  If `b<a`, then `s` is odd, so multiplication by
`s` permutes `Z/dZ`; if `b=a`, then `d=1`.  Hence the fibre is a complete
`d`-grid with every point repeated `2^b` times.

For `d<7`, the radius-`1/14` arcs around consecutive grid points are
disjoint.  Multiplication by `d` maps their union bijectively onto the strict
arc of radius `d/14`, proving `(6)` and the one-coset statement.

For the exact active-depth count, change coordinates to `n=st modulo d` and
put `u=s(2j+1)/(2R)`.  The selected integers modulo `d` are exactly the
cyclic interval

```text
distance(n+u,dZ)<d/14=q/2+c/14.                        (13a)
```

If `q` is even, this interval contains `q` points plus one exactly when
`||u||<c/14`.  If `q` is odd, it contains `q` points plus one exactly when
`||u+1/2||<c/14`.  But

```text
(s+qR)(2j+1)/(2R) = u+q(2j+1)/2 = u+q/2 modulo 1,
```

so the extra bit is precisely `(7a)`.  Restoring the `2^b` repeated copies
proves `(7b)`, including the dense-grid floor/ceiling consequence.

An active strict endpoint would give, for some integer `m`,

```text
14r(2ell+1)=2Q(14m plus-or-minus 1).
```

Its two sides have respective 2-adic valuations `b+1` and `a+1`, impossible
when `b<a`.  This proves the active endpoint claim; the inactive hostile in
the exact companion shows why it cannot be extended to `b=a`.  Equations
`(8)--(9)` are the `d=1` specialization.

For `(11)`, translation by `R` adds `r/2=1/2 modulo 1` to the phase.  On the
other hand, replacing `r` by `2R-r` changes the phase to `1/2` minus the old
phase.  These operations agree after taking distance to the nearest integer.
The two strict radius-`1/14` arcs are disjoint, and `(6)` with `d=2` identifies
their fibre support with `(10)`.  This proves `(11)`.

For completeness, fix a sheet in `(11a)`.  Relative to the seed phase, the
branch indexed by `h` is shifted by

```text
h(2ell+1)/2^k modulo 1.
```

The odd multiplier permutes the `2^k` equally spaced branch centres.  For
`k=1,2`, their open radius-`1/14` arcs are pairwise disjoint, and multiplication
of the phase by `2^k` fuses their union into the radius-`2^k/14` arc in
`(11b)`.  For `k=3`, spacing `1/8` is smaller than the arc diameter `1/7`, so
the eight arcs cover the circle; the diameter is smaller than `1/4`, so no
three meet.  Only adjacent labels in the sheet-local order can meet, and
odd multiplication preserves the alternating parity classes.  This proves
the Boolean branch claims. **QED.**

## 4. The naive odd-core descent is false

It is tempting to identify the odd-owner projection in `(10)` with the
ordinary half block on `R`.  That loses a factor of two in the radius.  The
smallest decisive hostile is

```text
Q=14: (1,3,4,5,9,11,13),                               (14)
```

an exact seven-block partition.  At `R=7`, however, the union of **every**
transverse ordinary half mask reaches only the reflection-fixed sheet `3`,
so no cover exists at any rank.  Therefore

```text
rho_H(2R)=rho_H(R)                                      (15)
```

is false.  The correct target is the mixed-radius grid carrier `(4)--(6)`.

This also separates the present map from its two closest proved mechanisms.
THM-3432 restricts a **fixed-zero** block on `2R` to one parity class of
sheets and obtains an ordinary half block on `R`; it retains an order-two
owner as a conditional sidecar.  Here we project a literal half block on all
sheets along the deck fibre.  An active odd owner then has the widened
radius `1/7` support and loses its orientation bit, so the THM-3432 transfer
cannot be used as an ownerwise inverse.  THM-3434 classifies ordinary
radius-`1/14` covers on odd moduli.  It therefore controls the deepest
pullback branch `(8)` and the odd proper joint periods used below, but it does
not classify the active widened-radius branches.

## 5. FINITE-EXACT even cap-seven census

This section is deliberately not an all-`Q` theorem.  Let

```text
A_0={8,9,10,11,12,13,15,23,25,29,51}.                  (16)
```

THM-3416 and audited THM-3434 prove that these are already-supported literal
half bases through rank seven in the even pullback and odd sectors.  In the
finite universe

```text
2<=Q<=362,       Q even,       no member of A_0 divides Q,          (17)
```

an ordinary transverse literal half cover by at most seven blocks exists
exactly when

```text
some member of {14,38,68,148} divides Q.               (18)
```

The four frozen positive witnesses are

```text
14:  (1,3,4,5,9,11,13),
38:  (1,9,17,20,21,29,37),
68:  (8,11,23,24,45,56,57),
148: (8,33,41,100,107,115,140).                        (19)
```

The `14` and `68` rows have odd multiplicity everywhere (the first is a
partition); the `38` and `148` rows have respectively four and eight
double-covered sheets.  These controls show that the dyadic grid mechanism
is compatible with both partition and overlap boundaries.

### Why the finite negative search is complete in `(17)`

For a family `C`, put

```text
m_Q(r)=Q/gcd(Q,r),       L=lcm_(r in C)m_Q(r),
h=Q/L=gcd(Q,r:r in C).                                 (20)
```

Then `h` divides every owner and direct cancellation gives

```text
B_(Q,hs)=pi_L^(-1)B_(L,s).                             (21)
```

Thus every cover descends to a joint-period cover on a divisor `L|Q` with
the same owner count.  The companion processes `(17)` in increasing order.
If a discovered atom divides `Q`, `(19)` pulls back and no search is needed.
Otherwise it searches the complete joint-period bank, adjoining one breaker
bit for every prime divisor of `Q`; all bits are covered exactly when
`gcd(Q,C)=1`, equivalently `L=Q`.  A proper even `L` was already decided by
induction, while an odd `L` is excluded by THM-3434 because `(17)` removes all
of its positive bases.  This proves the finite consequence `(18)` rather
than merely scanning full-period rows.

The census searches `60` new-period rows, of which `56` are negative, and
visits `1,799,207` memoized states and `1,834,885` branches.  The first even
modulus beyond the frozen range which avoids both `(16)` and `(18)` is

```text
Q=366.                                                   (22)
```

It is the explicit stopping boundary: this theorem neither searches nor
classifies it.  No inference beyond `Q=362` is licensed.

## 6. Connection and loss ledger

| field | exact content |
|---|---|
| source | a labelled literal half block on `Q=2^a R` |
| target | its odd-base smaller-radius bit and exact repeated dyadic-grid count in each fibre |
| map | `ell=j+Rt -> j`, with grid coordinate `t` |
| preserved with full sidecar | labels, strict endpoints, every sheet incidence, fibre counts, quotient two-part |
| destroyed by count projection | affine location of the cyclic interval, selected coset/orientation, layer overlaps and literal OR coverage |
| restoring sidecar | the subset `K_(r,j)` in `(4)`; equivalently its affine cyclic interval, and for `d<7` its unique coset modulo `d` |
| cheapest affine hostile | at `Q=8`, residues `1` and `3` have the same fibre count but different masks |
| cheapest descent hostile | the `Q=14` partition versus the uncoverable ordinary `Q=7` layer |
| finite boundary | exact only through `Q=362`; `Q=366` is unclassified |

There is no intrinsic pairwise orientation between owners, so a tournament
would be cosmetic.  The faithful object is a labelled fibre-grid clutter.

## 7. Exact companion and scope

Run

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_dyadic_fibre_grid_decomposition_thm3435.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_dyadic_fibre_grid_decomposition_thm3435.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_dyadic_fibre_grid_census_thm3435.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_dyadic_fibre_grid_census_thm3435.py
```

The first standard-library companion compares `(4)--(7b)` directly with
literal masks for every residue through `Q=192`: `4,755,520` sheet cells and
`449,888` active fibre rows.  It audits the active endpoint obstruction, the
inactive `Q=14,r=2` endpoint hostile, all three Boolean branch depths, and
the `Q=8` affine-loss hostile.  The second companion checks the earlier
fibre formulation through `Q=160` (`1,812,487` fibre rows), checks `5,700`
complementary two-sheet rows, replays `(14)` and `(19)`, and performs exactly
the bounded inductive census `(17)`.  Normal and optimized outputs are
byte-identical for both companions.

The all-modulus content is only the elementary fibre identity and its stated
corollaries.  The census is **FINITE-EXACT**, not evidence promoted to an
all-`Q` even rank-seven classification.  Nothing here supplies a physical
time, arbitrary common centre, LRC row, or decrement.  `LRC(14)` remains
open.
