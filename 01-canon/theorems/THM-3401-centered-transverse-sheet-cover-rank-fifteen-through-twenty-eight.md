---
id: THM-3401
title: "Centered transverse sheet-cover rank from fifteen through twenty-eight"
status: >
  PROVED exact rank formula for every q=15,...,28 plus sharp q=14 and q=29
  boundary theorems; VERIFIED-EXACT by a self-contained literal classifier,
  construction audit, and independent finite set-cover search;
  INDEPENDENTLY AUDITED.  This is a fixed-source-centre-zero theorem, a
  proper sub-slice of the zero-cochain locus after MISTAKE-384, and gives no
  LRC(14) decrement.
source: root-2608-centered-sheet-rank-2026-08-15
audit: self-contained arithmetic proof; 328 general quotient-block checks for q=14,...,29; all 287 transverse residue owners for q=15,...,28; 122 construction/irredundancy checks; 9828 literal set-cover subsets; q14 strict noncover, 24 q28 endpoint incidences, and q29 C14 edge-cover boundary; independent proof/type/strict-endpoint/cycle/replay/hash-seed audit complete
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
related:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3395-small-sheet-typed-cover-star-cochain
artifacts:
  - 04-computation/lrc15_first_effective_triphase_mode_probe_20260814.py
  - 05-knowledge/results/lrc15_first_effective_triphase_mode_probe_20260814.out
script: 04-computation/lrc_centered_transverse_sheet_cover_rank_thm3401.py
output: 05-knowledge/results/lrc_centered_transverse_sheet_cover_rank_thm3401.out
script_sha256: fff146868de4b1ec5993ed404fca4200e8e3eac47a7cb902ff51d559eef228e0
output_sha256: 12f2cf337c982068c8cfa0cb351a2772f05e55c1be18df4d0414f9a7251327dd
semantic_sha256: 41c31d22b1a1b6fa3637d129fad3a311565d21bf45ad9f322bd550d595c157cd
hash_basis: LF-normalized bytes
---

# THM-3401 -- centered transverse sheet-cover rank from fifteen through twenty-eight

**PROVED exact rank formula for every `q=15,...,28`, with sharp lower and
upper boundary theorems, VERIFIED-EXACT independent literal set-cover search,
and an independent audit.**

## 1. Statement

For `q>=2`, a positive speed `u` is **transverse** when `q` does not divide
`u`.  At common source centre zero put

```text
D_(q,u)(0)={ell in Z/qZ: ||u ell/q||<1/14}.           (1)
```

For every `15<=q<=28`, the minimum cardinality of a finite transverse speed
set `U` satisfying

```text
union_(u in U) D_(q,u)(0)=Z/qZ                       (2)
```

is

```text
r_0(q)=phi(q)/2 + #{p prime:p|q and p<q}
      =phi(q)/2+omega(q)-1_(q is prime).              (3)
```

The exact values are

```text
q:       15 16 17 18 19 20 21 22 23 24 25 26 27 28
r_0(q):   6  5  8  5  9  6  8  7 11  6 11  8 10  8. (4)
```

**Gauge correction (MISTAKE-384).**  Fixed source centre zero implies that
all THM-3398 affine gaps vanish, but the converse only puts all centre lifts
at one arbitrary common rational `c`.  It does not force `c=0`.  Translating
that centre to zero induces a common cyclic sheet relabelling only when `qc`
is integral.  Thus this theorem computes the fixed-zero rank `r_0`, not the
rank on the entire zero-cochain locus.  Its statement and proof below were
always fixed at zero and are unchanged by this typing repair.

The interval is sharp in both directions.  At `q=14` no transverse
common-centre cover exists.  At `q=29` a transverse common-centre cover does
exist with exact minimum rank seven; the continuation of `(3)` would instead
give fourteen.

## 2. Inheritance and connection contract

THM-3398 characterizes arbitrary cyclic sheet covers by selected consecutive
phase blocks and a complete affine cochain.  Restricting to centre zero kills
every affine gap, but it does not make all owners alike: the remaining
load-bearing coordinate is the quotient order

```text
g=gcd(u,q),              m=q/g.                       (5)
```

| field | exact connection |
|---|---|
| source | THM-3398's all-`q` selected-block cover at a common source time |
| target | an arithmetic set-cover rank on `Z/qZ` at time zero |
| map | send `u` to its literal set `(1)`, retaining `m=q/gcd(u,q)` |
| preserved | owner identity, gcd stratum, strict endpoint, and physical zero-cochain realization |
| destroyed | nonzero centres, affine tooth choices, and every nonzero cochain |
| required sidecar beyond `q=28` | the full short symmetric residue block, not merely unit/nonunit type |
| cheapest hostiles | the strict `q=14` endpoint and the `q=29` multiplication-by-two cycle |

The closest positive control is the canonical `q=15` rank-six edge from the
independent triphase probe.  The corrected near miss is the temptation to
continue `(3)` past 28: `q=29` identifies the first failed implication and its
replacement graph.

## 3. Exact owner classification

Write `g,m` as in `(5)` and put

```text
a=u/g mod m.
```

Then `a` is a unit modulo `m`, and direct cancellation gives

```text
||u ell/q||=||a ell/m||.                              (6)
```

Equivalently, if

```text
R_m={r in {0,...,m-1}:14 min(r,m-r)<m}                (6a)
```

and `pi_m:Z/qZ->Z/mZ` is reduction, then the general centered quotient-block
model is

```text
D_(q,u)(0)=pi_m^(-1)(a^(-1) R_m).                    (6b)
```

The classification below is the exact simplification of `(6b)` in the
stated range, not a replacement definition.

### Unit owners

If `g=1`, then `m=q`.  For `15<=q<=28`,

```text
1<q/14<=2.                                            (7)
```

Thus the only integral residues of least absolute value strictly smaller
than `q/14` are `0,+1,-1`; when `q=28`, the residues `+2,-2` are equality
points and remain excluded.  Hence

```text
D_(q,u)(0)={0,+u^(-1),-u^(-1)}                       (8)
```

for every unit speed `u`.  Each unit owner therefore covers zero and exactly
one sign-pair of unit sheets.

### Nonunit owners

If `g>1`, then `m` is a proper divisor of `q`, so throughout the stated range

```text
2<=m<=q/2<=14.                                        (9)
```

Every nonzero residue modulo `m` has circular distance at least `1/m`, hence
at least `1/14`.  Strictness in `(1)` leaves only residue zero in `(6)`, and
therefore

```text
D_(q,u)(0)={ell:m|ell}=m Z/qZ.                       (10)
```

This set contains no unit sheet.  Conversely, `(8)` contains no nonunit
sheet other than the shared zero.  The unit and nonunit cover obligations are
therefore disjoint.

## 4. Lower bound

The `phi(q)` unit sheets form `phi(q)/2` disjoint sign-pairs.  By `(8)`, each
unit owner covers exactly one of them, so every cover needs at least
`phi(q)/2` unit owners.

Now let `p<q` be a prime divisor of `q` and inspect the nonzero sheet
`ell=p`.  It cannot be covered by a unit owner.  If a nonunit owner covers it,
then `(10)` gives `m|p`.  Transversality gives `m>1`, so primality forces

```text
m=p.                                                  (11)
```

Distinct prime divisors force distinct quotient orders and hence distinct
owners.  This contributes at least one further owner for every proper prime
divisor of `q`, proving the lower bound in `(3)`.  Composite-kernel owners do
not evade it: a kernel `mZ/qZ` can hit the prime sheet `p` only when `m=p`.

## 5. Sharp construction

For every unit sign-pair `{+ell,-ell}`, choose the positive representative of
a speed satisfying

```text
u == ell^(-1) mod q.                                  (12)
```

Equation `(8)` covers that pair, and the chosen pairs partition all units.
For every proper prime divisor `p|q`, add the speed

```text
u_p=q/p.                                              (13)
```

Here `gcd(u_p,q)=q/p`, so its quotient order is `m=p` and `(10)` gives the
whole prime kernel `pZ/qZ`.  Every nonunit sheet has a common prime divisor
with `q`, hence belongs to at least one of these prime kernels.  The unit
owners and prime-kernel owners therefore cover all sheets using exactly the
number in `(3)`.

This also identifies the equality anatomy.  Every chosen unit owner has a
private unit sign-pair, while the sheet `p` is private to the quotient-order
`p` task among all nonunit kernels.  Removing any owner from the displayed
construction destroys the cover.

## 6. The `q=15` trimode edge

At `q=15`, the unit sign quotient has four classes, represented by

```text
1,2,4,7.
```

The four unit speeds `1,2,4,7` cover their four unit sign-pairs.  The speeds

```text
5=15/3,                 3=15/5                       (14)
```

cover respectively the kernels `3Z/15Z` and `5Z/15Z`, whose union is the set
of zero divisors.  Thus

```text
{1,2,3,4,5,7}                                         (15)
```

is a rank-six cover, and the lower bound proves it minimal.  At time zero its
six centres coincide, so every THM-3398 cochain gap is zero.  The four unit
owners fire the three-phase blocks `{0,+u^(-1),-u^(-1)}`, explaining
conceptually the canonical rank-six zero-cochain triphase edge and the
private-sheet obstruction found by the earlier finite probe.

## 7. Sharp boundaries

### The lower boundary `q=14`

For a unit speed, the first nonzero residue has distance exactly `1/14`, so
strictness makes

```text
D_(14,u)(0)={0}.                                      (16)
```

Every nonunit owner still has the kernel form `(10)` and therefore covers
only nonunit sheets.  Even the union of all transverse residue-speed types
misses the six units

```text
1,3,5,9,11,13.                                       (17)
```

Consequently no centered transverse cover exists at `q=14`.  The lower
endpoint of the theorem is genuine.

### The upper boundary `q=28`

For every unit speed, the two candidate sheets with phase residues `+2,-2`
satisfy

```text
||+2/28||=||-2/28||=1/14.                            (18)
```

They are excluded, so `(8)` and the rank formula remain valid at 28.

### The first failure `q=29`

At `q=29`, the allowed phase residues are now

```text
0,+1,-1,+2,-2.                                       (19)
```

On the fourteen unit sign-pairs

```text
G=(Z/29Z)^*/{+1,-1},                                 (20)
```

an owner with `x=[u^(-1)]` covers the two vertices `x` and `2x`.  Thus owner
types are the edges of the multiplication-by-two graph.  Exact modular
arithmetic gives

```text
2^14=-1 mod 29,
2^j is not +1 or -1 mod 29 for 1<=j<14.              (21)
```

Hence multiplication by two is one cycle `C_14` on `(20)`.  Seven edges are
necessary because one edge covers at most two vertices, and alternating
edges cover the entire even cycle.  Therefore

```text
r_0(29)=7,                                            (22)
```

whereas the obsolete one-sign-pair continuation of `(3)` gives fourteen.
The first failed implication is precisely `(8)`: a unit owner has acquired a
second sign-pair.  The cycle is the faithful quotient; orienting or completing
it to a tournament would discard its edge-cover predicate.

## 8. Exact verification and scope

The standard-library companion computes `(1)` directly from

```text
14 min(u ell mod q,-u ell mod q)<q,                  (23)
```

At centre zero the danger set depends only on `u mod q`, and transversality
excludes residue zero, so `u=1,...,q-1` is the complete owner-type universe.
The search uses neither fractions nor the classification formula.
It first checks `(6b)` on all `328` transverse residue owners for
`q=14,...,29`.  It then checks all `287` owners for `q=15,...,28` against
`(8)` and `(10)`, verifies the fourteen canonical constructions and every
single-owner deletion, and independently solves all fourteen set-cover
problems by exhaustive combinations of the distinct literal coverage masks.
It tests `9,828` subsets before the fourteen first covers; every lower rank is
therefore exhausted.  Discarding duplicate masks is lossless because a second
copy cannot enlarge a union.

The same literal engine proves the `q=14` noncover, checks all `24` signed
`q=28` endpoint incidences, and independently obtains both the literal
`q=29` cover rank and the `C_14` minimum edge cover.  There are no floating
literals or optimization-dependent assertions.  Reproduce with

```text
python 04-computation/lrc_centered_transverse_sheet_cover_rank_thm3401.py
python -O 04-computation/lrc_centered_transverse_sheet_cover_rank_thm3401.py
```

Both runs LF-normalized-byte-match the stored output.

## 9. Open continuation

The natural all-`q` question exposed by the two boundaries is whether, for
every `q>=15`, a centered transverse cover of rank at most six exists **iff**

```text
15|q or 16|q or 18|q or 20|q or 24|q.                (24)
```

The implication from divisibility in `(24)` to rank at most six follows from
the five base covers above and THM-3398's exact dilation/pullback.  The reverse
implication is **OPEN**.  A concurrent unpromoted scratch census reported no
reverse-implication exception through `q=200`; this theorem neither audits nor
depends on that census.  Proving or refuting the reverse implication requires
the full short-block quotient model `(6b)`, because the unit/nonunit collapse
used in `(8)`--`(10)` ends at 28.

This theorem concerns only common source centre zero, a fixed-gauge
sub-slice of THM-3398's zero-cochain locus.  It does not classify covers at a
mobile nonzero common centre, the full finite-mode clutter, or any
nonzero-cochain cover; it does not decrement the refined ledger or prove
LRC(14).
