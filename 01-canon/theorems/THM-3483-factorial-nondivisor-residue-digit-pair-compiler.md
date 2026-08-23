---
id: THM-3483
title: "Nondivisor residue digit polygons and pair-only factorial closure"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED.  For every
  prime p not dividing d, one finite residue rho_p(n,j,d) decides whether the
  raw Kummer--Legendre height of [v^j]A_n^d is exact.  If rho is nonzero at
  every extreme raw-hull vertex, the entire actual Newton polygon equals the
  coefficient-free raw digit hull.  Applying this pairwise at an adaptive
  nondivisor prime for each packet eliminates all five THM-3475 survivors,
  upgrading its 33/38 application to a pair-only 38/38 closure through d=2600.  Two
  exact implementations then apply the divisor and residue compilers to every
  d in 2606..4000: 934 rows take an inherited exit, divisor places close 420
  of the 461 residuals, and admissible nondivisor residues close the remaining
  41.  Hence the exact-support quadratic boundary is r<=3998.  This is a
  necessary-degree compiler and finite exact application, not a factor
  construction or arbitrary-support FC/SFC theorem.  Two further dual exact
  blocks close every row through d=10000: across 2606<=d<=10000, 4730/7395
  rows take inherited exits, divisor places close 2384/2665 residuals, and
  admissible nondivisor residues close the remaining 281.  A separately
  audited 51-row window through d=10030 then closes 9 rows by inherited exits,
  36 by divisor barcodes and 6 by admissible nondivisor barcodes.  Thus the
  current exact-support boundary is r<=10028; d=10031 is untested.
audit: >
  The residue congruence, extreme-vertex hull criterion, p=2 boundary, zero
  coefficients, slope-zero typing, and five large pair ledgers were
  independently proof- and hostile-audited.  A pinned exact companion checks
  23,395 coefficients in 746 profiles; all 665 rho-admissible profiles have
  the predicted actual hull, with no mismatch.  It freezes the first sharp
  rho=0 hostile and exact empty intersections for all five former packets.
  Separate primary and self-contained block companions agree on the complete
  1395-row record through d=4000, including 20 inadmissible profiles that are
  skipped rather than used.  A further independent window audit uses two
  no-import formula routes and 3,700 direct coefficient controls.  It repairs
  the scope of the d=9996 hostile: the THM-3475 all-divisor barcode
  intersection is incomplete there, while p=19 is only the first admissible
  killer in the declared bank.  Normal and optimized replays are byte-identical.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
  - THM-3478-factorial-dual-seven-exit-block-closure-through-r2603
related:
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
script: 04-computation/factorial_nondivisor_residue_digit_pair_compiler_thm3483.py
output: 05-knowledge/results/factorial_nondivisor_residue_digit_pair_compiler_thm3483.out
script_sha256: 9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723
output_sha256: 5bfc1b06cda024080d1c4c977511bb6f0f30ecdbf7ad6e8072a1052c072bebc9
semantic_sha256: f80c046942d62a8a6b6f3802d224cc47944568a7e9f3ef245d43baf91a4031c4
block_scripts:
  - 04-computation/factorial_adaptive_rho_block_4000_thm3483.py
  - 04-computation/factorial_adaptive_rho_block_4000_independent_audit_thm3483.py
block_outputs:
  - 05-knowledge/results/factorial_adaptive_rho_block_4000_thm3483.out
  - 05-knowledge/results/factorial_adaptive_rho_block_4000_independent_audit_thm3483.out
block_script_sha256:
  - b58fd73c96929fa287b09addfe86463d0ca9998a31e37dbd964d1456828396e4
  - 0b858f7b1154a3ee2dec43bf5238f7e6f24b524e9dc2b9f70f5b497fcef58934
block_output_sha256:
  - 48733decf5874c197b989d7731f0864082a65fe7b9056af67e0149d1e8a94896
  - 5c02ee06a8b14909b1b6677e1d853057eeb285ab5699cd354e5e30304533b329
block_semantic_sha256: 95d1c233d59d00c38ce456fa7c5f5e248414e01b5ba9dc2ae9f61725d6c19dbd
extension_scripts:
  - 04-computation/factorial_adaptive_rho_block_6000.py
  - 04-computation/factorial_adaptive_rho_block_6000_independent_audit.py
  - 04-computation/factorial_adaptive_rho_block_10000.py
  - 04-computation/factorial_adaptive_rho_boundary_10000_gate.py
extension_outputs:
  - 05-knowledge/results/factorial_adaptive_rho_block_6000.out
  - 05-knowledge/results/factorial_adaptive_rho_block_6000_independent_audit.out
  - 05-knowledge/results/factorial_adaptive_rho_block_10000.out
  - 05-knowledge/results/factorial_adaptive_rho_boundary_10000_gate.out
extension_script_sha256:
  - b65edcf2870714ca57456b8297afdd05284a09b302ec4b84d2e57829520c94d1
  - d416cb2955fd745394cf1043ac8c2eba28a6a97beb264dd9cbe9919ed8c96724
  - 105e62698d0c3cf0066a100e9d205a5c1f1c31e64cfdba0dc2fe23decd8f0eba
  - ee57eb8df0c25aa6bc0243f66f7f3b21affefa75f3b347b8ec19e10f88a2d5b3
extension_output_sha256:
  - ab629edc04e31d1889741688897bfe60f5249df60b64116391217393962b1ddf
  - 8d6adbcaa14c85d022f726db97a80a7bafc1288505ccde720b3f5fbf6ee2a922
  - 18b131aed2f380b1c1bace8beeb8488ced0e24599f4b7484d66a14e5869c0d22
  - 625e0e401e914a5bfa30d139e6309b375c4461d83c4dfe76755179f0d495fea5
extension_block_semantic_sha256:
  - 7f8ab74ae9fae027f32fd7eabaf0338c217319e274594bd603859a1bbcca28bd
  - d90179fdebd48dd82cd368b957c9602fbd287774287de0fecb73b4a84dca69f3
extension_combined_semantic_sha256: 4a364bbafdfef0dc6d905063c9570b987e00eca2506506b8160ce24e453878bf
window_scripts:
  - 04-computation/factorial_d10001_local_rule_window_20260823.py
  - 04-computation/factorial_d10001_local_rule_window_independent_audit_20260823.py
window_outputs:
  - 05-knowledge/results/factorial_d10001_local_rule_window_20260823.out
  - 05-knowledge/results/factorial_d10001_local_rule_window_independent_audit_20260823.out
window_script_sha256:
  - 03af13ffa39a785c044c24ee589c4e156375ea21d87178e72c69651de4241ffe
  - 6e467e8456d12d47d849d24bf7771ad330b1ec5a8e059b466c366650da443fe9
window_output_sha256:
  - ea7d26fea7a646e7c960938e0fc737d9c82a734df85edd125a2a8c9b61785263
  - dc958efde34b5d4c5ea0270ab36336310dcc6a90b85c2501d3e86fdeeea8526d
window_semantic_sha256:
  - f2d02ba3f51d6dd5b47a0976c775b307a3b73b64e01affaecf198f349a635213
  - a26de54d7bf92d5c7e02820f9ca46d5a50ab24729e180a388781ddafe9dd4e2f
hash_basis: raw bytes
---

# THM-3483 -- nondivisor residues lift raw factorial polygons

**PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED.**

## 1. Residue-lift theorem

Let `L(x^r)=r!` and define

```text
A_n^(d)(v)=L((d-x+v x^2)^n)=sum_(j=0)^n c_(n,j)v^j.       (1)
```

Fix a prime `p` with `p` not dividing `d`.  For `0<=j<=n`, put

```text
m=n-j,        m_0=m mod p,        j_0=j mod p,
d_0=d mod p in F_p^*.                                      (2)
```

Define the finite residue

```text
rho_p(n,j,d)
 =sum_(ell=0)^(m_0)
   C(m_0,ell)(-d_0^(-1))^ell (2j_0+1)_ell       in F_p.   (3)
```

Here `(y)_ell=y(y+1)...(y+ell-1)` and `(y)_0=1`.  Define the
raw coefficient height

```text
w_p(n,j)=nu_p C(n,j)+nu_p((2j)!).                          (4)
```

Then

```text
nu_p(c_(n,j)) >= w_p(n,j),                                (5)
```

with infinity assigned to a zero coefficient, and

```text
nu_p(c_(n,j))=w_p(n,j)  iff  rho_p(n,j,d)!=0.             (6)
```

Let `RawNP_p(A_n^d)` be the lower convex hull of the finite points
`(j,w_p(n,j))`.  If (3) is nonzero at every **extreme vertex** of this raw
hull, then

```text
NP_p(A_n^d)=RawNP_p(A_n^d).                               (7)
```

Conversely, if (3) vanishes at an extreme raw vertex, the height at that
index strictly rises or becomes infinite, so the raw polygon is not the
actual polygon.  The residue test does not by itself compute the amount of
the rise.

Thus a base-`p` digit hull plus a finite residue table replaces every large
integer coefficient whenever its extreme vertices are rho-admissible.

## 2. Proof of the residue formula

Choose the `j` factors contributing `v x^2` in (1), expand the remaining
`m=n-j` factors, and apply `L`.  Directly,

```text
c_(n,j)=C(n,j)(2j)! Z_(n,j),                              (8)

Z_(n,j)=sum_(ell=0)^m
 C(m,ell)d^(m-ell)(-1)^ell(2j+1)_ell.                    (9)
```

The integer `Z_(n,j)` proves (5).  Divide (9) by the unit `d^m` modulo `p`.
Every term with `ell>=p` vanishes because any `p` consecutive integers in
the rising product contain a multiple of `p`.  For `0<=ell<p`, Lucas gives

```text
C(m,ell)=C(m_0,ell) (mod p),                              (10)
```

which is zero when `ell>m_0`; the remaining rising product depends only on
`j_0`.  Therefore

```text
d^(-m) Z_(n,j)=rho_p(n,j,d) (mod p).                     (11)
```

Since `d` is a unit modulo `p`, (6) follows.

For (7), every actual point lies on or above its raw point by (5), while
each extreme raw vertex is an exact actual point by (6).  Those exact
endpoints pin every raw edge from above and below, so the actual lower
envelope is the same polygon.  If an extreme raw vertex has rho zero, its
unique point rises; extremality prevents the same raw polygon from surviving.

## 3. Pair-degree compiler

Apply (3)--(7) independently to

```text
F=A_(d-2)^(d),                 G=A_(d-1)^(d).             (12)
```

When both are rho-admissible, their exact Newton polygons require no
coefficient construction.  THM-3152 then converts every common finite slope
`a/b` in lowest terms into degree contributions in `b Z_(>=0)`, bounded by
the minimum horizontal capacity of the two slope faces.  The Minkowski sum
over all common slopes, together with the separately typed coordinate-root
capacity, is the exact local necessary degree barcode.

A common slope `0` is a unit-root block with denominator one.  It must be
retained in this sum and must not be confused with coordinate-root capacity.

THM-3161 and THM-3475 are rho-one specializations: their exact digit anchors
have either `m_0=0` or a rising-factor residue that forces (3) to be one.  The
new content is the arbitrary-`n,d`, nondivisor-prime residue sidecar and its
exact failure flag.

## 4. Five coefficient-free residual closures

THM-3475's divisor-only pair compiler closes `33` of the `38` seven-exit
residuals in `2501<=d<=2600` and leaves five exact necessary-degree packets.
At the adaptive nondivisor primes below, every extreme raw vertex of both
rows (12) is rho-admissible.  The displayed exact common blocks have format
`(slope,capacity,reduced denominator)`:

```text
d=2516, p=7:
 blocks=(0,1,1),(2/7,14,7),(16/49,98,49),(800/2401,2401,2401)
 old packet={503,1006,1509,2012};                         (13)

d=2564, p=13:
 blocks=(0,1,1),(2/13,26,13),(28/169,338,169),(366/2197,2197,2197)
 old packet={466,699,1165,1631,1864,2097,2330};          (14)

d=2571, p=7:
 blocks=(2/7,21,7),(16/49,147,49),(800/2401,2401,2401)
 old packet={2056};                                       (15)

d=2576, p=13:
 blocks=(2/13,39,13),(28/169,338,169),(366/2197,2197,2197)
 old packet={103,206,...,2472} (24 multiples of 103);     (16)

d=2586, p=7:
 blocks=(0,1,1),(2/7,21,7),(16/49,147,49),
        (800/2401,2401,2401),(3/7,14,7)
 old packet={47,141,188,235,282,329,
             2209,2256,2303,2350,2397,2444,2491,2538}.   (17)
```

In each row, the exact local degree set has empty intersection with the old
packet.  Including degree zero, the five local-set sizes are respectively

```text
36,36,32,24,96,                                           (18)
```

and both constant coefficients are rho-units, so the coordinate-root
capacity is zero in all five cases.  This proves the five pair-only closures,
without `R_1` and without constructing an integer coefficient array of degree
about `2500`.

Consequently THM-3475's finite application strengthens from `33/38` to

```text
38/38 pair-only closures through d=2600,                  (19)
```

or `r<=2598` after the inherited exits.  THM-3478 separately extends the
contiguous exact-support boundary through `r=2603`; this five-row application
improves the mechanism for its first difficult rows, while Section 5 below
uses the same residue theorem to advance the numerical endpoint.

## 5. Adaptive exact block through `d=4000`

THM-3478 closes the contiguous exact-support quadratic range through
`d=2605`, or `r=d-2<=2603`, and leaves at `d=2606` the divisor-only packet

```text
{521,1042,1563,2084}.                                   (20)
```

Two independent exact implementations now apply the same typed policy to
every row `2606<=d<=4000`: first inherit the seven elementary exits, then
intersect the complete THM-3475 pair barcodes at every divisor prime of
`d-1`, and finally try only rho-admissible THM-3483 pair polygons at ordered
small nondivisor primes.  Their common census is

```text
total rows                         1395
inherited seven-exit rows           934
post-exit residuals                  461
closed by divisor pair ledgers       420
closed by a rho-admissible prime       41
survivors                               0.             (21)
```

The rho-killer histogram is

```text
p=3:6, p=5:1, p=7:17, p=11:9, p=13:6, p=17:2.          (22)
```

No rho prime above `47` is invoked, and every actual rho killer is at most
`17` (divisor places may of course be larger).  The programs record `20`
inadmissible raw profiles and skip every one; a raw barcode is never used when
an extreme vertex has rho zero.

The former first residual (20) illustrates the mechanism without large
coefficients.  At `d=2606`, `p=3`, the exact raw/actual hulls are

```text
F=A_2604^2606:
 ((0,0),(3,2),(12,10),(93,90),(336,332),
  (2523,2518),(2550,2545),(2604,2600));

G=A_2605^2606:
 ((0,0),(1,0),(4,2),(13,10),(94,90),(337,332),
  (2524,2518),(2551,2545),(2605,2600)).                 (23)
```

Every `F` vertex has rho one; the constant `G` vertex has rho two and every
other `G` vertex has rho one.  The common blocks are

```text
(2/3,3,3), (8/9,9,9), (80/81,81,81),
(242/243,243,243), (2186/2187,2187,2187),
(1,27,1), (55/54,54,54).                               (24)
```

The component omitting the `2187` block has maximum degree `417`, while a
component using it has minimum degree `2187`; hence (20) has empty
intersection with the local barcode.

The first two block implementations serialize the same full semantic record with
SHA-256

```text
95d1c233d59d00c38ce456fa7c5f5e248414e01b5ba9dc2ae9f61725d6c19dbd. (25)
```

Two additional dual exact audits continue without a gap.  On
`4001<=d<=6000`, their common census is

```text
total rows                         2000
inherited seven-exit rows          1272
post-exit residuals                 728
closed by divisor pair ledgers      600
closed by a rho-admissible prime    128
survivors                             0.                    (26)
```

The block semantic SHA-256 is
`7f8ab74ae9fae027f32fd7eabaf0338c217319e274594bd603859a1bbcca28bd`.
On `6001<=d<=10000`, the exact census is

```text
total rows                         4000
inherited seven-exit rows          2524
post-exit residuals                1476
closed by divisor pair ledgers     1364
closed by a rho-admissible prime    112
survivors                             0.                   (26a)
```

Its semantic SHA-256 is
`d90179fdebd48dd82cd368b957c9602fbd287774287de0fecb73b4a84dca69f3`;
the lightweight cross-block gate pins both records with combined SHA-256
`4a364bbafdfef0dc6d905063c9570b987e00eca2506506b8160ce24e453878bf`.
Across the full audited range `2606<=d<=10000`, the totals are therefore

```text
7395 rows; 4730 inherited exits; 2665 residuals;
2384 divisor closures; 281 rho closures; 0 survivors.       (26b)
```

The actual rho-killer histogram in that full range is

```text
p=3:31, p=5:23, p=7:81, p=11:83, p=13:31,
p=17:18, p=19:9, p=23:4, p=29:1.                          (26c)
```

The unique `p=29` row is a useful boundary hostile.  At `d=6518`, divisor
places leave

```text
{3087,3430,4802,5145,5488,5831}.
```

Primes `13,17,19` retain the last candidate `5831`, and `p=23` is
rho-inadmissible; the admissible `p=29` barcode has a gap from `5669` to
`5887` and kills it.  Thus the earlier observed killer bound `p<=19` was only
finite-range evidence, and even `p<=23` is false on the present range.

A separate exact window covers `9980<=d<=10030`.  Its 51 rows split into
9 inherited exits and 42 residuals; THM-3475 divisor barcodes close 36 and
THM-3483 admissible nondivisor barcodes close the remaining 6.  No row
survives.  The anchor `d=10001` closes by the intersection of its `p=2` and
`p=5` divisor packets.  Five steps away, `d=9996` leaves degree `3998` after
all THM-3475 divisor barcodes, and the first admissible killer in the ordered
bank is `p=19`, by the exact gap

```text
3059 < 3998 < 6859.                                        (26d)
```

This refutes completeness of that specific all-divisor barcode criterion,
not every possible divisor-place argument.  Later admissible primes also kill
`3998`; no uniqueness claim for `p=19` is made.

Consequently every exact-support quadratic window with

```text
1<=r<=10028                                                  (26e)
```

has a nonzero member among its factorial moments at exponents
`r,r+1,r+2`.  This is a FINITE-EXACT application of the proved polygon
compilers.  The next untested row is `d=10031`, or `r=10029`; it is not a known
survivor.

Equivalently, if

```text
L(t^m)=m!,                   q(t)=a+bt+ct^2,             (26f)
```

with `abc!=0`, then for every `1<=r<=10028` the three values
`L(q^r),L(q^(r+1)),L(q^(r+2))` cannot all vanish.  This does not cover a
missing coefficient, translated or arbitrary support, all of `SFC(1)`,
`SFC(3)`, or `FC(3)`.

## 6. Exact audit and the sharp hostile

The companion's coefficient-control universe is

```text
4<=d<=60;
p in {2,3,5,7,11,13,17,19},       p not dividing d;
n in {d-2,d-1}.                                          (27)
```

It checks `23,395` exact coefficients in `746` polygon profiles.  Exactly
`665` profiles are rho-admissible at every extreme raw vertex; all `665`
actual polygons equal their raw polygons, with zero mismatch.  An independent
direct factorial expansion reproduces the congruence and every admissible
profile on the same universe.  The companion also freezes every raw hull,
vertex residue, common block, old packet, and empty intersection in
(13)--(17).

The first sharp hostile is

```text
d=4,       p=5,       F=A_2^(4),
coefficients=(10,4,24),
raw hull=((0,0),(2,0)),
vertex rhos=((0,0),(2,1)),
actual hull=((0,1),(1,0),(2,0)).                         (28)
```

The vanished residue at the raw constant vertex raises it and creates a new
negative edge.  Raw digit weights without rho therefore give a false polygon.

At `p=2`, necessarily `d_0=1`, and (3) reduces to

```text
rho_2(n,j,d)=1  iff  m=n-j is even.                       (29)
```

The companion checks this for odd `3<=d<=19` and every `0<=j<=n<40`.
This is a boundary rule, not an unqualified binary exactness theorem.

Reproduce the byte-identical normal and optimized transcripts with

```bash
python3 04-computation/factorial_nondivisor_residue_digit_pair_compiler_thm3483.py
python3 -O 04-computation/factorial_nondivisor_residue_digit_pair_compiler_thm3483.py
```

The block companions are reproduced by

```bash
python3 04-computation/factorial_adaptive_rho_block_4000_thm3483.py
python3 -O 04-computation/factorial_adaptive_rho_block_4000_thm3483.py
python3 04-computation/factorial_adaptive_rho_block_4000_independent_audit_thm3483.py
python3 -O 04-computation/factorial_adaptive_rho_block_4000_independent_audit_thm3483.py
```

The extension blocks and lightweight combined gate are reproduced by

```bash
python3 04-computation/factorial_adaptive_rho_block_6000.py
python3 -O 04-computation/factorial_adaptive_rho_block_6000.py
python3 04-computation/factorial_adaptive_rho_block_6000_independent_audit.py
python3 -O 04-computation/factorial_adaptive_rho_block_6000_independent_audit.py
python3 04-computation/factorial_adaptive_rho_block_10000.py
python3 -O 04-computation/factorial_adaptive_rho_block_10000.py
python3 04-computation/factorial_adaptive_rho_boundary_10000_gate.py
python3 -O 04-computation/factorial_adaptive_rho_boundary_10000_gate.py
```

The audited `d=10001` neighbourhood is reproduced by

```bash
python3 -B 04-computation/factorial_d10001_local_rule_window_20260823.py
python3 -B -O 04-computation/factorial_d10001_local_rule_window_20260823.py
python3 -B 04-computation/factorial_d10001_local_rule_window_independent_audit_20260823.py
python3 -B -O 04-computation/factorial_d10001_local_rule_window_independent_audit_20260823.py
```

## 7. Failure boundaries and information contract

- The condition `p` not dividing `d` is load-bearing: (11) divides by `d^m`.
- Rho nonvanishing is required only at extreme raw-hull vertices.  Vanishing
  at a nonextreme contact need not change the polygon.
- Rho zero proves a strict lift or zero coefficient, but gives neither the
  extra valuation nor the repaired hull.
- A local degree address is necessary, not a factor.  Only an empty
  intersection proves coprimality.
- Nonemptiness of every THM-3475 divisor barcode does not prove that all
  conceivable divisor-place refinements fail; the `d=9996` hostile is scoped
  to the displayed necessary-degree criterion.
- Zero slope and coordinate-root capacity are different data.  The former is
  retained in (13), (14), and (17); the latter is zero in all five rows.
- The result concerns the exact quadratic resonance pair.  It proves no
  arbitrary-support `SFC(1)`, `SFC(3)`, or `FC(3)` case.

The information contract is

```text
source:    raw Kummer--Legendre heights of A_n^d at a prime p not dividing d
target:    the actual Newton polygon, then a necessary pair-degree barcode
map:       extreme vertex -> finite residue rho -> unit lift -> polygon
preserved: vertex index, prime, slope, capacity, denominator, root-block type
lost:      higher valuation when rho=0, factor existence, root residues
sidecar:   (m mod p,j mod p,d mod p) and the complete extreme-vertex list
hostile:   A_2^4 at p=5, where rho=0 creates a new negative edge.             (30)
```

**QED.**
