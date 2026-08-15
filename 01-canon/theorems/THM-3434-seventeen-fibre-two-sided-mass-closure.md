---
id: THM-3434
title: "Seventeen-fibre two-sided mass closure and odd half-rank-seven support"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The residual 17^a and
  5*17^a half-twist towers are excluded by an order-sensitive two-sided fibre
  invoice; consequently the ordinary transverse literal half-twist rank is at
  most seven exactly on odd multiples of 9,11,13,15,23,25,29,or51.  This is
  not an arbitrary-time or physical LRC(14) result.
source: root seventeen-adic fibre session, 2026-08-15
audit: independent period descent, dependency and live-main THM-3429 repin delta, transverse universe, strict endpoints, quotient-density maxima, 2/3 fibre law, owner count, two-sided invoice, tower arithmetic, Q85 p17/p5 fixed-fibre budget, divisor support, conditional fixed-chart corollary, positive Q51 hostile, dual clean-room exact companions, normal/-O/stored replay, hashes, AST/security, documentation, and scope CLEAN; corrected an uncommitted sidecar that had failed to charge even full-order owners on the covered Q85 fibre
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3421-prime-half-twist-rank-seven-classification
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
  - THM-3432-order-two-fixed-half-parity-transplant
related:
  - THM-3425-half-twist-rank-six-primitive-breaker-profile-closure
  - THM-3428-rough-maximal-order-half-twist-rank-seven-exclusion
script: 04-computation/lrc_seventeen_fibre_two_sided_mass_closure_thm3434.py
output: 05-knowledge/results/lrc_seventeen_fibre_two_sided_mass_closure_thm3434.out
script_sha256: 910349a5c5c4178b068cfe5427d01d7b7efadc54875b01bac87b162d0eccf996
output_sha256: 3adbdc6c6937b13a8c6e7dfc1ab8b27ba001af6860341fb082e2a01a65dcf0dc
semantic_sha256: 13984be8cee76a72b5d1ff98011538e56eb9c5f86ff3db41ed80f7c536189f2d
independent_script: 04-computation/lrc_rank7_residual_17adic_tower_probe_20260815.py
independent_output: 05-knowledge/results/lrc_rank7_residual_17adic_tower_probe_20260815.out
independent_script_sha256: b70570b1baecb719399e7819fce772d2aad7d145e8e44811164046d623cff1fa
independent_output_sha256: 6b10dec50965d932b73d4c6f16d79c1aa00a4b733b425dd33b1a64595a5ce8f1
independent_semantic_sha256: 2939e6278cc64659599ec030d620cd52b12574cbfb2f531a3e34d8ebf8df1fcd
hash_basis: LF-normalized bytes
---

# THM-3434 -- seventeen-fibre two-sided mass closure and odd half-rank-seven support

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

For odd `Q>=3` and a residue `r` modulo `2Q`, put

```text
B_(Q,r)={ell in Z/QZ: ||r(2ell+1)/(2Q)||<1/14}.          (1)
```

Only **transverse** residues are admitted:

```text
Q does not divide r.                                    (2)
```

Thus the zero residue, whose mask would be all sheets, and the order-two
residue `Q`, whose mask is empty, are outside the owner universe.  Let
`rho_H(Q)` be the minimum cardinality of a transverse family whose masks in
`(1)` cover all sheets.  Then

```text
rho_H(Q)<=7
iff some d in {9,11,13,15,23,25,29,51} divides Q.       (3)
```

In the lower-base-free sector

```text
9,11,15,23,25 do not divide Q,                          (4)
```

the three positive divisor atoms `13,29,51` have exact rank seven.  In
particular, the two residual lanes left by THM-3429,

```text
Q=17^a       (a>=2),
Q=5*17^a     (a>=1),                                   (5)
```

admit no transverse seven-cover.

Let `rho_(Z,2)(2Q)` denote the minimum fixed-zero rank among families on
`2Q` which are required to contain the unique order-two residue `Q`.
THM-3432 now gives the exact conditional fixed-chart corollary

```text
rho_(Z,2)(2Q)<=8
iff some d in {9,11,13,15,23,25,29,51} divides Q.       (3a)
```

Formula `(3)` is an exact Boolean half-twist classification.  It supplies no
common physical time, runner-current certificate, arbitrary mode centre, or
LRC(14) decrement.

## 2. Inheritance and connection contract

THM-3416 classifies ordinary half-twist cap six.  THM-3421 closes the prime
cap-seven layer.  THM-3429's audited prime-fibre descent reduces every
remaining lower-base-free joint-period candidate, away from divisors
`13,29,51`, to `(5)`.  The present proof closes exactly those two towers.

| field | exact connection |
|---|---|
| source | a mixed-order transverse half-twist cover on `Q=17N` |
| target | one descended inactive mask on `N` plus its missed `17`-fibres; then the THM-3432 order-two fixed chart |
| map | prime-fibre projection `Z/QZ -> Z/NZ` and incidence summation over covered/missed fibres |
| preserved | strict endpoints, owner count, quotient order, block size, and every fibre incidence |
| destroyed | the affine lift character locating each active three-set inside a fibre |
| needed sidecar | the inactive quotient order, which controls its global density |
| corrected near miss | the false uniform estimate `|B|<=ceil(Q/7)` for lower-order pullbacks |
| positive boundary | `Q=51`, where coindex three makes the six active blocks large enough |
| cheapest hostile | `Q=85`, where the first mass bound is inconclusive but the five-fibre invoice closes |

The least-used coordinate is not the active lift cocycle itself but the
**two-sided mass** it cannot change: an active block must spend at least two
points on every covered base fibre as well as enough points on every missed
one.

## 3. Period descent and the global reduction

For a selected residue define its quotient order

```text
m_Q(r)=Q/gcd(Q,r).                                      (6)
```

Given any cover `R`, put

```text
L=lcm_(r in R) m_Q(r),       h=Q/L.                    (7)
```

Every `m_Q(r)` divides `L`, hence `h|r`; write `r=hs`.  Direct cancellation in
the phase word gives

```text
B_(Q,hs)=pi^(-1) B_(L,s),    pi:Z/QZ -> Z/LZ.          (8)
```

Moreover

```text
L/gcd(L,s)=m_Q(r),                                     (9)
```

so the descended family remains transverse, covers literally, and has joint
period `L`.  Thus it is enough to classify joint-period covers.

If one of the five integers in `(4)` divides `L`, it divides `Q`.  Otherwise
THM-3421 handles prime `L`, while THM-3429 reduces a composite `L` with no
divisor `13,29,51` to `(5)`.  Sections 5--7 exclude both lanes.  Conversely,
the rank-seven atoms

```text
13: (1,2,3,5,7,9,11),
29: (1,5,7,8,12,13,22),
51: (1,11,12,18,23,34,35)                             (10)
```

pull back along every multiple by scaling all residues.  THM-3416 supplies
the five lower divisor atoms in `(4)` and proves that the rows in `(10)` are
not rank six in the lower-base-free sector.  This proves `(3)` once `(5)` is
closed.

## 4. Three exact size lemmas

### 4.1 Quotient pullback

If a block on `Q` has quotient order `m`, its mask is the pullback of a block
on `m`, repeated `Q/m` times.  For every odd `m`, the danger interval has
length `1/7`, hence

```text
|B_(m,s)|<=ceil(m/7).                                  (11)
```

Two order-sensitive density consequences will be used:

```text
m=17^c, c>=1:       |B_(m,s)|/m <= 3/17,               (12)
m|5*17^c, m>1:     |B_(m,s)|/m <= 1/5.                (13)
```

For `(12)`, equality is possible at `m=17`; for `m>=289`, use
`ceil(m/7)<=(m+6)/7<=3m/17`.  For `(13)`, `m=5` is equality, while every
other possible divisor has `m>=17` and
`ceil(m/7)<=(m+6)/7<=m/5`.

### 4.2 A complete active seventeen-grid

Let `Q=17N` and let `r` be `17`-active.  On every fibre of
`Z/QZ -> Z/NZ`, equation `(1)` becomes a translate of the complete
seventeen-grid inside an open arc of length `1/7`.  The first two grid points
after either endpoint lie strictly inside because

```text
2/17<1/7,
```

whereas four points would span at least `3/17>1/7`.  Therefore every active
section has exactly two or three points:

```text
2<=|B_(Q,r) intersect F_y|<=3.                         (14)
```

For a five-fibre, the corresponding upper bound is one point because
`1/5>1/7`.

### 4.3 The two-sided invoice

More generally, let `p!=7` be prime, let `Q=pN`, and let the union of all
inactive descended masks have size `h`.  A translate of the `p`-grid in the
open danger arc contains between `floor(p/7)` and `ceil(p/7)` points.  If `a`
active blocks cover every fibre missed by the inactive union, incidence
summation gives the reusable inequality

```text
p(N-h)+a floor(p/7) h <= sum_(active i) |B_i|.          (15)
```

In either lane `(5)`, THM-3429 gives a mixed `17`-coordinate with an even
inactive owner.  The inactive descended masks cannot cover the base by at
most six owners, by THM-3416.  A missed fibre needs all seventeen points, so
the upper half of `(14)` forces at least six active owners.  Hence there are
exactly

```text
six 17-active owners and one 17-inactive owner.         (16)
```

Let `h` be the size of the unique descended inactive mask.  On the `N-h`
missed fibres, the active incidence mass is at least `17(N-h)`.  On the `h`
covered fibres, the lower half of `(14)` makes those same six blocks spend at
least `12h` further incidences.  Thus every putative cover satisfies

```text
17N-5h <= sum_(six active i) |B_i|.                    (17)
```

This inequality is insensitive to the affine lift cocycle that stopped
THM-3429.  The `Q=51` positive shows why the quotient-order sidecar is still
essential: its order-seventeen active pullbacks have size nine, larger than
`ceil(51/7)=8`.

## 5. The pure `17^a` tower

Let `Q=17^a=17N`, `a>=2`.  A `17`-active residue is coprime to `Q`, so each
active block is full-order and `(11)` gives

```text
sum_i |B_i| <= 6 ceil(17N/7).                          (18)
```

The inactive block downstairs has nontrivial quotient order `17^c`; by
`(12)`,

```text
h<=3N/17.                                              (19)
```

Combining `(17)--(19)` would force

```text
274N/17 <= 6 ceil(17N/7) <= 6(17N+6)/7.               (20)
```

After multiplication by `119`, `(20)` implies

```text
184N<=612,                                             (21)
```

which is false for `N>=17`.  This closes the first lane in `(5)`.

## 6. The five-times-`17^a` tail

Let

```text
Q=5*17^a=17N,       N=5M,       M=17^(a-1).           (22)
```

For a `17`-active owner, the only possible coindices are one and five.
Consequently its block size is at most

```text
5 ceil(17M/7),                                        (23)
```

and the six active blocks have total size at most thirty times that ceiling.
The inactive base block satisfies `(13)`, hence

```text
h<=N/5=M.                                              (24)
```

Equations `(17),(23),(24)` would force

```text
80M <= 30 ceil(17M/7) <= 30(17M+6)/7,                 (25)
```

and therefore

```text
50M<=180.                                              (26)
```

This is false for `M>=17`, closing every `a>=2` case.

## 7. The exceptional boundary `Q=85`

It remains to treat `a=1` in `(22)`.  The unique even `17`-inactive block is,
up to sign, residue `34` or `68`.  It has quotient order five and descends to
the singleton reflection-fixed sheet `2` modulo five.  Thus exactly four
seventeen-fibres, containing `68` sheets, are missed.

Let `b` of the six `17`-active owners be divisible by five.  Project now to
the five-coordinate.  Those `b` owners are inactive and descend to at most
six half masks on modulus seventeen; they cannot cover that base by
THM-3416.  On a missed five-fibre each active section has at most one point,
so the remaining `7-b` owners must cover five points.  Hence

```text
b<=2.                                                   (27)
```

On the four missed seventeen-fibres, a full-order owner contributes exactly
ten incidences: at `Q=85=14*6+1` its global size is `13` or `12` according to
even or odd parity, while its section over the covered fixed base sheet has
size `3` or `2`.  A five-divisible owner contributes at most twelve.  The
total available mass is therefore

```text
10(6-b)+12b = 60+2b <=64 <68.                          (28)
```

This contradiction closes `Q=85` and the second lane in `(5)`.

## 8. Exact companion and boundaries

Run

```bash
python3 04-computation/lrc_seventeen_fibre_two_sided_mass_closure_thm3434.py
python3 -O 04-computation/lrc_seventeen_fibre_two_sided_mass_closure_thm3434.py
```

The standard-library companion checks the quotient-density formulas through
eight tower levels; every point of the relevant seventeen- and five-fibre
controls; the integer inequalities through exponent ten; all eighty
`Q=85` active residue invoices; and `176,850` period-descent cells.  An
independent complete rare-coordinate search proves that `Q=85` has no
seven-cover in `760` memoized states and `844` branches.  It also reconstructs
the three positive atoms and twelve pullbacks, and independently reproduces
the divisor support `(3)` at every odd `Q<=101`.  Normal and optimized outputs
are intended to be byte-identical.

The independently written
`lrc_rank7_residual_17adic_tower_probe_20260815.py` rederives the argument in
third-point-event coordinates, checks `Q=51,85,289` with a joint-period
solver, and exhausts `1,006,880` forced-profile combinations at `Q=85`.  Its
audit caught and repaired a provisional bookkeeping error: every even active
owner, not only an order-seventeen owner, spends a third point on the covered
reflection-fixed fibre.  After that charge no coarse `Q=85` mass profile
survives; the pair-overlap bank is retained only as a hostile sidecar.

The strict transverse condition `(2)` is load-bearing.  Admitting residue
zero makes `(1)` the full set and destroys every rank statement.  Allowing a
mobile centre, forgetting the odd half-chart, or replacing literal union by a
fractional density likewise lies outside the theorem.  No LRC(14) consequence
is claimed.

**QED.**
