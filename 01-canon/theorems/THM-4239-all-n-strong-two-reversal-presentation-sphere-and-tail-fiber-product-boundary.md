---
id: THM-4239
title: "All-n strong two-reversal presentation sphere and tail-fiber product boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every n>=3, the
  strong fixed-gauge presentations obtained by reversing exactly two arcs
  of the labelled transitive tournament form two explicit families of total
  size (n-1)^2-3. Exact all-order formulas for every type-B tail owner isolate
  b=c=3 as the unique fiber where B can exceed the endpoint product, while
  25B<27P holds sharply throughout that fibered universe. Whole-sphere
  27/25, five-copy, and (OS+) conclusions remain finite-exact only.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4219-no-sink-endpoint-energy-floor-and-near-ordinal-sharpness
  - THM-4223-cyclic-cut-cover-boolean-mobius-hierarchy-and-two-bad-owner-obstruction
  - THM-4224-order-ten-minimal-plus-min-two-bad-owner-obstruction
  - THM-4225-bad-owner-upper-zeta-successor-rook-hierarchy
related:
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
scripts:
  - 04-computation/tournament_two_reversal_sphere_audit_thm4239.py
  - 04-computation/tournament_two_reversal_osplus_audit_thm4239.py
  - 04-computation/tournament_two_reversal_words_thm4239.py
  - 04-computation/tournament_two_reversal_order18_driver_thm4239.sh
outputs:
  - 05-knowledge/results/tournament_two_reversal_sphere_audit_thm4239.out
  - 05-knowledge/results/tournament_two_reversal_osplus_audit_thm4239.out
  - 05-knowledge/results/tournament_two_reversal_order15_18_candidate_thm4239.out
script_sha256:
  - 6686ed4315a5343f25514d5541edf7d2444c68aa739c0c3f68457b78e52ad0db
  - 29110f8a94995081056143270723ee06588e63ab9e7a230ef97c5cdb52d25cde
  - 450504c9cb05720f389cf02758c32bbd80fc01eb10c5c483bbe2960063405238
  - 305cfc06a6633ca1589ed27c593848f4f04c87af5f13f482a155f5ecfcc67132
output_sha256:
  - 2b0715daee2504b088422961917929a3e5621ecfd2f46a1407006b8f97590915
  - f44d10527508eb69cc43080992d79755063be6a8489480936028690d2f015eb4
  - 632baa31dc2e9f045fbf33531215c5cbf0d2d29489ced55fe9ae1bee811ef829
hash_basis: raw LF bytes
audit: >
  Hostile audit PASS after two repairs: eta=h also at (b,c)=(2,1), whose
  negative constant term preserves the unique b=c=3 obstruction, and the
  T(n,1) profile begins at n=5. Independent reachability, permutation,
  formula, and coefficient-sign probes agree with the exact 47,413-owner,
  350,550-remainder, and order-18 censuses. Presentation equality is kept
  separate from tournament isomorphism, and every finite-only scope is
  explicit.
---

# THM-4239 -- fixed-gauge two-reversal sphere and tail fiber

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Status separation used throughout:

1. **PROVED, all orders:** the classification and count of strong labelled
   presentations.
2. **PROVED, all orders:** the formulas below only for the type-B owner fiber
   `(x_t,z)`.
3. **FINITE-EXACT:** every presentation and every owner pair through order 14;
   a narrower `25B<=27P` scan through order 18; and the stated `(OS+)` grid.

Nothing here is an isomorphism-class census.  No all-order conclusion about
type A, non-tail owner pairs, arbitrary strong tournaments, HYP-9081, or
general `(OS+)` is claimed.

## 1. Notation

Let `P_n` be the labelled transitive tournament on

```text
0<1<...<z,                 z=n-1,
```

with `i->j` iff `i<j`.  A presentation is an unordered pair of distinct base
arcs; reversing those two arcs gives `P_n^R`.  Presentation equality means
literal equality of labelled reversal pairs, not tournament isomorphism.

For a tournament `T`, write:

```text
H(T)       number of directed Hamilton paths;
End_i      number of directed Hamilton paths ending at i;
c(T)       number of directed Hamilton cycles modulo rotation;
N_ij       number of spanning two-directed-path covers ending at {i,j};
B_ij       N_ij-End_i-End_j+c(T);
P_ij       End_i End_j;
mu_ij      min(End_i,End_j);
rho_ij     d_i^- d_j^- - |N^-(i) intersection N^-(j)|.
```

Thus `B_ij` is the exact two-bad-owner mass of THM-4223/4224 and `rho_ij` is
THM-4225's normalized two-successor-rook count.

For the global statistic, let `D_i` be one-defect mass owned by `i`,
`D=sum_i D_i`, and

```text
W=(n-1)H+D,
energy=sum_i(H+D_i)^2,
target=(W+H)(W+3H),
G=target-5 energy
 =D^2+2(n-4)HD+n(n-3)H^2-5sum_i D_i^2.
```

## 2. All-order classification of strong presentations

For every `n>=3`, the strong presentations are the following disjoint union:

```text
A_n = { {(0,z),e} :
        e is a base arc other than (0,z),(0,1),(z-1,z) },

B_n = { {(0,b),(c,z)} : 1<=c<=b<=z-1 }.
```

Consequently

```text
|A_n|=binom(n,2)-3,
|B_n|=(n-2)(n-1)/2,
|A_n union B_n|=(n-1)^2-3=n^2-2n-2.
```

### Necessity

Strongness forces vertex `0` to acquire an inneighbor, hence some reversed
arc is `(0,b)`.  It forces `z` to acquire an outneighbor, hence some reversed
arc is `(c,z)`.

If one reversal performs both jobs, it is `(0,z)` and the other reversal is
`e`.  Reversing `(0,1)` makes `1` a source; reversing `(z-1,z)` makes `z-1`
a sink.  These and repetition of `(0,z)` are exactly the three exclusions.

Otherwise the two reversals are `(0,b),(c,z)` with `b<z` and `c>0`.  If
`b<c`, the cut `{0,...,b}` dominates its complement: both reversals are
internal to one side, so no edge enters the cut.  Strongness therefore
requires `c<=b`.

### Sufficiency

In type A, if `e` is not an adjacent base arc, the directed Hamilton cycle

```text
0,1,...,z,0
```

survives.  If `e=(a,a+1)` with `1<=a<=z-2`, it is repaired to

```text
0,...,a-1,a+1,a,a+2,...,z,0.
```

In type B with `c<=b`, a directed Hamilton cycle is

```text
0,1,...,c-1,b+1,b+2,...,z,c,c+1,...,b,0,
```

with empty intervals omitted.  A tournament containing a directed Hamilton
cycle is strong.  This completes the all-order proof.

## 3. All-order type-B `(x_t,z)` tail formulas

Write the type-B presentation as `Y(b,c,m)` on

```text
0<1<...<b<x_1<...<x_m<z,
```

with reversals `(0,b),(c,z)`, where `1<=c<=b` and `m>=1`.  Its order is
`n=b+m+2`.  The strong type-B boundary `m=0` exists but has no tail owner and
is outside this fiber theorem.

Let `C_b=Y(b,c,m)[{0,...,b}]`; it is independent of `c,m`.  Define

```text
h_1=1,
h_b=2^(b-1)+1                                      (b>=2),

K_b=2*3^(b-1)+2^(b+1)-2                           (b>=1),

gamma_(b,c)=2^max(b-c-1,0).
```

For `b>=2` and `1<=d<=b`, put

```text
eta_b(d)=h_b                                         if d=b,
eta_b(d)=2^d 3^(b-d-1)+2^(b-d)-1_[d=1]             if d<b.
```

Then set

```text
eta=eta_b(c),
lambda=2^b-1                                         if c=b,
lambda=eta_b(b-c)                                    if c<b.
```

At the exceptional two-vertex core `b=c=1`, set

```text
(h,K,eta,gamma,lambda)=(1,4,2,1,1).
```

For every `m>=1`, the full path and cycle counts are

```text
H(Y)=eta*2^m+h,                    c(Y)=gamma,
End_z=h.
```

For `1<=t<=m`, put `u=2^(t-1)` and

```text
epsilon_t=0          if t<m,
epsilon_t=lambda     if t=m.
```

The complete `(x_t,z)` profile is

```text
End_(x_t)=eta*u,
N_(x_t,z)=K*u-epsilon_t,
B_(x_t,z)=(K-eta)u-epsilon_t-h+gamma,
P_(x_t,z)=eta*h*u,
mu_(x_t,z)=h,
rho_(x_t,z)=(b+t)(b+m)-(b+t-1).
```

### Derivation

For a subset of `C_b` not containing both `0,b`, the induced tournament is
transitive and has one Hamilton path.  If it contains both endpoints and `r`
of the `b-1` middle vertices, its Hamilton paths consist of the paths using
the block `b,0` (one for each split of those `r` vertices) and, when `r>0`,
the increasing path.  Hence `H(C_b)=h_b`.

Partitioning the core between two labelled paths gives

```text
K_b=sum_(A subseteq C_b) H(C_b[A])H(C_b[C_b\A]).
```

If `0,b` are split, the contribution is `2^b`.  If they lie together, the
two choices of component contribute
`2*(3^(b-1)+2^(b-1)-1)`.  Their sum is the displayed `K_b`.

The remaining constants have literal core meanings

```text
eta_(b,c)=sum_A H(C_b[A]) Start_c(C_b[C_b\A]),
lambda_(b,c)=sum_A H(C_b[A]) End_c(C_b[C_b\A]).
```

A disjoint increasing-run enumeration at the forced start `c` gives, for
`c<b`, contributions

```text
2^(b-c)  and  2^c 3^(b-c-1)-1_[c=1],
```

and gives `h_b` at `c=b`.  Reversing a core path and reflecting labels
`r |-> b-r` gives `lambda_(b,c)=eta_b(b-c)` for `c<b`; direct evaluation at
`c=b` gives `2^b-1`.  Rotating a cycle immediately after `z` and making the
same increasing-run split gives `gamma_(b,c)`.

For a Hamilton path ending at `x_t`, every later tail vertex is forced before
the block `z,c`; each earlier tail vertex has a binary side choice.  This
gives `eta*2^(t-1)`.  For a two-path cover ending at `(x_t,z)`, every later
tail vertex is forced into the `z` component and the core partition sum is
`K_b`.  If `t<m`, that component has a later-tail vertex and no boundary
exclusion.  If `t=m`, the covers whose bare core path ends at `c` expose the
reversed edge `c,z`; exactly `lambda` covers must be removed.  Substitution
in `B=N-End_(x_t)-End_z+c(Y)` gives the formulas.

### Exact local boundaries inside this fiber

Let

```text
Q_(b,c)=K_b-(h_b+1)eta_(b,c).
```

For `b>=2`, `eta_(b,c)>=h_b`, with equality if and only if `c=b` or
`(b,c)=(2,1)`, so

```text
Q_(b,c) <= q_b
q_b=2*3^(b-1)+2^(b-1)-4^(b-1)-4.
```

Indeed, if `c<b` and `k=b-c`, then

```text
eta_(b,c)-h_b
 =2^c(3^(k-1)-2^(k-1))+2^k-1_[c=1]-1.
```

For `k=1` this is `1-1_[c=1]`, and for `k>=2` it is positive.  This proves
the corrected equality statement, including its exceptional `(2,1)` row.

Here `q_2=0`, `q_3=2`, and `q_b<0` for every `b>=4`; the `b=1` coefficient
is zero.  The base `b=4` is `-6`, and the inequality persists by the elementary
induction on
`4^r-2*3^r-2^r+4` (`r=b-1`).  The extra equality `(b,c)=(2,1)` is harmless:
for `b=1`, direct substitution gives `B-P=-epsilon_t<=0`; for `b=2`,
`Q=0` but `B-P=-2-epsilon_t<0`; and for `b=3,c<3`, integrality and
`eta>h` give `Q<=-4`.  Therefore the only type-B tail fiber in which `B>P`
is possible is `b=c=3`, namely the already-known `X_m` family.

Specializing the general formula recovers THM-4224's inherited boundary (it
is not a new proof claim about owner pairs outside this fiber):

```text
interior t<m: B-P=2^t-4,              B-P-mu=2^t-9;
terminal t=m: B-P=2^m-11,             B-P-mu=2^m-16.
```

Thus, inside the fiber, product failures occur exactly at interior `t>=3`
and terminal `m>=4`; `+min` failures occur exactly at interior `t>=4` and
terminal `m>=5`.  For fixed order `n=m+5`, the counts are `n-7` for `n>=9`
and `n-8` for `n>=10`, respectively (zero below those ranges).

For the multiplicative diagnostic,

```text
27P-25B = R_(b,c)u+25(epsilon_t+h-gamma),
R_(b,c)=(27h+25)eta-25K.
```

At `b=1`, `R=4`.  For `b>=2`, it is bounded below by

```text
r_b=27*4^(b-1)-21*2^(b-1)+102-50*3^(b-1).
```

The values at `b=2,3,4` are `18,0,312`; for `b>=4`, positivity follows from

```text
r_(b+1)-3r_b=27*4^(b-1)+21*2^(b-1)-204>0.
```

Equality in the slope occurs only at `b=c=3`, where the remaining slack is
exactly `100` in an interior fiber and `275` at the terminal.  Hence

```text
25B_(x_t,z) < 27P_(x_t,z)
```

for every type-B tail row in every order, with `27/25` the sharp limiting
constant along `X_m`.

Finally, on `X_m`, the exact rook margins are

```text
interior:
 P+rho-B     =(t+3)m+2t+11-2^t,
 P+rho+mu-B  =(t+3)m+2t+16-2^t;

terminal:
 P+rho-B     =m^2+5m+18-2^m,
 P+rho+mu-B  =m^2+5m+23-2^m.
```

The terminal margins are nonnegative exactly for `m<=6` and both first fail
at `m=7` (order 12).  Interior failures are exactly the integer pairs
satisfying the corresponding displayed strict exponential inequality.

## 4. Full-sphere finite exact results

Literal reachability and subset DP were run for every two-reversal candidate,
every classified strong presentation, and every unordered owner pair for
`3<=n<=14`:

```text
candidate presentations tested: 13,104
strong presentations profiled:      782
all-owner pair rows:              47,413
type-B tail rows matching every symbolic formula above: 1,001
```

The exact rows by order are frozen in
`05-knowledge/results/tournament_two_reversal_sphere_audit_thm4239.out`.
Across this finite universe:

```text
n                 3  4  5  6  7  8  9 10 11 12 13 14
product failures  0  0  0  0  0  0  2  3  4  5  6  7
+min failures     0  0  0  0  0  0  0  2  3  4  5  6
rho failures      0  0  0  0  0  0  0  0  0  1  2  3
25B>27P           0  0  0  0  0  0  0  0  0  0  0  0
G<0               0  0  0  0  0  0  0  0  0  0  0  0
```

Every local failure in this **finite** range is an `X_m` `(x_t,z)` row.  This
last sentence is not promoted to an all-order statement about all owner pairs.

The maximum `B/P` is `0` at `n=3`, is `1` for `4<=n<=8`, and for `9<=n<=14`
is uniquely

```text
(27u-4)/(25u),          u=2^(n-7),
```

at presentation `{(0,3),(3,z)}` and owners `(n-3,z)`.  An optimized exact
subset-DP scan extends only `25B<=27P` and this maximum through `n=18`, with

```text
n=15: 6908/6400,  n=16: 13820/12800,
n=17: 27644/25600, n=18: 55292/51200,
```

and zero candidate failures on another 125,071 owner-pair rows.

### Five-copy minima (finite sphere)

The exact minimum gaps for `n=3,...,14` are

```text
0, 44, 510, 3186, 19842, 122778, 753810, 4597146,
27894642, 168653178, 1017163410, 6124242906.
```

Equality occurs only at the labelled `C3`.  For `n>=6` in this finite range,
the minimum has presentation multiplicity three:

```text
{(0,z-2),(z-2,z)}, {(0,z-1),(z-2,z)}, {(0,z),(z-2,z)}.
```

These are three labelled presentations of one tournament isomorphism type,
THM-4219's `T(n,1)`.  At `n=5` the same minimum has five presentations; at
`n=4` all six presentations represent the unique strong class and have gap
`44`.

For `n>=5`, the minimum candidate itself has the inherited all-order
profile below; this is not an all-order sphere-minimality assertion.  Set
`x=2^(n-4)` and `y=3^(n-4)`:

```text
H=3x+3,
W=14y+18x-12,
energy=(609x^2+570xy-330x+196y^2-180y+45)/5,
target=(W+H)(W+3H),
G=6(17xy+2y-7x^2+4x-3)>0.
```

The transitive order-seven tournament is a filter-hostile control: it is not
strong and its exact gap is `-10922`.

## 5. Finite `(OS+)` boundary

Using THM-4184's exact capacity/remainder engine,

```text
R_+(A,C)=G_+(A triangleright C)-H(C)^2G_+(A)-H(A)^2G_+(C),
```

the audit tested all 475 strong two-reversal presentations of orders `3..12`
against every labelled no-sink right tournament of orders `3..5`:

```text
right counts: q=3:2, q=4:32, q=5:704, total 738
exact remainders: 475*738=350,550
nonpositive remainders: 0.
```

The minimum remainders by left order are

```text
3204, 27336, 204996, 1497060, 10355052,
70305924, 477746076, 3385587984, 23959160964, 179696428164.
```

In every order a labelled `C3` is a minimizing right factor.  Dropping the
no-sink condition is a genuine hostile boundary: for left `C3` and transitive
right factors of orders `1,...,5`, the exact remainders are

```text
-72, -288, -756, -1800, -4104.
```

This is FINITE-EXACT only; it does not prove `(OS+)` for arbitrary right order
or for all-order two-reversal left factors.

## 6. Replay and artifacts

From the repository root, reproduce the two Python censuses at both interpreter
optimization levels and compare them with the frozen outputs:

```bash
python3 -B 04-computation/tournament_two_reversal_sphere_audit_thm4239.py \
  > /tmp/tournament-two-reversal-sphere.out
python3 -B -O 04-computation/tournament_two_reversal_sphere_audit_thm4239.py \
  > /tmp/tournament-two-reversal-sphere-opt.out
cmp /tmp/tournament-two-reversal-sphere.out \
  /tmp/tournament-two-reversal-sphere-opt.out
cmp /tmp/tournament-two-reversal-sphere.out \
  05-knowledge/results/tournament_two_reversal_sphere_audit_thm4239.out

python3 -B 04-computation/tournament_two_reversal_osplus_audit_thm4239.py \
  > /tmp/tournament-two-reversal-osplus.out
python3 -B -O 04-computation/tournament_two_reversal_osplus_audit_thm4239.py \
  > /tmp/tournament-two-reversal-osplus-opt.out
cmp /tmp/tournament-two-reversal-osplus.out \
  /tmp/tournament-two-reversal-osplus-opt.out
cmp /tmp/tournament-two-reversal-osplus.out \
  05-knowledge/results/tournament_two_reversal_osplus_audit_thm4239.out

sh 04-computation/tournament_two_reversal_order18_driver_thm4239.sh "$PWD" \
  > /tmp/tournament-two-reversal-order15-18.out
cmp /tmp/tournament-two-reversal-order15-18.out \
  05-knowledge/results/tournament_two_reversal_order15_18_candidate_thm4239.out
```

The stable C++ driver hash-locks THM-4224's inherited subset-DP source, raises
only its census order cap, corrects the inherited presentation label, and
feeds it every word produced by the independent presentation generator.
The checked-in artifacts are:

```text
04-computation/tournament_two_reversal_sphere_audit_thm4239.py
04-computation/tournament_two_reversal_osplus_audit_thm4239.py
04-computation/tournament_two_reversal_words_thm4239.py
04-computation/tournament_two_reversal_order18_driver_thm4239.sh
05-knowledge/results/tournament_two_reversal_sphere_audit_thm4239.out
05-knowledge/results/tournament_two_reversal_osplus_audit_thm4239.out
05-knowledge/results/tournament_two_reversal_order15_18_candidate_thm4239.out
```

Normal and optimized Python runs byte-match.  The full-sphere output explicitly
freezes the three minimizing presentations for every `n>=6` and every
`rho`-failure row, rather than relying on printed counts alone.  The
independent hostile audit additionally rebuilt reachability through order
eight, every core constant through `b=8`, all 70 tail rows through order
eight, and scanned the coefficient signs through `b=200`.

## 7. Explicit nonclaims

- The strong-presentation count is a count of distinct labelled tournaments
  in this one fixed transitive gauge.  It is not invariant under relabelling,
  is not an isomorphism-class count, and does not describe a radius-two ball
  around every transitive presentation.
- The all-order local formulas and inequalities cover only type-B owner pairs
  `(x_t,z)` with `m>=1`.
- The assertion that every full-sphere product failure is an `X_m` row is
  only FINITE-EXACT through order 14.
- `25B<=27P` on the whole sphere is FINITE-EXACT through order 18, not proved
  for arbitrary order.
- Sphere-wide HYP-9081 positivity/minimality is FINITE-EXACT through order 14.
- The `(OS+)` statement is only the displayed finite left/right grid.
- None of these finite statements applies to arbitrary strong tournaments.
