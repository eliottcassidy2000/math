# Synchronized half-grid physical covers are a finite divisor-chart problem

**Date:** 2026-08-15  
**Status:** PROVED-ELEMENTARY all-`q` affine-chart reduction for synchronized
half-grid physical times; PROVED-ELEMENTARY exact half-grid rank
families on every even `q>=8`, every odd multiple of three `q>=9`, and every
odd multiple of five not divisible by three `q>=25`; PROVED-ELEMENTARY exact
rank table for `15<=q<=28`, with VERIFIED-EXACT independent set-cover and
normalization audits; PROVED-ELEMENTARY exact even-rank-two leakage ladder
`P=-a q^2/2` for odd `1<=a<q/7`.  MISTAKE-389 retracts the original zero-mode-cochain
interpretation: proved [THM-3405](../01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md)
supplies an additional mode-divisibility condition.  The exact half-grid
results survive and have no LRC(14) ledger consequence.

## 1. Inheritance and the corrected question

[THM-3398](../01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md)
turns a physical cover into a selected-mode block cover carrying a complete
affine cochain.  [THM-3402](../01-canon/theorems/THM-3402-atomized-sheet-covers-and-constructive-cochain-locus.md)
independently refines this into repeated-owner atoms, reconstructs the
complete physical cover locus, and compresses exactly back to one consecutive
mode per owner.  THM-3405 then proves that a vanishing **mode** cochain has an
extra owner-gcd divisibility law.  The present result classifies the larger
physical slice where every owner phase lies on the common half-grid; it does
not assert that the chosen physical time is any selected mode's centre.
[THM-3401](../01-canon/theorems/THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md)
computes the rank at fixed source centre zero.  MISTAKE-384 separates that
slice from the genuine zero-mode-cochain locus; MISTAKE-389 in turn separates
the genuine zero-mode-cochain locus from the still larger half-grid physical
slice.

For a positive transverse owner `u`, define its exact danger set at physical
time `c` by

```text
D_(q,u)(c)={ell in Z/qZ: ||u(c+ell/q)||<1/14}.          (1)
```

Let `r_hg^+(q)` be the least number of **distinct positive owners**, with no
upper bound on their sizes, whose exact sets `(1)` cover `Z/qZ` at one time
`c` satisfying

```text
2quc in Z                  for every active owner u.   (1a)
```

This is the synchronized half-grid rank.  Every zero-mode-cochain certificate
obeys `(1a)`, but the converse is false.

The inheritance pass is:

- closest proved mechanism: THM-3398's mode/cochain iff and quotient-block
  formula;
- canonical hostile: the `q=15,c=1/150` physical partition whose half-grid
  scalar fails THM-3405's mode divisibility;
- corrected near misses: MISTAKE-384's lost additive gauge and MISTAKE-389's
  half-grid/mode-centre conflation;
- least-used sidecar: `gcd(q,u)|h` for each selected mode residue `h`.

The connection contract for the new reduction is:

| field | content |
|---|---|
| source | rational physical times satisfying `(1a)` and all positive transverse owners |
| target | finitely many divisor-indexed set-cover charts `O(q,g),E(q,g)` |
| map | reduce `theta=qc=a/b`, divide the forced factor of each owner, then apply one affine sheet permutation and one unit owner reindexing |
| preserved | exact strict danger sets, half-grid cover rank, transversality, physical realization, and common sheet-affine gauge |
| destroyed | literal time numerator/denominator, literal owner size, chosen sheet origin, and mode-centre data |
| required sidecar | odd/even denominator type and `g=gcd(floor(b/2),q)`; THM-3405 mode divisibility is additionally required to infer zero cochain |

## 2. The infinite centre locus collapses to finite charts

Write the common twist in lowest terms as

```text
theta=qc=a/b,                 gcd(a,b)=1.              (2)
```

By the synchronized half-grid definition `(1a)`,

```text
2quc in Z.
```

Equation `(2)` therefore forces

```text
b divides 2u.                                           (3)
```

This elementary divisibility is the finite-state coordinate for the physical
half-grid slice.  It is only a necessary coordinate for a common mode centre.

### Odd denominator

If `b` is odd, write `u=bv`.  Then `(1)` is decided exactly by

```text
v(a+b ell) mod q.                                      (4)
```

Put `g=gcd(b,q)`.  Transversality makes `g<q`, and `g` is odd.  A sheet
translation can replace `a` by `A=a+bk` coprime to `q`: primes dividing `g`
already miss `a`, and for every other prime divisor of `q` one avoids one
residue class of `k`.  Write `b=g b_1` and `q=gq_1`.  Since `b_1` is a unit
modulo `q_1`, choose a unit sheet multiplier `s mod q` satisfying

```text
b_1 s = A mod q_1.                                    (5)
```

After `ell -> k+s ell` and the unit owner reindex `v -> Av`, `(4)` becomes

```text
v(1+g ell) mod q.                                     (6)
```

This is the odd chart `O(q,g)`.  Its owner types are `v mod q` with
`q/g` not dividing `v`.

### Even denominator

If `b=2d`, then `a` is odd and `(3)` gives `u=dv`.  Exact danger is decided by

```text
v(a+2d ell) mod 2q.                                   (7)
```

Put `g=gcd(d,q)<q`.  The same prime-avoidance argument chooses
`A=a+2dk` coprime to `2q`, and a unit `s mod q` with

```text
(d/g)s=A mod q/g.                                     (8)
```

The same affine sheet change and owner reindex turn `(7)` into

```text
v(1+2g ell) mod 2q.                                   (9)
```

This is the even chart `E(q,g)`, with `q/g` not dividing `v`.

Conversely, every listed chart is physical: take

```text
O(q,g): c=1/(gq),   u=gv,   g an odd proper divisor of q;
E(q,g): c=1/(2gq), u=gv,   g any proper divisor of q. (10)
```

Equations `(6)` and `(9)` are then the literal phases in `(1)`.  Thus the
apparently infinite rational-centre problem has exactly

```text
#{g:g|q,g<q,g odd} + #{g:g|q,g<q}                     (11)
```

chart types.

### The missing mode-centre divisibility

The chart theorem stops at synchronized physical time.  THM-3405 writes the
active owner gcd as `d_U`, puts `g_U=gcd(q,d_U)`, and sets

```text
a_U=2q d_U c.
```

A zero mode cochain requires `g_U|a_U`.  This does not follow from `(3)`.
At `q=15`, chart `E(15,5)` has the exact partition

```text
c=1/150,   U=(5,40,50),   d_U=g_U=5,   a_U=1.        (11a)
```

All three owners satisfy the half-grid condition and their danger sets are
the three residue classes modulo three, but `5` does not divide `1`.
Equivalently, owner five has half-grid residue `h=1`, while a lawful mode
centre residue must be divisible by `gcd(15,5)=5`.  Thus `(11a)` is a
physical partition with nonzero relative mode-centre data, not a mobile
common-mode-centre certificate.

The containing mode is unique for each owner.  In owner order `(5,40,50)`,
its centres and THM-3398 pair cochain are

```text
(x_5,x_40,x_50)=(0,1/120,1/150),
(P_(5,40),P_(5,50),P_(40,50))=(-50,-50,100),          (11b)
```

where `P_ij=2q u_i u_j(x_i-x_j)`.  Triangle closure holds exactly and
`(||P||_1,||P||_infinity)=(200,100)`.  This turns the typing correction into
a first explicit cochain-leakage control: half-grid synchronization does not
make the drift small, let alone zero.

### Owner replicas make dominance safe

In `O(q,g)` a danger set depends only on `v mod q`; in `E(q,g)` it depends
only on `v mod 2q`.  Every type therefore has infinitely many positive owner
replicas.  If one danger set is contained in another, it may be replaced by
the larger set without violating owner distinctness: use a fresh period
replica if necessary.  Hence deduplication and inclusion-maximalization
preserve the minimum all-owner rank.  This point would be false under the
literal owner cap `{1,...,14}`.

## 3. Universal capacity and three infinite exact rank families

For any transverse owner, put

```text
m=q/gcd(u,q)>1.
```

Its danger set is the inverse image of a strict arc on an `m`-point cyclic
grid.  Such an arc contains at most `ceil(m/7)` grid points, so every block
has at most

```text
B(q)=max_(m|q,m>1) (q/m) ceil(m/7)                    (12)
```

sheets.  Consequently

```text
r_hg^+(q) >= ceil(q/B(q)).                             (13)
```

The following constructions attain `(13)` on three infinite divisor
families.  Put `q=pd`, choose the chart `E(q,d)`, use

```text
c=1/(2dq),                  u=dv.                     (14)
```

The displayed owners partition the sheets into the `p` residue classes.

### Even sheets

For `p=2`, `q>=8`, take

```text
V_2(d)={1,2d-1}.                                      (15)
```

The two danger sets are precisely the even and odd sheets.  Here `B(q)=q/2`,
so

```text
r_hg^+(q)=2                     (q even, q>=8).        (16)
```

### The same parity partition carries a quadratic leakage ladder

The chart quotient in `(16)` retains the even/odd blocks but forgets an
unbounded exact mode sidecar.  Let `q=2d>=8` and choose any odd integer `a`
with `1<=a<q/7`.  At

```text
c_a=a/q^2,             U_q=(q/2,q(q-1)/2),            (16a)
```

the first owner has phase `a/(2q)+ell/2`, while the second has phase
`(q-1)(a/q+ell)/2`.  The strict inequality `7a<q` therefore makes their
danger sets exactly the even and odd sheets.  Thus every `a` in `(16a)`
realizes the **same** rank-two Boolean partition.

Their unique containing THM-3398 single-atom modes have centres

```text
x_0=0,                 x_1=a/[q(q-1)].                (16b)
```

Indeed `c_a` lies within the first source radius precisely when `7a<q`, and

```text
|c_a-x_1|=a/[q^2(q-1)] < 1/[7q(q-1)]
```

is the corresponding strict inequality for the second owner.  Substitution
in the exact pair cochain gives the closed formula

```text
P_01=2q u_0 u_1(x_0-x_1)=-a q^2/2.                  (16c)
```

This is a genuine leakage tariff, not a zero-cochain representative.  The
active owner gcd and THM-3405 gauge are

```text
d_U=gcd(U_q)=q/2,      g_U=gcd(q,d_U)=q/2,
a_U=2q d_U c_a=a.                                      (16d)
```

Since `0<a<q/7<q/2`, the necessary divisibility `g_U|a_U` fails for every
member of the ladder.  Hence one fixed two-block clutter supports linearly
many inequivalent physical scalars and cochain magnitudes.  XOR, an
unweighted two-vertex tournament, or the normalized divisor chart sees only
the parity split and erases exactly the integer `a` recovered by `(16c)`.

The tariff is quadratically homogeneous in sheet degree: replacing `q` by
`kq` at fixed `a` multiplies `|P_01|` by `k^2`.  This agrees with the
degree-two cochain homogeneity in THM-3398.  It is an agreement of the exact
sidecar law; the canonical witnesses `(16a)` themselves are not being
identified with literal pure dilates.

There are `q/14+O(1)` admissible odd scalars at each even `q`.  Consequently,
if `A(N)` counts these positive leakage witnesses through `q<=N`, then

```text
A(N)=N^2/56+O(N),
sum_(q<=N,q even) #{a:1<=a<q/7,a odd}/q
    =N/28+O(log N).                                   (16e)
```

This is the harmonic-series distinction made concrete.  The support set of
sheet degrees in `(16)` has logarithmic reciprocal mass, but restoring the
discarded scalar fibre turns it into a two-dimensional lattice whose
reciprocal-degree mass is linear.

### Odd multiples of three

For `p=3`, the boundary `q=9` uses `V_3(3)={1,5,7}`.  For odd `d>=5`, take

```text
V_3(d)=
 {1,2d-1,2d-2},  d=0 mod 3;
 {1,2d,  2d-1},  d=1 mod 3;
 {1,2d-2,2d  },  d=2 mod 3.                           (17)
```

Substitution in `(9)` gives one complete residue class modulo three per
owner.  For every odd divisor `m>=3`,
`ceil(m/7)/m<=1/3`, so `B(q)=q/3`.  Hence

```text
r_hg^+(q)=3        (q odd, 3|q, q>=9).                (18)
```

### Odd multiples of five away from three

Let `p=5`, with `q` odd, `3` not dividing `q`, and `q>=25`.  Set

```text
V_5(d)={1} union ({2d-2,...,2d+2} \ 5Z).              (19)
```

The interval in `(19)` contains one representative of every class modulo
five; removing the multiple of five leaves four unit classes.  The owner
`v=1` fires exactly class zero.  For a remaining `v`, the unique solution of
`v ell=4 mod 5` is dangerous with residue distance `|v-2d|<=2`, while all
other classes have distance greater than `5d/7`.  Thus the five blocks are
the five residue classes.  If `m|q,m>1`, then `m` is odd, not divisible by
three, and

```text
ceil(m/7)/m <= 1/5
```

(check `m=5,7,11,13`; the inequality is immediate for `m>=15`).  Therefore
`B(q)=q/5` and

```text
r_hg^+(q)=5        (q odd, 5|q, 3 not|q, q>=25).      (20)
```

These are divisor pullbacks rather than mere copies of the fixed-zero
construction.  Multiplying `d` multiplies every residue-class fibre degree,
scales owner speeds linearly, and scales the canonical centre quadratically.
The partition is functorial under reduction modulo `p`, while literal owner
labels are not.  This is the precise degree-graded monoid shadow: fibre
degrees compose multiplicatively, but an owner-word sidecar is needed inside
each grade.

## 4. Exact all-owner table from fifteen through twenty-eight

The capacity bound `(13)` is sharp throughout the window except at `q=17,19`.
For those two primes, the only charts are `O(q,1)` and `E(q,1)`.  In either
chart all maximal danger sets have size three, share one central sheet, and
their noncentral pairs are mutually disjoint.  Covering all noncentral sheets
therefore requires every pair, namely `(q-1)/2` owners.  This proves the two
remaining lower bounds without search.

| `q` | `r_hg^+(q)` | minimizing chart(s) | one exact owner witness |
|---:|---:|---|---|
| 15 | 3 | `O(15,5), E(15,5)` | `(5,20,25)` |
| 16 | 2 | `E(16,4), E(16,8)` | `(8,56)` |
| 17 | 8 | `O(17,1), E(17,1)` | `(1,2,3,4,5,6,7,8)` |
| 18 | 2 | `O(18,9), E(18,9)` | `(9,81)` |
| 19 | 9 | `O(19,1), E(19,1)` | `(1,2,3,4,5,6,7,8,9)` |
| 20 | 2 | `E(20,5), E(20,10)` | `(10,90)` |
| 21 | 3 | `O(21,7), E(21,7)` | `(7,49,56)` |
| 22 | 2 | `O(22,11), E(22,11)` | `(11,121)` |
| 23 | 6 | `E(23,1)` | `(1,4,5,7,9,11)` |
| 24 | 2 | `E(24,4), E(24,6), E(24,12)` | `(12,84)` |
| 25 | 5 | `O(25,5), E(25,5)` | `(5,20,30,45,55)` |
| 26 | 2 | `O(26,13), E(26,13)` | `(13,169)` |
| 27 | 3 | `E(27,3), O(27,9), E(27,9)` | `(9,45,63)` |
| 28 | 2 | `E(28,7), E(28,14)` | `(14,182)` |

The sequence is therefore

```text
3,2,8,2,9,2,3,2,6,2,5,2,3,2.                        (21)
```

The original draft compared `(21)` with the owner-14 mobile common-mode-centre
sequence.  MISTAKE-389 withdraws that comparison: the two ranks optimize
different predicates.  Since every zero-mode-cochain cover is half-grid,
`(21)` supplies lower bounds for the genuine all-owner ranks, but no equality
follows without imposing THM-3405's gcd gauge on every candidate.  That
separate computation is now complete in
[the divisor-ancestry reflection](lrc-zero-mode-cochain-divisor-ancestry-and-boolean-realization-codex-20260815.md):
the genuine sequence is `(6,4,8,4,9,5,8,6,6,4,6,7,4,7)`.

## 5. The periodic harmonic support of the infinite theorem

The three families `(16),(18),(20)` are disjoint and, after their finite
boundaries, depend only on `q mod 30`:

```text
S = {q even}
    union {q odd:3|q}
    union {q odd:5|q,3 not|q}.                         (25)
```

They occupy respectively `15`, `5`, and `2` residue classes modulo 30.  Thus
the exact synchronized half-grid rank is now known on `22/30=11/15` of all
sheet numbers, and
[THM-3359](../01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md)
gives

```text
#(S intersect [1,N])=(11/15)N+O(1),
sum_(q<=N,q in S) 1/q=(11/15)log N+O(1).              (26)
```

Weighting by the proved ranks gives the logarithmic coefficient

```text
2(1/2)+3(1/6)+5(1/15)=11/6.                           (27)
```

This is a genuine subset-of-the-harmonic-series statement.  It is not a
claim about the unknown complementary `4/15` of sheet numbers.

## 6. Pullback to the Berggren spine and Fibonacci recurrence

[THM-3334](../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
identifies the distinguished Berggren `U`-spine value

```text
Q_n=(2n+3)^2+2=4n^2+12n+11.                           (28)
```

Modulo three,

```text
Q_n=n^2+2 mod 3,
```

so `3|Q_n` exactly when `3` does not divide `n`.  Every `Q_n` is odd.  Also
`Q_n=4n^2+2n+1 mod 5`, whose discriminant is `3 mod 5`, a nonsquare; hence
`5` never divides `Q_n`.  It follows from `(18)` that

```text
r_hg^+(Q_n)=3                   whenever n not=0 mod 3. (29)
```

Thus the synchronized half-grid rank is exactly three on two thirds of the
**parameter classes** along the hidden Berggren branch.  The value subseries
`sum 1/Q_n` nevertheless converges because `Q_n` is quadratic.  This is the
index/value distinction emphasized by THM-3359: periodic parameter support
does not force a divergent harmonic series on sparse represented values.

For Fibonacci numbers `F_n`, the elementary modular recurrences give

```text
2|F_n iff 3|n,      3|F_n iff 4|n,      5|F_n iff 5|n. (30)
```

After the finite small-value boundaries, `(16),(18),(20)` therefore imply

```text
r_hg^+(F_n)=2  if 3|n;
r_hg^+(F_n)=3  if 4|n and 3 not|n;
r_hg^+(F_n)=5  if 5|n, 3 not|n, and 4 not|n.          (31)
```

On the Pisano-60 index cycle these occupy `20+10+6=36` classes, so the exact
rank is classified on an eventual density `3/5` of Fibonacci indices.  The
index harmonic support has coefficient `3/5`, while
`sum 1/F_n` over the corresponding **values** converges exponentially.  A
linear recurrence, its modular index language, and its represented-value
harmonic subset are three different carriers.

## 7. Why this is not a tournament or XOR theorem

The affine charts retain physical danger blocks but forget the selected mode
centres, so they do not determine any complete pair observable.  MISTAKE-389's
q15 witness proves in particular that the observable need not be zero.  The
intrinsic retained object is a cyclic block partition or cover clutter, not a
tournament.  In particular:

- the two-owner family is a parity bipartition;
- the three-owner family is a ternary residue partition, not a directed
  triangle;
- the five-owner family is a five-colour partition;
- owner counts do not acquire orientations merely because their cardinalities
  are four or six;
- an orientation would require a separately reconstructed antisymmetric mode
  cochain, which the chart quotient destroys.

The ternary family does support a rooted branch-transplant picture: increasing
the fibre degree `d` pulls the same three quotient colours back to larger
fibres.  But this is a divisor ray with a residue-dependent owner sidecar,
not the free ternary ancestry tree of Berggren words.

The useful analogy with Rybin--Zhang--Luo's
[*XX^t Can Be Faster*](https://arxiv.org/abs/2505.09814) is algorithmic and
typed.  Their structured symmetric target admits a candidate-generation,
exact-relation, and minimum-cover pipeline.  Here arithmetic normalization
generates exact divisor charts and a minimum set-cover stage selects blocks.
RXTX's linear-span relation is absent here; conversely, LRC's strict phase,
owner, and centre gauge are absent there.  Reducing either problem to XOR
would additionally collapse signs in characteristic two.  The common grammar
is **structured target plus exact sidecar plus minimum cover**, not a theorem
transfer.

## 8. Verification and next frontiers

Run

```text
python 04-computation/lrc_all_owner_halfgrid_physical_chart_probe_20260815.py
python -O 04-computation/lrc_all_owner_halfgrid_physical_chart_probe_20260815.py
```

The companion:

1. pins THM-3398, THM-3402, and THM-3405;
2. audits every admissible affine numerator/slope and every transverse owner
   row for `15<=q<=28` against its normalized chart;
3. solves every chart by branch-and-bound and independently by exhaustive
   combinations;
4. checks the universal capacity lower bounds and the `q=17,19` central-pair
   anatomy, plus the q15 mode-divisibility hostile from MISTAKE-389;
5. checks every explicit divisor partition through `q=500` (`247` even,
   `82` ternary, and `32` five-colour instances);
6. reconstructs all `4,482` strict even-parity leakage instances through
   `q=500`, including their centres, gauge failure, quadratic cochain, and
   first excluded odd scalar;
7. contains no floating literals or `assert`, and pins a semantic digest that
   survives optimized execution.

The cheapest decisive continuations are:

1. classify the complementary odd `q` coprime to 30 by the prime-chart cover
   function, beginning with the hostile contrast `q=23,29,31`;
2. prove a composition law for general `E(q,g)` and `O(q,g)` chart ranks, not
   only the residue partitions of orders two, three, and five;
3. extend the now-complete THM-3405 divisor-ancestry computation beyond
   `Q=28` and prove or refute its apparent primitive rank-four floor;
4. classify the full nonzero cochain spectra over the ternary and five-colour
   half-grid minima; the even branch is now the exact ladder `(16a)--(16d)`;
5. transport a chart certificate into an actual refined LRC row.  Until that
   final typed map exists, none of these rank collapses decrements the live
   LRC(14) ledger.
