# Reflected-pair transfer to projected `k=2,3`: exact obstruction and the affine-typed survivor

Status: **PROVED scoped transfer + FINITE-EXACT + VERIFIED-EXACT**.
[THM-3360](../01-canon/theorems/THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails.md)
supplies the reflected high-pair floor for the `k=2` composition;
[THM-3382](../01-canon/theorems/THM-3382-independent-superunit-affine-tail-and-reflected-residue-closure.md)
is an independent weaker-floor route. Direct transfer to the projected wall
is **REFUTED** before projection. This sidecar changes no live ledger and is
not an arbitrary projected-sector closure.

## 1. The two observable types

Fix a six-body root `E`, ruler `L=14*lcm(E)`, and an integer body-safe cell
`j`.  For an arbitrary projected drift `z`, put

```text
A_z(j)={u in [0,1]: ||z(j+u)/L||<1/14},
S_j(Z)=[0,1] minus union_(z in Z) A_z(j).
```

Modulo the finitely many tooth endpoints, THM-2941 `(25g)--(25i)` says

```text
P_(E,Z)=phi_L(C_E minus union_(z in Z)D_z)=union_(j in J_E) S_j(Z).   (1)
```

Consequently `mu(P_(E,Z))>=mu(S_j(Z))` for every one retained cell.  A
projected completion with `k` aligned labels requires

```text
k=2: mu(P_(E,Z))<25/91,
k=3: mu(P_(E,Z))<36/91.                                  (2)
```

The matched-residue source of THM-3360/3376 is much narrower. Its vertices are
six pairs
`(e_i,q_i)` with distinct body labels `e_i in E` and positive levels `q_i`,
mapped to the reflected drifts

```text
z_i=Lq_i-e_i.                                             (3)
```

At the one common upper-median body-safe cell its vertex observable is

```text
mu(A_(Lq-e)(j))=1/7+e/[7(Lq-e)],                          (4)
```

its edge observable is the exact located overlap
`w_ij=mu(A_i(j) intersect A_j(j))`, and its intrinsic low/high relation is
formed from the reduced ratio of the **levels** `(q_i,q_j)`, not from the raw
drifts `(z_i,z_j)`. On the 649-body domain, THM-3360 proves every high edge has

```text
w_ij > DMAX/5,
DMAX=186636088362/11773143757375.                         (5)
```

Hunter is lawful only for edges on one common cell:

```text
mu(S_j(Z)) >= 1-sum_i mu(A_i(j))+sum_(ij in T) w_ij.      (6)
```

For the reflected `k=1` source there are six drifts, `T` has five edges, and
the target is `1/7`; this is exactly why five copies of `(5)` pay the whole
six-level debt.

For the projected sectors the edge budget is different:

```text
k=2: five drifts, four-edge tree, baseline safe mass 2/7, target 25/91;
k=3: four drifts, three-edge tree, baseline safe mass 3/7, target 36/91. (7)
```

Thus a literal five-edge transplant is not even combinatorially typed.

## 2. The first failure occurs before projection

The exact THM-3351 queue reconstruction left `13` complete families and
`113` rows at `z1=216`. THM-3361 has since removed three of them; this
historical `113`-row universe is an audited superset of the current `110`.
On that superset:

```text
ruler range                                      129360..5045040,
intersection with the THM-3360/3376 649 bodies                  0,
reflected robust-edge count                            15 on 113/113,
representations 216=qL-e, q>=1, e in E                          0. (8)
```

The body mismatch is structural: the 649-body bank has at most ten robust
reflected edges, whereas every remaining projected row has robust graph
`K6`.  The affine mismatch is stronger: `qL-e>=L-14>216` on every row.
Therefore the matched-residue source does not embed into even one of these
projected rows, including the current `110`. This is the first failed map;
no quotient theorem has yet been invoked.

The cheapest attempted repair is also false.  Let
`H=floor(13L/132)+1` and call the raw pair `(216,H)` high by its reduced
integer ratio.  Every one of the 113 pairs is high in that surrogate sense,
but at the upper-median body-safe cell its exact overlap histogram is

```text
overlap 0: 108 rows,
overlap 1:   5 rows.                                    (9)
```

The first frozen witness is atlas row `94`:

```text
E=(1,3,8,10,11,14), L=129360, j=65835, H=12741,
(216,H)/gcd=(72,4247), 72+4247=4319,
A_216(j)=[0,1], A_12741(j)=empty,
mu(A_216(j) intersect A_12741(j))=0.                    (10)
```

So a positive high-pair floor does not survive replacement of the reflected
affine rays by arbitrary projected drifts.

## 3. The first genuine quotient loss

Suppose one nevertheless generalizes the source to arbitrary drifts.  The
projected denominator passport records

```text
d=L/gcd(L,z),   z=(L/d)u+hL,   gcd(u,d)=1.               (11)
```

A located pair weight is not a function of `d`, or even of `(d,u)`.  In the
same row and cell as `(10)`, take

```text
z=12741, z'=12741+L=142101,
d(z)=d(z')=43120, u(z)=u(z')=4247 mod 43120.             (12)
```

Exact local intervals give

```text
mu(A_216 intersect A_12741)=0,
A_142101=(9625/15789,35035/47367),
mu(A_216 intersect A_142101)=6160/47367.                 (13)
```

Thus ray height is already load-bearing after denominator and unit have been
fixed.  The passport quotient destroys the pair observable and the low/high
graph.  A lawful projected Hunter state must retain at least `(d,u,h,j)` (or
the literal drift and cell).  Equation `(1)` then takes a union over `j` and
forgets cell address; independently maximized pair weights from different
cells cannot be added.  A common-cell or every-local-coordinate sidecar is
mandatory.  THM-3351's translated-band and denominator-two terminals retain
exactly this missing affine information.

## 4. Strongest survivor: four reflected residues close projected `k=3`

Restore the affine type `(3)`, but use only four distinct body labels and
allow their positive levels to repeat.  On any body-safe cell, subadditivity
and `(4)` give

```text
mu(S_j)>=3/7-D4,
D4=sum_(i=1)^4 e_i/[7(q_i L-e_i)].                       (14)
```

An exact all-`3003` body census (with the monotonicity `q_i>=1`) gives

```text
D4 <= 123896/5540535
```

at `E=(1,2,3,4,6,12)`, `L=168`, using labels `12,6,4,3` at level one.
Therefore

```text
mu(S_j)-36/91 >= 3/91-D4 >= 58759/5540535 > 0.          (15)
```

This proves, independently of any pair-floor theorem, that every projected
`k=3` packet
whose four drifts occupy four distinct reflected residue rays `Lq-e` closes
on every body. It removes none of the historical `113` or current `110` rows
because `(8)` shows that their first drift is outside this affine image.

## 5. Strongest survivor: five reflected residues close projected `k=2`

Now work on the 649 upper-median bodies and take five distinct body labels,
with drifts `Lq_i-e_i`.

If two levels repeat, the existing THM-2941 same-level graph is `K6` on all
649 bodies (the two chromatic exceptions are disjoint from this bank).  Its
good-pair floor exceeds the full six-label level-one debt.  One repeated
edge therefore gives `mu(S_j)>2/7>25/91`.

Assume the five levels are distinct. Sorted distinct levels satisfy
`q_i>=i`. Moreover, for `q<r`,

```text
f(e,q)-f(e,r)=eL(r-q)/[7(qL-e)(rL-e)]
```

increases with `e`, so the rearrangement maximum pairs larger labels with
smaller levels. Termwise monotonicity makes common scale two dominate every
larger scale. This proves the universal five-label debt bound

```text
D5 <= 183680141/11691304625                             (16)
```

at `E=(1,2,3,4,6,12)`, with labels `(12,6,4,3,2)` assigned to levels
`(1,2,3,4,5)`.  If the intrinsic high graph has two edges, those two edges
form a forest and extend to a `K5` tree.  Using `(5)--(7)`,

```text
mu(S_j)-25/91
 > 1/91 + 2(DMAX/5)-D5
 = 95318697414/58865718786875 > 0.                      (17)
```

In particular every disconnected-low five-level packet closes, but `(17)`
also handles connected-low packets having at least two high edges.

With at most one high edge the low graph is connected.  Duplicate-free
reverse search gives the complete five-level primitive atlas

```text
high edges  0    1    2    3     4     5      6
shapes      2   14   82  363  1400  4969  12559,
total 19389.                                             (18)
```

Thus only `16` primitive rays remain.  The debt plus zero/one high-floor gate
closes all of them from common scale `2` onward (weakest scale-two margin
`11573483686/3423705030837`), and closes `11` already at scale one.  The five
scale-one heads are

```text
(1,2,3,4,6), (1,2,3,4,12), (1,2,3,6,12),
(1,2,4,6,8), (1,2,4,6,12).                              (19)
```

The exact physical terminal scans every

```text
649 bodies * 5 shapes * 720 injections = 2,336,400      (20)
```

packet, evaluates all ten pair masses with the canonical THM-3352 engine,
and takes the maximum common-cell `K5` Hunter tree.  Failures are zero.  The
weakest exact margin over `25/91` is

```text
1107040893783422917/37818341551113357325>0              (21)
```

at

```text
E=(2,6,7,8,9,12), L=7056, j=3780,
shape=(1,2,3,4,6), slots=(1,4,0,3,2).                   (22)
```

At `(22)`, a literal rational interval merge reproduces all ten fast pair
masses and gives

```text
actual safe mass = Hunter lower bound
 =126337309517223212/415586170891355575.                 (23)
```

Consequently, composing with the proved floor `(5)`, **every projected `k=2`
packet on the 649-body bank whose five drifts occupy five distinct reflected
body-residue rays closes**, with arbitrary positive levels.  This is a
residue-subsector theorem, not a closure of arbitrary projected `k=2`.

## 6. Lawful next test for the current 110-row wall

The reflected high graph is the wrong carrier for the current wall.  The
next exact state should retain one common cell and literal ray tokens
`(d,u,h)`.  At fixed local coordinate, high drifts act as translated cyclic
bands on the body-safe cell set; their lawful pair sidecar is the
denominator/gcd CRT needle overlap of THM-2941 `(25p)--(25q)`, not the raw
drift-ratio graph. THM-3361 already closed the former `gcd72/L720720` rows
`191,228,332`. The next bank must be reconstructed from the live post-THM-3361
queue, retaining units and heights before comparing a common-cell or
every-coordinate Hunter/translated-band certificate against `36/91`.

## 7. Reproduction and frozen artifacts

```bash
python3 04-computation/lrc14_reflected_pair_projected_transfer_audit_20260814.py
python3 -O 04-computation/lrc14_reflected_pair_projected_transfer_audit_20260814.py
python3 04-computation/lrc14_reflected_k2_five_residue_heads_scan_20260814.py --processes 4
python3 -O 04-computation/lrc14_reflected_k2_five_residue_heads_scan_20260814.py --processes 4
```

Current SHA-256 values (LF-normalized bytes):

```text
3e2ab5103b211608b4d962fda3edccd75406c8907d139654ff762699051fcc91  transfer audit
17417e80deb5bdb9819a8071068f70f8595f2ebce100690ef6dd1369b53cda62  transfer output
f92234d91026e93111d81e4514b876c87cb09ce34f88c06b24a6a7fb74cf853c  five-head scan
6bb99f478e247814bf93aa3557e2df5e91176d067452494aef8bc876324eb2af  five-head output
```

Both artifacts have byte-identical ordinary and `python3 -O` outputs.  The
five-head scan also freezes semantic digest
`9a3a9ab7b779c990d41f908469159f52c02a5e09053330b45269fdf03b5f7685`.
