# LRC14 THM-2388 blocker-cage capacity escape

**Status:** FINITE-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.

THM-2388,
`THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law`,
is proved, exact, and independently hostile-audited. THM-1166,
`THM-1166-seven-wall-fano-gcd-discrepancy`, and THM-2263,
`THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor`,
are proved and independently audited. The same-parent tensor refinement
also uses proved and audited THM-2391,
`THM-2391-blocker-caged-septimal-single-layer-address-reduction`.

The result repairs an initially false inference. THM-2388's exact
`36/343` is an integrated multiplicity excess, not the measure of the
off-blocker toothpick set. A blocker-cage capacity estimate nevertheless
forces the following exact alternative on all `165` profiles:

```text
one of ten oriented cross-ancestor ratios at binary depth M=1,

or

clean-hole mass at least 1/26754,
with a canonical top-labelled charged cell of mass at least 1/1391208
and an owner-resolved C_7 x C_13 cell of mass at least 1/9042852.
```

Thirteen-adic pair caps alone first settle `132/165`; live septimal
typing sharpens this to `150/165`; coupling the last ancestor ratio and
computing its exact four-event union leaves the ten finite types.

## 1. The cage and the exact missing mass

Retain THM-2388's notation

```text
c_j=13 C_j,

B=D_(C_1) union D_(C_2) union D_(C_3),

K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i)),

Z={K=0},

X=D_(q_*)^c intersection D_(c_3)^c.                 (1)
```

THM-2388 proves

```text
mu(Z intersection X)>=36/343,                       (2)

Z intersection X
 subset T^(-1)(B minus D_(C_3)),       T(x)=13x.    (3)
```

Let

```text
mathcal C
 =B intersection T^(-1)(B minus D_(C_3))            (4)
```

be the blocker cage. Then

```text
(Z intersection X) minus B
 subset S_0:={y notin B:K(y)=0},                    (5)
```

and (2)--(4) give

```text
mu(S_0)>=36/343-mu(mathcal C).                       (6)
```

This is the correct dichotomy: either the cage has enough capacity to
absorb all holes, or a positive off-cage toothpick set exists.

## 2. Six pair caps bound the cage

Since

```text
T^(-1)(B minus D_(C_3))
 subset D_(c_1) union D_(c_2),
```

the union bound gives

```text
mu(mathcal C)
 <=sum_(i=1)^3 sum_(j=1)^2
      mu(D_(C_i) intersection D_(c_j)).              (7)
```

For the two same-owner pairs,

```text
c_j=13C_j
```

and Haar invariance under `x -> C_j x` reduce the intersection to
`D_1 intersection D_13`. The only positive-length overlap is the
central tooth, so

```text
mu(D_(C_j) intersection D_(c_j))=1/91.              (8)
```

Every other pair with distinct `13`-adic valuations is controlled by
THM-2263. At valuation gap `d>=1`, its sharp cap is

```text
d even:
  1/49+6/(49*13^d);

d odd:
  1/49+5/(588*13^d).                                (9)
```

## 3. Every repeated-first profile escapes

Using only thirteen-adic information, for a repeated-first profile

```text
(lambda_1,lambda_2,lambda_3)=(1,1,c),

5<=c<=19,                                           (10)
```

the two shallow cross pairs in (7) have gap one, and the two pairs from
`C_3` have gap

```text
d=c-2>=3.
```

Among `3<=d<=18`, the largest cap in (9) occurs at the even gap `d=4`:

```text
rho_gap(d)<=583/28561.                              (11)
```

Indeed the even and odd defects decrease separately, and the ratio of
the `d=4` defect to the `d=3` defect is `72/65`.

Equations (7)--(11) yield

```text
mu(mathcal C)
 <=2/91+2(23/1092)+2(583/28561)
 =17981/171366

 <36/343.                                           (12)
```

The exact margin is

```text
36/343-17981/171366
 =1693/58778538.                                    (13)
```

Thus, already before using the septimal lane,

```text
mu(S_0)>=1693/58778538>0                            (14)
```

for all `15` repeated-first profiles.

## 4. The thirteen-adic strict-profile sweep

For strict profiles

```text
(lambda_1,lambda_2,lambda_3)=(1,b,c),

2<=b<c,                    5<=c<=19,                (15)
```

the same exact ledger settles `117` of the `150` profiles.

Thirty profiles contain a cross pair with equal valuation:

```text
b=2

or

c=b+1.                                              (16)
```

The distinct-valuation cap (9) does not apply to that pair. Among the
remaining `120` profiles, only

```text
(1,3,6), (1,4,6), (1,4,7)                          (17)
```

miss the coarse six-pair threshold. The smallest positive strict margin
is

```text
67/2260713
```

at `(1,3,5)`. Therefore (14), with its slightly smaller margin, is
uniform over

```text
15+117=132
```

valuation profiles.

At this stage the other `33` strict profiles remain. The `30` rows in
(16) need quotient-lock/same-layer information, while (17) appears to
need an improvement over the six-pair union bound or use of the omitted
`D_(c_3)^c` factor.

There is a sharp warning against trying to settle (17) by
**thirteen-adic** pair-equality compatibility alone. The following
thirteen-adically compatible unit-part ladders have actual
four-cross-pair sums whose six-pair invoices still exceed `36/343`:

```text
profile       (u_1,u_2,u_3)   six-pair invoice       excess

(1,3,6)       (1,12,12)       251891/2399124         4307/117557076
(1,4,6)       (1,1,12)        20991/199927           363/9796423
(1,4,7)       (1,1,1)         273066/2599051         13686/127353499.
```

Here `c_i=13^(lambda_i)u_i`. Thus the fact that simultaneous even-gap
maxima force common unit parts is insufficient: an odd-gap upper endpoint
can coexist with one sacrificed even endpoint and keep the invoice above
threshold. These ladders deliberately omit THM-2388's live septimal
separation. They are hostiles only to a **thirteen-adic pair-sum**
shortcut, not to the result below.

## 5. Septimal cancellation settles all but one exact lane

THM-1166's folded formula has a useful cross-prime specialization:

```text
nu_7(a)!=nu_7(b)  implies
mu(D_a intersection D_b)=1/49.                      (18)
```

Indeed, after dividing by `gcd(a,b)`, one reduced coefficient is `0` or
`7` modulo `14`, while the other is a seven-unit. Therefore the two
folded arguments in THM-1166 coincide modulo `14`, and its correction
term vanishes exactly.

In the live THM-2388 lane,

```text
nu_7(C_3)=nu_7(c_3)>M>
nu_7(c_1),nu_7(c_2).                                (19)
```

Both `C_3/c_j` cage pairs in (7) are consequently exactly `1/49`.
For a strict profile with `b>=3`, the live invoice is

```text
2/91 + rho_gap(b) + rho_gap(b-2) + 2/49.            (20)
```

It is below `36/343` for all `135` such profiles. The worst family is
every profile with `b=4`; its exact off-cage margin is

```text
6042/9796423.                                       (21)
```

For the `15` repeated-first profiles, the exact live invoice and margin
are

```text
2/91+2(23/1092)+2/49 =401/3822,

36/343-401/3822 =1/26754.                           (22)
```

Thus positive off-cage toothpick mass is unconditional on

```text
15 repeated-first + 135 strict =150 profiles,        (23)
```

with uniform floor `1/26754`.

The remaining `15` profiles have `b=2`, so `C_2` and `c_1` have the same
thirteen-adic valuation. If their septimal valuations differ, (18) gives
the stronger invoice and margin

```text
2/91+25/1183+3/49 =864/8281,

36/343-864/8281 =36/57967.                          (24)
```

The exact residual lane is therefore

```text
b=2 and nu_7(C_2)=nu_7(c_1).                        (25)
```

This is an ancestor-equality/same-septimal-layer problem, not another
generic gap-cap problem. It is already finite in reduced ratio. Write

```text
(C_2,c_1)=g(a,b),

gcd(a,b)=1,       a<=b,       7*13 does not divide ab.
```

All cage terms except this pair total `695/8281`. Thus failure of the
off-cage inequality forces

```text
rho(a,b)>=1219/57967,

Delta(a,b)/(ab)>=144/1183,
```

where `Delta` is THM-1166's folded numerator. Since `Delta<=49`,

```text
ab<=floor(49*1183/144)=402.
```

This first, deliberately pairwise enumeration leaves `124` unordered
reduced ratios. In fact

```text
ab<=345,       a<=11,       b<=197
```

on this bank, and the census by `a` is

```text
a:count =
1:49, 2:20, 3:22, 4:11, 5:11,
6:1, 8:2, 9:3, 10:2, 11:3.
```

Two compatibility steps make the bank much smaller.

First, the gap-two cross pair is determined by the same ratio. In oriented
coordinates write

```text
(C_2,c_1)=13h(a,b).
```

Then

```text
(C_1,C_2,c_1,c_2)=h(b,13a,13b,169a).
```

Because `169=1 mod 14`, THM-1166 gives the same folded defect in
`rho(a,b)` and `rho(b,169a)`, with the latter correction divided by
`169`. The absorption threshold becomes

```text
Delta(a,b)/(ab)>=156/595,
```

so `ab<=186` analytically. Exact enumeration leaves `52` unordered
ratios, with the sharper bounds

```text
ab<=177,       a<=10,       b<=85.
```

Second, the four low-low pair events should not be union-bounded
separately. Their union is exactly

```text
U_(a,b)
 =(D_b union D_(13a))
   intersection
   (D_(13b) union D_(169a)).                        (26)
```

The `C_3` part of the cage costs at most `2/49` by septimal
separation. Hence absorption requires

```text
mu(U_(a,b))>=36/343-2/49=22/343.                    (27)
```

A direct exact interval-union calculation on the common endpoint grid
of order `14*169*a*b` leaves only ten oriented ratios:

```text
(C_2:c_1)     mu(U)

1:1           193/1183
1:2           114/1183
2:1           239/2366
1:3           263/3549
3:1           95/1183
4:1           331/4732
2:3           43/546
3:2           95/1183
3:4           491/7098
4:3           331/4732.
```

Thus the final capacity dichotomy is a charged off-cage cell or one of
only six unordered ratios

```text
{1:1,1:2,1:3,1:4,2:3,3:4},
```

with orientation retained as above. These are a tiny
Farey/consecutive-ratio cage. They still need the septimal heavy-word and
endpoint/owner transport sidecars, but the former `33`-profile continuum
has collapsed to ten explicit oriented address types.

The dichotomy is uniform away from those ten types. Among all admissible
oriented pairs with `ab<=186` below the threshold (27), the smallest
direct-union margin is

```text
9/22295
```

at the `4:5/5:4` boundary. Beyond that finite range, the largest possible
folded quotient is

```text
Delta/(ab)<=16/73,
```

with equality at `(a,b)=(3,73)`; the `49/224` analytic tail is already
smaller. This gives off-cage margin at least

```text
934/4231591
```

for every nonbank `b=2` ratio. Both margins exceed the repeated-profile
floor `1/26754`. Therefore, across all `165` profiles:

```text
one of ten oriented cross-ancestor types

or

off-cage toothpick mass at least 1/26754.            (28)
```

The next subsection uses THM-2391 to eliminate this bank for `M>=2`,
so the first alternative ultimately survives only at `M=1`.

The least-used sidecar now becomes completely explicit. Apply THM-2377's
balanced Bockstein descent to the two same-layer originals:

```text
c_1=13hb,      c_2=169ha,

(c_2-rho c_1)/7=13h d,

rho in {-3,-2,-1,1,2,3}.
```

On the ten oriented types, `(rho,d)` is respectively

```text
1:1 -> (-1,2)     1:2 -> ( 3,1)
2:1 -> (-2,4)     1:3 -> ( 2,1)
3:1 -> (-3,6)     4:1 -> ( 3,7)
2:3 -> (-3,5)     3:2 -> ( 2,5)
3:4 -> ( 1,5)     4:3 -> ( 1,7).
```

Thus the abstract carry has collapsed to `d in {1,2,4,5,6,7}`: two
types descend exactly to `C_2`, and two climb one septimal digit. This
does not itself empty a type. It is the cheapest exact interface to
THM-2390's labelled `W=7/8` heavy words.

### 5.1 The top comb closes the sole `M=2` cage

THM-2391 adds a depth filter to the ten types. For `M>=2`, its two
low-blocker progressions have a common absolute step. Since

```text
C_2/C_1=13a/b,
```

the oriented address must obey

```text
13a=+/-b mod 7^M.                                   (29)
```

At modulus `49`, only `(a,b)=(4,3)` survives the ten-row table. At
modulus `343`, none survives. Thus `M>=3` is already empty, while `M=2`
has the exact low union

```text
U
 =(D_3 union D_52)
   intersection (D_39 union D_676),

mu(U)=331/4732.                                     (30)
```

The top-comb exclusion in `X` removes more mass than the union-bound
slack in (30). Write the common low scale and the top speed, after
dividing their gcd, as coprime integers `p,n`. The live typing gives

```text
7 does not divide p,       49 divides n,

13 does not divide pn.                              (31)
```

Put `F=1_U`. Its exact endpoint grid has `139` linear pieces, with the
first and last joined through the circle endpoint. Hence `U` has `138`
circular components and

```text
Var(F)=276.
```

For the centered danger interval `d=1_D`, normalized Fourier
coefficients obey

```text
|Fhat(k)|<=138/(pi |k|),

|dhat(l)|<=1/(pi |l|).                              (32)
```

Since `gcd(p,n)=1`, the common Fourier frequencies are
`(k,l)=(nt,-pt)`. Therefore

```text
|mu(U_p intersection D_n)-mu(U)/7|
 <=sum_(t!=0) 138/(pi n|t|) * 1/(pi p|t|)
 =46/(np).                                         (33)
```

The low union plus the two high septimally separated cage pieces
exceeds the hole floor by only

```text
331/4732+2/49-36/343
 =1347/231868.                                      (34)
```

Meanwhile

```text
mu(U)/7-1347/231868
 =485/115934.                                       (35)
```

Thus (33) closes every `np>=10996`; at the first tail integer its exact
reserve is `12/159351283`.

It remains only the finite core. Write `n=49m`. Then `np<10996` implies
`pm<=224`. An exact interval audit of all `753` coprime primitive pairs
subject to (31), deliberately including higher-seven-depth controls,
has the unique minimum

```text
mu(U_p intersection D_n)
 >=1849/231868,

equality only at (p,n)=(1,49).                       (36)
```

The actual cage is contained in

```text
(U_p intersection D_n^c)
 union the two high cross pieces.
```

Equations (34)--(36) consequently give the uniform clean-hole floor

```text
delta>=1849/231868-1347/231868
      =251/115934>0.                                (37)
```

This closes the sole `M=2` address. Its top-labelled `/52` cell has
mass at least `251/6028568`, and its owner-resolved `/338` tensor cell
has mass at least `251/39185692`. The only cage boundary left by this
entire argument is therefore the ten oriented table at binary depth
`M=1`.

## 6. Canonical top-labelled charge and a same-parent bi-root tensor

On `S_0 intersection X`, every quotient blocker is safe, so all original
blockers are safe on the thirteen inverse roots. Scalar-cover exactness
and THM-2388's transfer identity give:

```text
guard root count =4,

each ordinary q root count=2,

one inverse root has unit multiplicity two,

the other twelve have multiplicity one.             (29)
```

A general deterministic ordinary-label deletion may choose the first
`q` outside the overlap pair, giving the adjacent two-root word used in
the first version of this note. There is a sharper canonical choice:
delete the distinguished septimal top label `q_*` itself. Let `W` be the
union of the other five unit masks on the inverse-root fibre and put

```text
A=1-W.                                              (30)
```

If `q_*` is outside the unique double pair, both of its roots are
exclusive and `A` is its adjacent two-root word. If `q_*` belongs to the
double pair, one root is covered by the partner and the other is
exclusive, so `A` is a singleton. Thus in both cases

```text
W=1-A,

A is singleton or an adjacent q_*-scaled two-set.   (31)
```

With normalized `C_13` Fourier transform, every nonzero mode of either
word is nonzero and `What(a)=-Ahat(a)`. The exact ledgers are

```text
singleton:
  sum_(a!=0)|Ahat(a)|^2=12/169,
  sum_(a!=0)|Ahat(a)|^4=12/28561;

adjacent:
|Ahat(a)|^2
 =(2+2cos(2 pi a delta/13))/169
 >0                         for every a!=0,

  sum_(a!=0)|Ahat(a)|^2=22/169,

  sum_(a!=0)|Ahat(a)|^4=62/28561.                   (32)
```

Partition a positive off-cage set by:

- a deterministic lower blocker owner, at most `2` choices;
- singleton/adjacent status, `2` choices; and
- the support of `A`, `13` choices in either status.

Thus one canonical top-labelled cell has mass at least

```text
(1/26754)/(2*2*13)
 =1/1391208.                                        (33)
```

The older arbitrary-label construction remains a positive control with
`2*15*13=390` cells and floor `1/10434060`, but is never needed for the
sharper conclusion. On the top-labelled cell the phase and support in
(32) are constant, so integration cannot cancel any of the twelve
nonzero target colours.

THM-2391,
`THM-2391-blocker-caged-septimal-single-layer-address-reduction`, adds a
second root direction on the **same parent**. Every clean parent is
`q_*`-safe and `c_3`-safe, and all seven translates `y+s/7` retain those
two high-safe conditions. The primitive lower load

```text
L_7(s)
 =1_(E_H)(y+s/7)
  +sum_(q_i!=q_*)1_(D_(q_i))(y+s/7)
  +1_(D_(c_1))(y+s/7)+1_(D_(c_2))(y+s/7)
```

covers all seven siblings and has total incidence eight. Consequently

```text
L_7(s)-1=1_(s=d(y))                                 (34)
```

for a unique `d(y) in F_7`. At `s=0`, the six unit masks vanish because
`K(y)=0`, while the scalar cover forces at least one low blocker.
Therefore

```text
d=0     iff both low blockers are active;

d!=0    iff exactly one low blocker is active.       (35)
```

There are exactly thirteen owner/address categories: one simultaneous
category at `d=0`, and `6*2` nonzero-address/exclusive-owner categories.
Together with the `26` singleton/adjacent supports in (31), one
owner-resolved same-parent `C_7 x C_13` cell has mass at least

```text
(1/26754)/(13*26)=1/9042852.                         (36)
```

On a fixed such cell, the collapsed bi-root tensor is literally

```text
1_(s=d) A(r).
```

Its normalized Fourier coefficient is

```text
(rho/7) zeta_7^(-ell d) Ahat(k),
```

so it is nonzero for every `ell in F_7` and every `k!=0`. This is a
coefficient-derived owner/address/target tensor on one common parent,
not merely two unrelated positive marginals.

## 7. Scope and remaining transport

The source-to-target map in this note is

```text
THM-2388 hole mass and cage inclusion
 + THM-2263 gap-sensitive pair caps
 + THM-1166 septimal folded cancellation
 + exact low-union geometry
 -> ten oriented cross-ancestor types at M=1
    or clean-hole mass at least 1/26754
 -> canonical q_* deletion charge
 -> THM-2391 owner-resolved C_7 x C_13 tensor.        (37)
```

On the charged side it preserves the distinguished top label, exact
low-owner word, septimal double address, root support/status, and all
twelve nonzero target colours. It discards the later terminal-word
endpoint and does not identify the tensor with one canonical THM-2305
current or fixed `(m,X)` triangle. The ten cage types retain only their
oriented ratio and Bockstein descendant. A target-shift-covariant base
filtration is still required on both sides.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open. This note must not be cited as a proved dependency until
the present implication is independently audited and promoted.

## Exact reproduction

Run

```bash
python3 04-computation/lrc14_thm2388_cage_capacity_escape_candidate.py
python3 -O 04-computation/lrc14_thm2388_cage_capacity_escape_candidate.py
```

and compare both transcripts after LF normalization with

```text
05-knowledge/results/lrc14_thm2388_cage_capacity_escape_candidate.out
```

The companion imports no interval geometry beyond the exact formulas
stated above. It checks every repeated and strict valuation profile,
the thirteen-adic `132/165` first pass, the live septimal `150/165`
split, the exact residual lane, its `124` pairwise bank, its `52`
coupled bank, and the final ten-orientation exact-union bank, the three
thirteen-adic pair-sum hostiles outside the live septimal typing, the
THM-2391 depth filter, `138`-component whole-union BV tail, and all
`753` primitive top-comb pairs closing the sole `M=2` address, the
fallback `/390`, canonical top-labelled `/52`, and owner-resolved
bi-root `/338` cell invoices, and all singleton/two-root energy
constants.

LF-normalized SHA256:

```text
script 5eb069ff9a2f1df2c572251999b58150af7f08aa62915ceb653acf673c3838f3
output 3b3cb4c418d81c596ac2cfce798200a69d818b6e37842551ebb945abc3dede7e
```

Independent audit is pending.
