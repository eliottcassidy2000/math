---
id: THM-2135
title: "Root-profile invoices and unit-annulus closure for the scalar flood tails"
status: >
  PROVED. The two scalar tails isolated by THM-2133 have sharper thirteen-root
  profiles. In the six-unit/one-deep case, four unit singleton phases are
  impossible everywhere and three are impossible on the guard-danger phase;
  this gives a triple cyclic-contraction invoice d|420v. Exact unit-annulus
  certificates force 13^4|v. In the five-unit/two-deep case, all ten unit-unit
  edges and the five guard-unit edges carry a two-colour
  gcd cover. Goodman counting forces at least two monochromatic divisibility triangles. A root-needle
  multiplicity reformulation moreover forces a depth-at-least-six flood on
  measure at least 19/546 inside the two deep combs. Finite exact capacity
  arguments close deep valuation pairs (1,1) and (1,2). This theorem originally
  left higher deep valuations open; THM-2138 now closes both scalar tails at
  every depth. The fivefold guard pencil and higher ranks remain, so LRC(14)
  is still open.
source: codex-2026-07-22-LRC-scalar-flood-root-profile
depends_on:
  - THM-2080
  - THM-1166
  - THM-2133
related:
  - THM-2128
  - THM-2138
script:
  - 04-computation/lrc14_scalar_flood_mod169_codex_20260722.py
  - 04-computation/lrc14_five_plus_two_unit_annulus_depth2_codex_20260722.py
  - 04-computation/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.py
output:
  - 05-knowledge/results/lrc14_scalar_flood_mod169_codex_20260722.out
  - 05-knowledge/results/lrc14_five_plus_two_unit_annulus_depth2_codex_20260722.out
  - 05-knowledge/results/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.out
script_sha256:
  - bda80cbbe785fcda5a1c57409b70a1564697a11cccb3f2ac5106dc33ac4d846f
  - bfdfde772c1dcc7f76805f9c9416aba35744bab3287f12c8e24244d37a9d10fa
  - 9d8774e3179ce34d8bb39f9bb3c607b7fee897e5f3429106b5043f637a4289bc
output_sha256:
  - 1720b141077ba5c546ca44b6c2df0991c52c47b69a34be4e6be815d172676cde
  - 42d08d7b604def211f7f930889bec3ddcbbb4fe8f6b82d595a7834a46d8637a1
  - c62e35c90e51b82620d7cd977d1cfa1ec7f85226356ad2e6cd28fd9a44542dae
hash_basis: working-tree bytes with LF line endings
---

# THM-2135 -- root profiles of the two scalar flood tails

Put

```text
C_H={t:||Ht||>1/7},          E_H={t:||Ht||<1/7},
D_r={t:||rt||<1/14}.                                  (1)
```

Retain either scalar containment isolated in THM-2133:

```text
(I)  C_H subset D_(13v) union union_(i=1)^6 D_(q_i),

(II) C_H subset D_(13v_1) union D_(13v_2)
             union union_(i=1)^5 D_(q_i),             (2)
```

up to a null set, where all displayed coefficients are positive, `H` is odd,
`H` and every `q_i` are thirteen-units, the `q_i` are distinct, and in case
(II) `v_1!=v_2`. Endpoints are discarded whenever a root count is made. Any strict violation then persists on an interval, so the almost-
everywhere convention loses none of the containments below.

## 1. The exact thirteen-root profile

Fix a phase `y` and let `R_y={x:13x=y}`. For a thirteen-unit `a`, direct
inspection of thirteen equally spaced points gives

```text
#{x in R_y:||Hx||<=1/7}
 =3 if y in E_H, and 4 if y in C_H,                  (3)

#{x in R_y:||ax||<1/14}
 =1 if y in D_a, and 2 if y notin closure(D_a).      (4)
```

The endpoints omitted above are finite. A deep tooth is constant on the
needle:

```text
13v x=v y                         for every x in R_y. (5)
```

In case (I), suppose `y` is safe for `v` and let

```text
S(y)=#{i:y in D_(q_i)}.
```

The six unit masks have total root capacity `12-S(y)`. They must cover ten
guard-safe roots over `E_H` and nine over `C_H`. Therefore

```text
y in E_H  implies S(y)<=2,
y in C_H  implies S(y)<=3.                           (6)
```

Equivalently, with closures only at the target boundary,

```text
intersection_(i in I) D_(q_i) subset closure(D_v)     |I|=4,
E_H intersection intersection_(i in J)D_(q_i)
                         subset closure(D_v)          |J|=3. (7)
```

In case (II), if `y` is safe for both `v_1,v_2`, the five unit masks have
capacity `10-S(y)`. Equations (3)--(4) now force

```text
E_H intersection D_(q_i)
       subset closure(D_(v_1)) union closure(D_(v_2)),          (8)

D_(q_i) intersection D_(q_j)
       subset closure(D_(v_1)) union closure(D_(v_2))   (i!=j). (9)
```

For the single target in (7), equality at a finite common-kernel point can be
perturbed outward while all source phases remain strict. Thus that target may
be treated as strict in Section 2. The two closed targets in (8)--(9) must
remain closed: opposite boundary contacts can cover opposite sides of a source
neighborhood. THM-2133's cross lemma was proved with exactly those closed arcs.

## 2. The six-plus-one triple contraction

The fourfold part of (7) recovers THM-2133's invoice

```text
gcd(q_i:i in I) divides v                  for every |I|=4. (10)
```

There is a stronger guard-coloured triple invoice. Fix a triple `J` and put

```text
d=gcd(q_i:i in J),                 g=gcd(H,d).         (11)
```

Restrict the second containment in (7) to the order-`d` subgroup of the
circle on which the three `q_i` vanish. It says

```text
{k mod d:||Hk/d||<1/7}
 subset {k mod d:||vk/d||<1/14}.                     (12)
```

We use the following finite cyclic contraction.

> **Cyclic contraction lemma.** If (12) holds, then `g|v`, and either
> `d|v` or `d/g<=7`.

Indeed, the kernel of multiplication by `H` in `Z/dZ` has order `g`.
If `v` did not kill it, its nontrivial cyclic image would contain a phase of
norm at least `1/3`, contradicting (12). Hence `g|v`. Write

```text
d=gm,             H=gH_0,             v=gv_0,
gcd(H_0,m)=1.                                        (13)
```

After multiplying indices by `H_0^(-1) mod m`, (12) becomes

```text
||j/m||<1/7       implies       ||a j/m||<1/14,       (14)
```

where `a=v_0 H_0^(-1) mod m`. If `8<=m<=14`, the input `j=1` already
contradicts (14) unless `a=0`. If `m>=15`, choose the least absolute
representative `1<=A<=m/2` of a nonzero `a`. The input `j=1` first forces
`A<m/14`. Now put

```text
j=ceil(m/(14A)).                                     (15)
```

Then `j<m/7`, while

```text
m/14<=jA<m/14+A<m/7<m/2.                             (16)
```

This violates (14). Thus `a=0`, so `d|v`, whenever `m>=8`, proving the
lemma. Since every integer at most seven divides `420`, its weaker uniform
form is

```text
gcd(q_i,q_j,q_k) divides 420v             for every triple. (17)
```

In particular, every prime `p>=11` that divides three unit coefficients to
orders `a_i,a_j,a_k` divides `v` to at least their minimum order. This is a
genuine higher-arity valuation invoice; the four-uniform ledger alone does
not see it.

There is also a safe local-width consequence because (7) has only one target.
Order the unit coefficients as `q_(1)<...<q_(6)`. The central component of
the four smallest unit bands has halfwidth `1/(14q_(4))`, so (7) forces

```text
v<=q_(4).                                             (17a)
```

The central component of `E_H` and the three smallest unit bands has halfwidth
`min(1/(7H),1/(14q_(3)))`. If `H<2v` and `v>q_(3)`, both source widths exceed
`1/(14v)`, contradicting the single-target containment. Hence

```text
H<2v       implies       v<=q_(3).                    (17b)
```

The single-target hypothesis is essential; Section 3 records why two
high-frequency target combs can bridge one another's central tails.

## 3. The five-plus-two divisibility-coloured `K_6`

Apply the compact-subgroup cross lemma from THM-2133 to (8)--(9). It gives

```text
gcd(H,q_i) divides v_1 or v_2,                         (18)
gcd(q_i,q_j) divides v_1 or v_2.                      (19)
```

Thus every edge can be assigned a colour `s in {1,2}` satisfying

```text
edge {q_i,q_j}: gcd(q_i,q_j)|v_s,
edge {H,q_i}:   gcd(H,q_i)|v_s.                       (20)
```

These are all edges of the complete graph on

```text
{H,q_1,...,q_5}.                                     (21)
```

Choose one qualified colour on every edge. In fact at least two triangles are
monochromatic. If `r_x,b_x` are the red and blue degrees at vertex `x`, every
nonmonochromatic triangle is counted at two mixed vertices, so their number is

```text
(1/2) sum_x r_x b_x <=(1/2)*6*6=18.                  (22)
```

There are twenty triangles in `K_6`. Hence at least two pay the following
invoice. If a monochromatic triangle has labels `a,b,c`, then for one `s`

```text
lcm(gcd(a,b),gcd(a,c),gcd(b,c)) divides v_s.          (23)
```

For every prime, the valuation on the left is the median of the three vertex
valuations. This is a higher-order divisibility invoice, not a magnitude
bound: high-frequency target teeth can bridge the two tails of a small source
component, so no inequality such as `v_s<=max(a,b)` follows from (8)--(9).
The minimal hostile witness is

```text
D_1 intersection D_7=(-1/98,1/98) subset D_8 union D_105, (23a)
```

because the positive `D_105` tooth is `(13/1470,15/1470)` while
`13/1470<1/112` and `15/1470=1/98`. Both target coefficients exceed seven.

## 4. The faithful Kakeya needle in the five-plus-two tail

Remove the five unit danger combs from the guard-safe set and put

```text
L=C_H minus union_(i=1)^5 D_(q_i),
N(alpha)=#{x in L:13x=alpha}.                         (24)
```

The scalar cover (II) is equivalent, up to endpoints and null sets, to

```text
support(N) subset D_(v_1) union D_(v_2).              (25)
```

Indeed, both deep values are constant on the root needle by (5). Conversely,
outside `L` a unit tooth or the guard already settles the point. The function
`N`, rather than its support alone, is the multiplicity sidecar lost by the
quotient `[13]`.

Hunter's guard-centred star and THM-2080 give

```text
measure(L)>=5/42,              integral N=13 measure(L)>=65/42. (26)
```

THM-1166's pair floor gives

```text
measure(D_(v_1) union D_(v_2))
 <=2/7-1/91=25/91.                                  (27)
```

Thus the conditional mean of `N` on its support is at least

```text
(65/42)/(25/91)=169/30.                              (28)
```

Also `N<=10`, because every thirteen-root needle has at least three
guard-unsafe points. If

```text
h_6=measure{alpha:N(alpha)>=6},                       (29)
```

then (25)--(27) imply

```text
65/42<=5(25/91-h_6)+10h_6,
h_6>=19/546.                                         (30)
```

Every survivor therefore needs a positive-measure depth-at-least-six flood
inside two deep combs. This is stronger than the pair-gcd graph and identifies
the next object to attack: the high-multiplicity root needles, with their
affine masks retained.

## 5. Exact mod-169 closure of the first deep layer

Suppose first that case (I) has `13 not|v`. Evaluate the cover on
`169`-torsion and multiply the numerator by `H^(-1) mod 169`. The deep tooth
is safe whenever the normalized numerator is nonzero modulo thirteen. The
remaining guard-safe universe is

```text
U={z mod 169:13 not|z and ||z/169||>1/7},             |U|=110. (31)
```

A unit coefficient `r` gives the exact mask

```text
S_r={z in U:min(rz mod 169,-rz mod 169)<=12}.         (32)
```

Since every modulus used below is a power of thirteen, neither a guard
boundary `7||z||_N=N` nor a terminal boundary `14||az||_N=N` can occur.
Every uncovered torsion point is therefore strict and thickens to an open
failure of the almost-everywhere cover.

There are `78` unit classes modulo sign and `77` distinct restricted masks.
The frozen exact include/exclude search proves

```text
max_(r_1,...,r_5)|union S_(r_i)|=88,
max_(r_1,...,r_6)|union S_(r_i)|=96.                 (33)
```

Its branch bound adds the largest remaining individual marginal gains and
therefore overestimates every continuation. It enumerates every mask subset;
repeated residue masks cannot enlarge a union. Both maxima are below `110`,
so (33) contradicts case (I). Hence its remaining deep coefficient obeys

```text
13 divides v.                                        (34)
```

The unit-annulus obstruction continues for four valuation levels. Write

```text
v=13^s w,             13 not|w,             N=13^(s+2). (35)
```

On the torsion grid `t=z/N`, normalization by the unit `H` preserves
`13 not|z`, and

```text
(13v)t=wz/13.                                        (36)
```

Thus the deep tooth is safe throughout

```text
U_N={z mod N:13 not|z and 7||z||_N>N},               (37)
```

and a cover would force six multiplicative translates

```text
S_a={z in U_N:14||az||_N<N}                          (38)
```

to cover `U_N`. Repeated residue masks only weaken the union.

The companion counts every sign class by an exact Euclidean floor sum. If
`L=floor(N/7)` and `R=floor((N-1)/14)`, the full-grid count on
`L<z<N-L` is a sum of the two indicators

```text
az mod N <=R,                 az mod N >=N-R.         (39)
```

Each indicator is a difference of sums
`sum floor((az+b)/N)`. Removing the multiples `z=13z'` subtracts the
identical full-grid count at `N/13`. The Euclidean recurrence is checked
against its definition in 512 hostile cases, exhaustively against literal
masks at depths two and three, and on every extremizer, runner-up, and a
deterministic hostile sample at depths four and five. The exact census is

```text
m    N       |U_N|   sign classes   maximum   maximizers       runner-up
2    169       110        78           20     eleven classes       18
3    2197     1450      1014          240     6,183                236
4    28561   18830     13182         3140     6,2380              3084
5    371293 244810    171366        40800     6,30941            40058. (40)
```

At `m=2`, the stronger union search gives `96<110`. At `m=3` and `m=5`,
six times the maximum is respectively `1440<1450` and `244800<244810`.
At `m=4`, any packet containing a nonmaximal mask has total capacity at most

```text
5*3140+3084=18784<18830.                             (41)
```

An all-maximal packet uses only the two displayed masks, whose union has size
`5758`. Hence `s=0,1,2,3` are all impossible, and case (I) satisfies the much
stronger invoice

```text
13^4 divides v.                                      (42)
```

For case (II), the five-mask conclusion has a short finite exact capacity
certificate. A unit
danger mask has `25` points on `Z/169Z`. At least five lie in the guard-unsafe
interval: if `b` is the least absolute inverse of its coefficient, two
distinct integers `k_1,k_2 in {1,...,12}` have

```text
||b k_j||_169<=24;                                   (43)
```

the finite ranges are checked by the companion and summarized by the mask-
size histogram

```text
usable size:  0   8   12   16   18   20
number:       1   1    3   21   40   11.             (44)
```

Negation supplies the other two unsafe points, together with zero. Thus each
unit mask uses at most twenty points of `U`, and five use at most `100<110`.
If both `v_1,v_2` are thirteen-units, both deep teeth are safe throughout
`U`, which is impossible. Therefore

```text
13 divides v_1 v_2.                                  (45)
```

The unequal next layer also fails. Suppose after reordering that

```text
13 not|v_1,                 v_2=13w,       13 not|w.  (46)
```

On `Z/2197Z`, the tooth `13v_2=169w` is safe off the zero column modulo
thirteen. For each of the `78` sign classes of `v_1 mod 169`, also remove the
danger mask of `13v_1`. The surviving universes have the exact histogram

```text
size:       1210  1218  1222  1226  1228  1230
classes:       2     1     3    21    40    11.       (47)
```

Scanning all `1014` unit sign classes modulo `2197`, the minimum over the 78
rows of

```text
|U|-5 max_q |U intersection D_q|                    (48)
```

is `30`. The unique worst row is `v_1=14`, with `|U|=1230` and maximum unit
mask `240` at `q=183`. Thus five unit masks cannot cover, even before their
overlaps are charged. The deep valuation pair `(1,2)` is impossible.

The three companions recompute the masks, capacity tables, floor sums, and
exact maxima. Normal and optimized runs agree byte for byte.

## 6. Scope and Tournament Analysis

For the divided scalar coefficients, this theorem forces the `6+1` profile
to depth at least five and closes `5+2` profiles `(0^5,1,1)` and
`(0^5,1,2)`. At this stage the remaining valuation profiles were

```text
(0^6,a),       a>=5,
(0^5,1,b),     b>=3,
(0^5,a,b),     2<=a<=b.                              (49)
```

THM-2138 subsequently proves the all-depth unit-annulus extremal law, bounds
every positive-valuation mask by the unit maximum with the needed even-depth
gap, and eliminates every profile in (49). The triple contraction, two Ramsey
triangle invoices, and root-flood constraint (30) remain valid structural
information, but (49) is no longer an open branch.

For case (I), `U_N` and `S_a` live in the multiplicative unit group
`(Z/13^m Z)^*`, and `S_a=a^(-1)S_1`. Equivalently, a mask count is the number
of thirteen-primitive points of the congruence lattice

```text
L_a={(x,y):y=a x mod 13^m}                           (50)
```

inside the rectangle with radii `N/7` and `N/14`. At depths `m=3,4,5`, the
observed extremizers are the short lattice directions `(1,6)` and `(12,-1)`,
namely `a=6` and `12a=-1`; the unique runner-up obeys `11a=1`. THM-2138 proves
this extremal law at every depth and uses it to close the entire `6+1` tail.

Candidate vertices were the seven scalar teeth, the six labels
`{H,q_1,...,q_5}`, the 78 residue masks, the 110 torsion points, and the
thirteen-root needles. The pair observable in Section 3 is a two-valued gcd
absorption relation. An edge may carry both colours, so it is not a
tournament; choosing one colour is safe only for the Ramsey consequence and
forgets the alternate divisor. The mod-169 observable is symmetric mask
intersection size. Orienting ties by labels supplies a Hamiltonian path but
destroys the zero cuts and union predicate. The faithful objects are the
two-coloured complete graph with divisor sidecars, the undirected weighted
mask graph with its bitsets, and the Kakeya needle multiplicity `N`.
