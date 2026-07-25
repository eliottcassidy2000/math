---
id: THM-2273
title: "Shallow-owner flow and deep-successor gap spread"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT DIRECT-INTERVAL/ROOT REFEREE.
  In every strict first-depth-one scalar profile (1,b,c), removing the
  deepest blocker leaves two labelled shallow exclusive-owner flows of
  total mass at least 5696989/367580070. Their exact time-two split and
  the middle owner's expiration combine to give a common-time image at
  time b+1 of mass at least 5696989/76962600>1/14, wholly outside the
  named deepest-successor danger comb. If its successor speed is
  s=13^(c-b-1)u_3, that image occupies at least
  ceil(39878923*s/461775600) distinct positive-mass safe gaps. Hence all
  135 nonadjacent strict rows force at least two gaps; the fifteen b=2
  rows force at least fifteen. On the fifteen adjacent rows the image is
  uniformly larger than 1/13, and either u_3<=11 or it hits at least two
  gaps. This supplies a quantitative ancestry/gluing sidecar but no
  successor cover, excludes no profile, and does not prove LRC(14).
source: codex-2026-07-25-shallow-owner-gap-spread
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
related:
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
script:
  - 04-computation/lrc14_shallow_owner_deep_gap_spread_thm2273.py
  - 04-computation/lrc14_shallow_owner_deep_gap_spread_referee_thm2273.py
output:
  - 05-knowledge/results/lrc14_shallow_owner_deep_gap_spread_thm2273.out
  - 05-knowledge/results/lrc14_shallow_owner_deep_gap_spread_referee_thm2273.out
script_sha256:
  - 6c2b5e07eb07ea50421ae0a21adf655b26a89adfaef879241e17f082b91f5315
  - c0b0f7a6dfdd67da169ba1adc7cce554b12b662eff591466656b23e6d68f3908
output_sha256:
  - e10b58d1971c31efcf32c4f40874f56f1fef180c06635eb6ea2ce02bc498f9ba
  - 0cf2e27ca7c6acd564ecae06f54e7c5ee1a35d1a426103cd3725ab7b8d61f935
hash_basis: working-tree bytes (LF)
---

# THM-2273 -- two shallow flows spread across the deepest safe gaps

Use the scalar five-unit/three-blocker notation

```text
D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i),                 (1)
```

where `H,q_1,...,q_5` are thirteen-units, `H` is odd, and the strict
first-depth-one blocker profile is

```text
c_1=13u_1,          c_2=13^b u_2,          c_3=13^c u_3,

2<=b<c,                         5<=c<=19,             (2)
```

with all `u_j` thirteen-units. The scalar cover and THM-2198 give

```text
A_0 subset D_(c_1) union D_(c_2) union D_(c_3),

measure(A_0)>=delta_5:=961/6930.                      (3)
```

The theorem removes the deepest owner before measuring multiplicity. This
creates a two-label flow which retains a named exclusion through the common
time `b+1`.

## 1. A parity-sharp guard/deep cap

Put

```text
E_H={x:||Hx||<1/7}.
```

For any deep blocker `q` of exact thirteen-adic valuation `d`, divide the
pair by its gcd and write

```text
A=H/gcd(H,q),             B=q/gcd(H,q)=13^d k.       (4)
```

Then `A` is odd,

```text
gcd(A,k)=1,                    13 does not divide Ak.
```

THM-2080's exact fold is

```text
measure(E_H intersection D_q)
 =2/49+(2/(AB))F(x,y),

x=(A mod 14)/14,             y=(B mod 7)/7,

F(x,y)=min(x,y)+(x+y-1)_+-2xy.                       (5)
```

Since `D_q` has mass `1/7`,

```text
measure(C_H intersection D_q)
 =5/49+13^(-d) R_d(A,k),

R_d(A,k)=-2F(x,y)/(Ak).                              (6)
```

The parity-sharp upper endpoint is

```text
d odd:       R_d(A,k)<=5/49,
             equality exactly at (A,k)=(1,1);

d even:      R_d(A,k)<=5/294,
             equality exactly at (A,k)=(1,6).        (7)
```

Here is the complete proof of the only finite step. The elementary square
bound `F>=-1/8` gives

```text
R_d(A,k)<=1/(4Ak).                                   (8)
```

For odd `d`, every cell with `Ak>=3` is strictly below `5/49`; the two
admissible cells `Ak<=2` have unique maximum `5/49` at `(1,1)`. For even
`d`, every cell with `Ak>=15` is strictly below `5/294`; the complete
twenty-two-cell bank

```text
Ak<=14,       A odd,       gcd(A,k)=1,       13 not|Ak
```

has unique maximum `5/294` at `(1,6)`. Thus define the exact cap

```text
G_d=
  5/49+5/(49*13^d),            d odd,
  5/49+5/(294*13^d),           d even.               (9)
```

Then

```text
measure(C_H intersection D_q)<=G_d.                 (10)
```

This improves the continuous-envelope debit
`5/49+1/(4*13^d)`. Both equality families in (7) are realized after any
admissible common dilation.

## 2. Deepest removal forces exact-one shallow mass

Remove the deepest blocker:

```text
S=A_0 minus D_(c_3).                                 (11)
```

By (3),

```text
S subset D_(c_1) union D_(c_2).                     (12)
```

Define the two shallow exclusive-owner pieces

```text
E_1=S intersection D_(c_1) minus D_(c_2),

E_2=S intersection D_(c_2) minus D_(c_1).           (13)
```

They are disjoint, lie in `A_0`, avoid `D_(c_3)`, and exhaust the part of
`S` with shallow multiplicity one. Consequently

```text
measure(E_1)+measure(E_2)
 >=measure(A_0)
   -measure(C_H intersection D_(c_3))
   -measure(D_(c_1) intersection D_(c_2)).           (14)
```

Let `C_d` denote THM-2263's sharp upper danger-pair endpoint:

```text
C_d=
  1/49+6/(49*13^d),             d even,
  1/49+5/(588*13^d),            d odd.               (15)
```

The shallow valuation gap is `b-1`, so (9), (14), and (15) give

```text
measure(E_1)+measure(E_2)>=L(b,c),

L(b,c)=delta_5-G_c-C_(b-1).                          (16)
```

Over all `150` strict profiles, the unique minimum is

```text
(1,b,c)=(1,3,5),

L(b,c)>=L_0:=5696989/367580070.                      (17)
```

The uniqueness is structural. Among positive gaps, `C_d` has its unique
maximum at the smallest even gap `d=2`, forcing `b=3`. Among `c>=5`, the
correction in `G_c` is uniquely largest at the smallest odd depth `c=5`.
The two equality rows are compatible on one ladder:

```text
H=1,       (c_1,c_2,c_3)=(13,13^3,13^5).            (18)
```

Thus (17) is sharp for the two pair ledgers used in (14), though equality
in the full residual estimate is not claimed.

For comparison, retaining the continuous debit in (8) would give only

```text
558290567/36022846860;
```

the sharp gain in (17) is `29/72773428`.

## 3. The common-time transport minimax

Put

```text
T(x)=13x mod 1,          e_j=measure(E_j).            (19)
```

Because each `E_j` lies in `C_H`, the ten-root guard cap gives

```text
measure(T(E_j))>=(13/10)e_j.                         (20)
```

The first piece is owned by `c_1=13u_1`, so

```text
T(E_1) subset D_(u_1).
```

A unit danger comb occupies at most two roots of a thirteen-root fibre.
Applying that cap at the next step gives

```text
measure(T^2(E_1))>=(169/20)e_1.                     (21)
```

The second piece avoids `D_(c_1)`, hence

```text
T(E_2) subset D_(u_1)^c.
```

The complement of a unit danger comb occupies at most twelve roots.
Therefore

```text
measure(T^2(E_2))>=(169/120)e_2.                    (22)
```

Membership in `D_(c_2)` is constant through the first two images because
`b>=2`. Equations (13) imply

```text
T^2(E_1) subset D_(c_2/169)^c,

T^2(E_2) subset D_(c_2/169).                         (23)
```

The two sets in (23) are disjoint. Circle multiplication is
noncontracting on image mass, so for

```text
Y=T^(b+1)(E_1 union E_2)                             (24)
```

equations (21)--(23) give

```text
measure(Y)>=(169/20)e_1+(169/120)e_2.                (25)
```

There is a second lower bound. At time `b`, the second flow lies in the
unit owner comb:

```text
T^b(E_2) subset D_(u_2).
```

Noncontraction from time two and the final two-root cap give

```text
measure(T^(b+1)(E_2))>=(2197/240)e_2,

measure(Y)>=(2197/240)e_2.                           (26)
```

Minimize the maximum of the right sides of (25)--(26) subject to
`e_1+e_2>=L`. The two affine branches meet at

```text
e_1:e_2=11:12,
```

and give the exact minimax coefficient

```text
measure(Y)>=(2197/460)(e_1+e_2).                    (27)
```

Combining (17) and (27),

```text
measure(Y)>=Y_0:=5696989/76962600

                  =1/14+1397623/538738200

                  =6/91+4357723/538738200.          (28)
```

The coefficient `2197/460` is stronger than choosing a half-mass owner and
applying THM-2255's `169/20` expansion. Improving it from the displayed
local data would require a restriction on the split `(e_1,e_2)` or on the
overlap created after time two.

## 4. The image avoids a named deepest successor

Every point of `E_1 union E_2` avoids `D_(c_3)`. Since `b+1<=c`,
nonmembership transports exactly:

```text
x notin D_(c_3)
 iff T^(b+1)x notin D_s,

s=c_3/13^(b+1)=13^m u_3,
m=c-b-1>=0.                                          (29)
```

Therefore

```text
Y subset D_s^c.                                      (30)
```

This is the named target exclusion missing from the raw one-owner carrier
of THM-2261. It is an exclusion, not a successor cover.

The complement of `D_s` consists, up to null endpoints, of exactly `s`
cyclic safe gaps

```text
I_r=[
 (r+1/14)/s,
 (r+13/14)/s
],                        r in Z/sZ,                 (31)
```

each of length

```text
6/(7s).                                              (32)
```

Let `N_s(Y)` be the number of these gaps whose intersection with `Y` has
positive Haar measure. Equations (28), (30), and (32) give the actual-core
bound

```text
N_s(Y)
 >=ceil((7s/6)Y_0)

 =ceil(
   39878923 * 13^(c-b-1) * u_3
   / 461775600
 ).                                                  (33)
```

In particular:

```text
c=b+1:  at least one gap;

c=b+2:  at least two gaps;

c=b+3:  at least fifteen gaps;

c=b+4:  at least 190 gaps.                           (34)
```

Thus all `135` nonadjacent strict rows `c>=b+2` force a genuine
multi-gap image. This class consists of:

```text
15 rows with b=2,

120 rows with 3<=b<=c-2.                             (35)
```

For the `b=2` boundary, `m=c-3>=2`; hence all fifteen rows force at least
fifteen deepest-successor gaps. This is precisely a boundary on which
THM-2266's centered deep-pair packet degenerates, so the two carriers are
complementary.

For the fifteen adjacent rows `(1,b,b+1)`, the exact profile floor is
slightly stronger. Its unique worst row is `(1,5,6)`, and

```text
measure(Y)
 >=271263857/3501798300

 =1/13+1894757/3501798300.                           (36)
```

Here `s=u_3`. Equation (33) gives the useful finite alternative

```text
u_3<=11,

or Y occupies at least two distinct successor-safe gaps.              (37)
```

The cutoff is sharp for the mass-only argument: at `u_3=11`, one gap has
mass `6/77`, which is still larger than the floor in (36).

## 5. Composition with the bounded-pair atlas

THM-2266 independently treats exactly the `120` rows in the second line of
(35). It forces one of six reduced products, pairing the shallow core
`u_1` with `H` or one of the five `q_i`, to be at most `757`. Combining
the two theorems therefore gives, on those rows,

```text
one of 3,643 ordered primitive pair shapes,

and at least the gap count (33) at the deepest successor.              (38)
```

This conjunction is useful but is not yet a finite interval atlas. The
bounded pair does not involve

```text
s=13^(c-b-1)u_3.
```

Keeping the triggered primitive pair shape fixed while varying `u_3`
leaves its reduced product unchanged but makes both `s` and the gap count
in (33) unbounded. It also changes the cyclic locations of those gaps.
Thus no bound or localization of the deepest gap indices follows from the
number `757` alone.

There is nevertheless a precise finite decision problem once a full
coefficient row is fixed, or once it is taken from the existing finite
rank-eleven/rank-twelve atlases:

1. form the rational endpoint arrangement of `C_H`, the five unit combs,
   and the three blocker combs;
2. extract the interval components of `E_1,E_2`;
3. split their images under `T^(b+1)` at every wrap point;
4. intersect the images with the indexed gaps `I_r` in (31), retaining
   source owner, root sheet, image multiplicity, and positive edge mass;
5. attach the actual successor-owner/carry eligibility sets and compute the
   binary cut energies of THM-2267.

The resulting weighted bipartite object is the **gap-ancestry graph**:
source interval components on one side and deepest safe-gap indices on the
other. Every operation above is exact rational interval arithmetic. A
candidate row fails if this graph admits no legal successor section
consistent with its cover/carry labels.

To make (38) alone into a finite uniform atlas, one still needs a bounded
reduced ratio from `s` to the triggered pair scale, together with bounds or
finite charts for the remaining unit masks. Equivalently, one may execute
the gap-ancestry test directly inside the already proved global finite
rank charts. THM-2273 supplies the compulsory number of gap vertices and
their ancestry labels; it does not supply those missing coefficient bounds.

## 6. Scope, equality, and hostile controls

What is new is the conjunction

```text
positive two-shallow exact-one mass
  + exact time-two owner split
  + common-time image expansion
  + named deepest-comb exclusion
  + quantitative safe-gap spread.                  (39)
```

It addresses all `150` strict rows. It does not address the fifteen
repeated-first profiles `(1,1,c)`.

Three hostile boundaries are load-bearing.

1. The row `(1,3,5)` and ladder (18) show that both pair debits in the
   uniform shallow-mass ledger can be sharp simultaneously.
2. The mass split `e_1:e_2=11:12` makes both branches of the transport
   minimax equal. The local constants alone cannot improve `2197/460`.
3. On adjacent rows with `u_3<=11`, one deepest safe gap has enough
   capacity for the guaranteed image. Thus a uniform two-gap claim there
   would be false at the level of retained information.

Most importantly, (30) is not an image of the original scalar cover. It
says where the transported shallow flow cannot go, but does not name a
family which must cover it after time `b+1`. Therefore neither the
`1/14`, `1/13`, nor multi-gap crossing contradicts the scalar cover.
THM-2273 excludes no valuation profile and LRC(14) remains open. The next
obligation is to attach actual successor eligibility and ancestry to the
gap-ancestry graph.

## 7. Reproduction

Run

```bash
python3 04-computation/lrc14_shallow_owner_deep_gap_spread_thm2273.py
python3 -O 04-computation/lrc14_shallow_owner_deep_gap_spread_thm2273.py

python3 04-computation/lrc14_shallow_owner_deep_gap_spread_referee_thm2273.py
python3 -O 04-computation/lrc14_shallow_owner_deep_gap_spread_referee_thm2273.py
```

The primary companion checks both parity-sharp guard/deep banks, every
strict profile, the unique uniform and boundary minimizers, equality
compatibility, the exact transport minimax, and all profile-gap counts.
The independent referee does not call either folded formula: it constructs
the relevant combs as exact rational interval unions, intersects both
finite endpoint banks directly, verifies the root laws on all `169`
endpoint-free parent phases, and checks the safe-gap geometry. Normal and
optimized transcripts are byte-identical to their stored outputs. QED.
