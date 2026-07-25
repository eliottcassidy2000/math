---
id: THM-2234
title: "First-depth-one private two-owner mass and two-step image expansion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch,
  peeling any actual blocker of exact 13-adic depth one leaves a residual
  covered by the other two blockers with Haar mass at least 2593/90090.
  The improvement over the generic six-tooth Hunter floor comes from the
  exact odd-guard charge restriction: the peeled blocker has reduced
  denominator divisible by thirteen, whose worst negative charge is 5/637
  at (1,13). The private residual expands by 13/10 on its first image; that
  image avoids the peeled unit comb and therefore expands by another 13/12
  on its next image. The resulting exact floors are 2593/69300 and
  33709/831600. Both fibre caps are locally sharp, so a further uniform
  gain needs labelled winding, owner/carry information, or another
  post-deletion constraint. This strengthens the low-first-depth scalar
  carrier but excludes no valuation profile by itself and does not prove
  LRC(14).
source: klein-2026-07-24-first-depth-one-private-mass
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2224-transfer-owner-word-temporal-union-bound
  - THM-2226-three-checkpoint-bellman-sieve-and-eight-profile-residue
script: 04-computation/lrc14_first_depth_one_private_mass_thm2234.py
output: 05-knowledge/results/lrc14_first_depth_one_private_mass_thm2234.out
script_sha256: 7a9c0b842b38c81b619c6fa505ba526b6539c18b2c75a2d77b88f96542c5e46e
output_sha256: 9ea42461e17f79c6cf2cc726e58df413b9d032e7bc40f3b357e2645e53ce70f3
hash_basis: working-tree bytes (LF)
---

# THM-2234 -- the depth-one private two-owner remainder

Retain the scalar-cover notation

```text
D_b={x in R/Z:||bx||<1/14},
C_H={x in R/Z:||Hx||>1/7}.                           (1)
```

Suppose `H,q_1,...,q_5` are positive thirteen-units, `H` is odd, the
`q_i` are distinct, and the three distinct positive actual blockers
`c_1,c_2,c_3` are divisible by thirteen. Assume

```text
C_H subset union_(i=1)^5 D_(q_i)
             union D_(c_1) union D_(c_2) union D_(c_3)     (2)
```

almost everywhere. This is THM-2198's scalar `5+3` branch. Suppose one
displayed blocker has exact depth one and write

```text
c_1=13a,                         13 does not divide a. (3)
```

Define the unit residual and its piece private from `c_1` by

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

P=A_0 minus D_(c_1)
 =C_H minus (union_(i=1)^5D_(q_i) union D_(c_1)).     (4)
```

Then the following statements hold.

```text
measure(P)>=2593/90090,                              (5)

P subset D_(c_2) union D_(c_3)             a.e.,    (6)

measure(T(P))>=2593/69300,                           (7)

measure(T^2(P))>=33709/831600,                       (8)
```

where `T(x)=13x mod 1`. Moreover

```text
T(P) subset T(A_0) minus D_a,                        (9)

T(P) subset D_(c_2/13) union D_(c_3/13)    a.e.     (10)
```

If both remaining blockers have depth at least two, then also

```text
T^2(P)
 subset D_(c_2/13^2) union D_(c_3/13^2)    a.e.     (11)
```

The first assertion supplies the quantitative private-owner floor that was
left only topologically positive in THM-2198. The second image supplies a
new strict improvement over that theorem's later noncontraction floor.

## 1. The specialized odd-guard charge ledger

Put

```text
E_H={x:||Hx||<1/7}.
```

For any positive coefficient `r`, reduce

```text
A=H/gcd(H,r),                 B=r/gcd(H,r).           (12)
```

Then `A` is odd. THM-2080 writes

```text
measure(E_H intersection D_r)=2/49+e(A,B)             (13)
```

and proves

```text
e(A,B)>=-1/(4AB).                                    (14)
```

THM-2137 records the five worst distinct negative charges available to the
five `q_i`:

```text
5/294, 8/539, 3/245, 3/245, 4/441.                  (15)
```

The depth-one blocker has a stricter row type. Since `13` does not divide
`H` but `13` divides `c_1`, its reduced denominator in (12) satisfies

```text
13 divides B,                 13 does not divide A.   (16)
```

The sharp negative-charge bound on this subfamily is

```text
-e(A,B)<=5/637,                                      (17)
```

with equality at `(A,B)=(1,13)`.

Indeed, if `AB>=32`, equation (14) gives

```text
-e(A,B)<=1/128<5/637.                                (18)
```

For `AB<=31`, coprimality, odd `A`, and `13|B` leave only

```text
(A,B)=(1,13),(1,26).                                 (19)
```

Direct substitution in THM-2080's exact fold gives respectively

```text
-e(1,13)=5/637,             -e(1,26)=3/1274.         (20)
```

This proves (17).

Apply Hunter's union bound with the star centred at `E_H` to the six danger
combs in (4). Since `measure(E_H)=2/7` and every danger comb has measure
`1/7`, equations (13), (15), and (17) give

```text
measure(P)
 >=5/49
   -(5/294+8/539+3/245+3/245+4/441)
   -5/637
 =2593/90090.                                        (21)
```

The modular debit in (21) is not an incompatible collection of pair rows.
All six worst rows occur simultaneously for

```text
H=99,
(q_1,...,q_5)=(594,9,165,495,22),
c_1=1287.                                            (22)
```

Their reduced pairs are

```text
(1,6),(11,1),(3,5),(1,5),(9,2),(1,13).              (23)
```

Thus the independent-charge ledger itself cannot improve (21). This does
not assert equality in Hunter's full union bound, where higher intersections
remain.

Finally, (6) follows directly from the cover: on `P`, neither a unit tooth
nor `c_1` is active, so one of `c_2,c_3` must own the point. Images of the
null exceptional set remain null.

## 2. The first image: the guard cap

For a Borel set `B`, define the root occupancy

```text
n_B(y)=#{x:T(x)=y and x in B}.                       (24)
```

Haar disintegration gives

```text
13 measure(B)=integral n_B(y)dy.                     (25)
```

Because `P subset C_H`, THM-2198's guard-root count gives

```text
n_P(y)<=10                                           (26)
```

away from finitely many endpoints. Its essential support is `T(P)`, so

```text
13 measure(P)<=10 measure(T(P)).                     (27)
```

Equations (21) and (27) prove

```text
measure(T(P))
 >=(13/10)(2593/90090)
 =2593/69300,                                        (28)
```

which is (7).

The shallow tooth is constant on every root fibre:

```text
1_(D_(13a))(x)=1_(D_a)(T(x)).                        (29)
```

Consequently a point in `T(P)` cannot lie in `D_a`, proving (9). For the
same reason, (6) pushes forward to (10). We use only these inclusions; no
generic image-of-a-set-difference rule is needed.

## 3. The second image: the peeled-unit complement cap

For a thirteen-unit `a`, THM-2198's exact root count is

```text
#{x:T(x)=y and x in D_a}=2-1_(D_a)(y)                (30)
```

almost everywhere. Hence the complementary root count is

```text
#{x:T(x)=y and x notin D_a}
 =11+1_(D_a)(y)
 <=12.                                               (31)
```

Put `E=T(P)`. Equation (9) says `E subset D_a^c`. Applying (25) and (31)
to `E` gives

```text
13 measure(E)<=12 measure(T(E)).                     (32)
```

Combining (28) and (32) yields

```text
measure(T^2(P))
 >=(13/12)(2593/69300)
 =33709/831600,                                      (33)
```

proving (8). If `13^2` divides each of `c_2,c_3`, pushing (6) through two
images also proves (11).

There is a useful strict-ladder restatement. If

```text
1=lambda_1<lambda_2<lambda_3
```

and `R_1=T(A_0)\D_a` is THM-2198's first peeled remainder, then (9) gives
`T(P) subset R_1`. Noncontraction of circle multiplication and (33) imply

```text
measure(T^u(R_1))>=33709/831600
                 for every 1<=u<=lambda_2-1.         (34)
```

In particular, THM-2198's endpoint remainder `R_2` at the second valuation
has the improved lower floor (33), replacing `2593/69300` there.

## 4. Both root caps are sharp

Neither `10` in (26) nor `12` in (31) can be lowered from the displayed
local information.

For the first cap, work modulo `13^3=2197` over the phase

```text
y=1/169,
H=1,          (q_1,...,q_5)=(1,2,3,4,5),
a=14.                                                  (35)
```

The ten guard sheets are `2,3,...,11`. Every one of the five distinct unit
masks is the singleton sheet `0`, while `y notin D_a`; hence all ten guard
sheets lie in the private residual.

For the second cap, take

```text
H=a=1,          (q_1,...,q_5)=(1,2,3,4,5),           (36)
```

again on the `13^3` torsion grid, with `c_1=13`. The image `T(P)` contains
all twelve nonzero roots

```text
1/13,2/13,...,12/13                                  (37)
```

of zero under `T`, and omits the zero root. Thus one second-stage fibre has
occupancy exactly twelve. All inequalities in these controls are strict
because a power of thirteen is coprime to seven and fourteen.

The companion checks the charge censuses, the compatible packet (22), every
fraction in (21), (28), and (33), and both torsion controls using exact
integer arithmetic. Reproduce with

```bash
python3 04-computation/lrc14_first_depth_one_private_mass_thm2234.py
python3 -O 04-computation/lrc14_first_depth_one_private_mass_thm2234.py
```

Normal and optimized transcripts are byte-identical.

## 5. Scope and next obstruction

The new second-image floor is still far below the uniform two-comb capacity

```text
measure(D_u union D_v)<=25/91.                       (38)
```

Therefore (33) alone excludes no valuation profile. The sharp controls show
why another scalar occupancy cap is not the missing move. A continuation
must retain at least one of:

```text
the labelled anti-defect root,
the owner of the two-deep union,
cross-core phase/carry,
or a genuinely new post-deletion floor after the second owner expires. (39)
```

This is a quantitative strengthening of the depth-one lane, not a closure
of the scalar branch or of LRC(14). QED.
