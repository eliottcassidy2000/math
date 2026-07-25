---
id: THM-2255
title: "Valuation-separated pair cap and exclusive-owner mass"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY REFEREED. If two danger-comb
  speeds have unequal 13-adic valuations, their intersection has Haar mass
  at most 25/1183; equality holds exactly when their reduced unordered ratio
  is 1:169. The proof combines THM-1166's folded formula, an analytic
  ab>=346 tail, and a complete 70-pair coprime bank for ab<=345, independently
  replayed by exact rational interval intersection. Consequently every one
  of the 165 remaining first-depth-one scalar profiles has a quantitative
  multiplicity-one blocker-owner stratum. For the 150 strict rows (1,b,c),
  1<b<c, its total mass is at least 88159/1171170 and some labelled owner
  has mass at least 88159/3513510. For the 15 rows (1,1,c), the corresponding
  floors are 14627/585585 and 14627/1756755. Transporting the guaranteed
  labelled stratum through its owner's expiration gives expansion by at
  least 169/20; in the strict class the resulting image has mass at least
  88159/415800>1/7 and therefore cannot fit in one successor danger comb.
  This forces a genuine post-expiration split/handoff in any consumer, but
  supplies no such covering law itself and excludes no valuation profile.
source: codex-2026-07-25-depth-one-exclusive-owner
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-1191-four-comb-pair-floor-refuted-by-thirteen-adic-ladder
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
script:
  - 04-computation/lrc14_valuation_separated_pair_owner_mass_thm2255.py
  - 04-computation/lrc14_valuation_separated_pair_owner_mass_referee_thm2255.py
output:
  - 05-knowledge/results/lrc14_valuation_separated_pair_owner_mass_thm2255.out
  - 05-knowledge/results/lrc14_valuation_separated_pair_owner_mass_referee_thm2255.out
script_sha256:
  - 3fdef0f1760e06e3590f509bb2bf648931ec1efa8904a0f99c51bc2637d6cf3f
  - 46326deabb035a0c0d4a3ef2e47a974e06ef3916deaaaf9dc28abea357a44f81
output_sha256:
  - 1ad0bcac88b9dab2b415d7556c80274425d70aa82abd2116c4120137317ffb65
  - dbd93cd6460d6cc0b9cd64c167389cbcd1708f433432b8d0f1f404f790b35f47
hash_basis: working-tree bytes (LF)
---

# THM-2255 -- valuation separation forces exclusive owner mass

Put

```text
D_s={x in R/Z:||sx||<1/14},
rho(s,t)=measure(D_s intersection D_t).               (1)
```

There are two conclusions. First, unequal thirteen-adic valuations force
the sharp pair cap

```text
nu_13(s)!=nu_13(t)  implies  rho(s,t)<=25/1183.       (2)
```

Second, applying (2) to the three actual blockers in the remaining
first-depth-one scalar branch produces a uniform mass of points having
exactly one blocker owner. This is stronger than the merely topological
private-owner conclusion and lives on the positive residual itself.

## 1. The sharp ramified pair cap

Divide `s,t` by their gcd and write the reduced unordered pair as

```text
1<=a<b,                    gcd(a,b)=1.                (3)
```

If the original valuations are unequal, exactly one of `a,b` is divisible
by thirteen. THM-1166 gives the exact folded formula

```text
rho(a,b)
 =1/49+[F((a+b) mod 14)-F((b-a) mod 14)]/(196ab),

F(r)=r(14-r),                    0<=r<14.             (4)
```

Since `0<=F<=49`, equation (4) implies

```text
rho(a,b)<=1/49+1/(4ab).                              (5)
```

For `ab>=346`, the right side is at most

```text
1/49+1/(4*346)
 =1433/67816
 <25/1183.                                           (6)
```

It remains to inspect the finite universe

```text
1<=a<b,  gcd(a,b)=1,  ab<=345,
exactly one of a,b divisible by thirteen.             (7)
```

There are exactly `70` pairs in (7). Exact evaluation of (4) gives

```text
max_(7) rho(a,b)=25/1183,                             (8)
```

with the unique reduced equality pair

```text
(a,b)=(1,169).                                       (9)
```

An independent referee constructs `D_a,D_b` as unions of exact rational
intervals and intersects them directly for all `70` pairs. It obtains the
same maximum and equality locus without calling (4).

At the equality pair, THM-1191's thirteen-adic ladder formula also gives
directly

```text
rho(1,13^2)
 =1/49+(6/49)13^(-2)
 =25/1183.                                           (10)
```

Common dilation preserves Haar pair mass. Thus equality in (2) holds
exactly for the unordered speed pairs

```text
{s,t}={g,169g},                       g>=1.           (11)
```

This proves (2), including its sharp boundary.

## 2. The positive residual and its owner multiplicity

Use THM-2198's scalar `5+3` notation:

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

c_j=13^(lambda_j)u_j,              13 does not divide u_j,

1=lambda_1<=lambda_2<lambda_3<=19.                  (12)
```

The three actual blockers `c_1,c_2,c_3` are distinct, and the scalar cover
gives

```text
A_0 subset D_(c_1) union D_(c_2) union D_(c_3)       (13)
```

almost everywhere. The sharp five-unit residual floor is

```text
measure(A_0)>=delta_5:=961/6930.                     (14)
```

For each labelled blocker, define its exclusive-owner stratum

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(k!=j)D_(c_k).                    (15)
```

The `E_j` are pairwise disjoint. By (13), their complement in `A_0`
consists of points having at least two blocker owners. Hence

```text
A_0 minus union_j E_j
 subset union_(1<=i<j<=3)(D_(c_i) intersection D_(c_j)), (16)
```

and therefore

```text
sum_j measure(E_j)
 >=delta_5-sum_(i<j)rho(c_i,c_j).                    (17)
```

This is the key change of object. The pair caps are not being used as an
unconditional Bellman penalty. They are applied on the actual positive
residual to quantify the part where the blocker-owner word has
multiplicity exactly one.

## 3. The 150 strict-depth profiles

Suppose

```text
1<lambda_2<lambda_3.                                 (18)
```

Every blocker pair then has unequal thirteen-adic valuations. Equations
(2), (14), and (17) give

```text
sum_j measure(E_j)
 >=961/6930-3(25/1183)
 =88159/1171170
 =0.0752742983512... .                               (19)
```

By the labelled pigeonhole principle, some `j` satisfies

```text
measure(E_j)
 >=88159/3513510
 =0.0250914327837... .                               (20)
```

The current `165`-row first-depth-one ledger has exactly `150` rows of
this type:

```text
(1,b,c),        5<=c<=19,        2<=b<c.             (21)
```

No equality relation among the normalized cores `u_j` was used.

## 4. The 15 repeated-first-depth profiles

Now suppose

```text
lambda_1=lambda_2=1<lambda_3.                         (22)
```

The pair `(c_1,c_2)` consists of distinct speeds but has equal valuation, so
THM-1166's generic distinct-pair cap gives

```text
rho(c_1,c_2)<=1/14.                                  (23)
```

The other two pairs have unequal valuations and obey (2). Thus

```text
sum_j measure(E_j)
 >=961/6930-1/14-2(25/1183)
 =14627/585585
 =0.0249784403631... .                               (24)
```

Some labelled owner satisfies

```text
measure(E_j)
 >=14627/1756755
 =0.00832614678769... .                              (25)
```

There are exactly `15` rows of this type:

```text
(1,1,c),                       5<=c<=19.             (26)
```

Equations (21) and (26) give the exact split

```text
150+15=165.                                          (27)
```

## 5. The expansion-at-expiration observable

The labelled floor has a concrete dynamical consumer. Put

```text
T(x)=13x mod 1.                                      (28)
```

Fix a label `j` supplied by (20) or (25), and write

```text
lambda=lambda_j,                 c_j=13^lambda u_j.  (29)
```

Because `E_j subset A_0 subset C_H`, the exact guard-root count gives

```text
#{x:T(x)=y, x in E_j}<=10
```

almost everywhere. Haar disintegration therefore yields

```text
measure(T(E_j))>=(13/10)measure(E_j).                 (30)
```

Every multiplication map is noncontracting on Borel image mass:

```text
measure(T(B))>=measure(B).                            (31)
```

Exact division support gives

```text
T^lambda(E_j) subset D_(u_j).                         (32)
```

Finally, a unit danger comb occupies at most two roots of each
thirteen-root fibre:

```text
#{x:T(x)=y, x in D_(u_j)}
 =2-1_(D_(u_j))(y)
 <=2.                                                (33)
```

Apply (33) to `B=T^lambda(E_j)` and combine it with (30)--(31):

```text
measure(T^(lambda+1)(E_j))
 >=(13/2)measure(T^lambda(E_j))
 >=(169/20)measure(E_j).                             (34)
```

For every strict-depth profile, (20) and (34) imply

```text
measure(T^(lambda+1)(E_j))
 >=88159/415800
 =1/7+28759/415800
 =0.212022607023... .                                (35)
```

Every danger comb has mass exactly `1/7`. Consequently the marked image in
(35) cannot be contained, even almost everywhere, in any single successor
danger comb. A consumer which transports this exclusive-owner flow through
its owner's expiration must exhibit a genuine split, collision, or handoff.

For the repeated-first-depth class, the same argument gives the still
subcritical but uniform floor

```text
measure(T^(lambda+1)(E_j))
 >=14627/207900
 =1/7-15073/207900
 =0.0703559403560... .                               (36)
```

Thus the two valuation classes have a sharp qualitative difference at the
one-comb threshold.

There is also exact owner-word persistence before expiration. If two labels
`i,j` distinguish two exclusive strata, then their membership and
nonmembership bits transport without loss through every level

```text
r<=min(lambda_i,lambda_j),                            (37)
```

because

```text
x in D_(13^r v) iff T^r(x) in D_v.                   (38)
```

Their images are therefore disjoint through that common depth. An
“owner switch” before expiration is not the right object; the first place
where information can be lost is the unit-owner fibre at expiration.

## 6. Scope and tournament reading

The theorem proves a large exclusive-owner carrier and a forced
post-expiration spill. It does **not** prove that the spill is covered by a
named family of successor combs. Equation (35) therefore does not contradict
the scalar cover and excludes none of the `165` profiles.

THM-2257 clarifies the exact missing implication by positive contrast. Its
depth-three image lower bound closes because a transported negative-support
inclusion puts the whole image inside a named union of known mass. Equation
(35) supplies the analogous lower image in depth one, but there is presently
no valid post-expiration inclusion. The next theorem must create that target
set; comparing the two image sizes without it would repeat an invalid
image-of-cover transfer.

The honest finite relation has the three labelled owner strata as vertices.
In the strict class, orient `i -> j` when `lambda_i<lambda_j`; this valuation
tournament is transitive. Its useful sidecar is not the orientation alone
but the vertex mass and the expiration time. Equation (35) says that at
least one weighted vertex produces more than one comb of image mass when it
expires. The next proof obligation is to retain which root sheets receive
that excess and whether two owner flows collide there. A tournament on bare
blocker labels loses exactly this information.

The equality control `(1,169)` is also structural. It is the even-gap
positive-correlation edge of THM-1191's thirteen-adic ladder, so the
constant `25/1183` cannot be lowered by another pairwise argument. Any
improvement to (19) must use compatibility among all three blocker pairs,
the five-unit residual, or the labelled image/carry geometry.

## 7. Exact verification

Run

```bash
python3 04-computation/lrc14_valuation_separated_pair_owner_mass_thm2255.py
python3 -O 04-computation/lrc14_valuation_separated_pair_owner_mass_thm2255.py

python3 04-computation/lrc14_valuation_separated_pair_owner_mass_referee_thm2255.py
python3 -O 04-computation/lrc14_valuation_separated_pair_owner_mass_referee_thm2255.py
```

The primary companion checks the analytic tail comparison, all `70` pairs
in (7), the unique equality class, common-dilate hostile controls, and every
fraction in (19)--(25). The independent referee constructs and intersects
exact rational interval unions for the whole finite bank and several
unreduced equality controls. Normal and optimized transcripts are
byte-identical for both paths. QED.
