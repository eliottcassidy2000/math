---
id: THM-2561
title: "Positive-interval physical six-comb blind-necklace hostile"
status: >
  PROVED + VERIFIED-EXACT.  On the direct F_13 predecessor chart, the
  primitive physical coefficients H=183 and
  (q_1,...,q_5)=(95,93,114,198,304) realize THM-2558's blind safe mask
  {0,1,4} throughout the open base interval
  (13/1281,29/2772), of exact mass 53/169092.  Their one guard-failure
  mask and five ordinary danger masks are pairwise disjoint and partition
  the ten empty roots.  Adding blockers (13,169,13^5) realizes the strict
  profile (1,2,5) and makes the same mask a positive exclusive-owner fibre
  on a nested interval of mass 6/199927.  With k_a=95 and q_*=k_b=93,
  root 3 carries an
  isolated literal old k_a failure while q_* and all four neutral roles
  are safe, yet every one of the twelve lexicographic necklace selectors
  misses root 3.  Thus physical six-comb incidence, q_* safety, and removal
  of cofailure do not eliminate the 202 blind necklaces; a global
  scalar-cover/deep-blocker input or a different selector is necessary.
  The target-informed chord 0->3 of slope 3 does see the isolated failure
  and has Cayley gradient -1, so the hostile is specific to the
  lexicographic selector and does not obstruct the target-informed bypass.
  No global scalar cover, later target root, row exclusion, or LRC(14)
  conclusion follows.
source: radon-kakeya-2026-07-27-physical-blind-necklace
depends_on:
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2558-sparse-owner-fibre-all-slope-target-role-forcing-boundary
related:
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2559-target-informed-chord-and-universal-old-repair-packet
script: 04-computation/lrc14_physical_blind_necklace_hostile_thm2561.py
output: 05-knowledge/results/lrc14_physical_blind_necklace_hostile_thm2561.out
script_sha256: 1a5975d5981750e29165bd0a4b390269759b9466f50ff4f9fed08035d3962632
output_sha256: cb23ef65b9e86706a4141248804b8866125b01bfe219712e7d45765e8c110cf8
hash_basis: working-tree bytes (LF)
---

# THM-2561 -- the blind necklace occurs on a positive physical interval

**PROVED + VERIFIED-EXACT.**

THM-2558 leaves two possible readings of its 202 blind necklaces:

```text
abstract support hostile only,

or

an obstruction which survives the actual guard/unit comb geometry.       (1)
```

The second reading is correct.  In fact the canonical three-root blind
mask occurs on a positive interval with a stronger sidecar: the missed
target-active failure is isolated, while the designated deepest graft role
is safe.  The obstruction is therefore not caused by arbitrary Boolean
masks, a wall endpoint, a neutral cofailure, or loss of `q_*` safety.

## 1. Direct-fibre failure masks and their walls

Put `p=13` and use the direct predecessor chart

```text
iota_r(z)=(z+r)/13,                r in F_13, z in [0,1).    (2)
```

For an ordinary role `w` and a guard `H`, write

```text
D_w={x:||wx||<1/14},

E_H={x:||Hx||<1/7},                C_H=E_H^c a.e.,           (3)

A_0=C_H minus union_(i=1)^5 D_(q_i).                         (4)
```

The choice of strict or weak membership at the finitely many walls is
irrelevant below because the base interval will be open and wall-free.

It is useful to treat an ordinary radius as `a=1` and the guard radius as
`a=2`.  The role-failure mask is

```text
F_(w,a)(z)
 ={r:||w(z+r)/13||<a/14}.                                   (5)
```

It can change only when, for some integer `m`, sign, and root `r`,

```text
w(z+r)/13=m plus_or_minus a/14,

z=(13/w)(m plus_or_minus a/14)-r.                            (6)
```

Thus constancy on a proposed base interval is a finite exact wall
calculation, not a numerical continuity inference.

## 2. A primitive physical realization of the blind necklace

Choose

```text
H=183,

(q_1,q_2,q_3,q_4,q_5)=(95,93,114,198,304),                  (7)

z_0=1/97.
```

The six coefficients are positive, pairwise distinct thirteen-units with
distinct residues

```text
(1,4,2,10,3,5) mod 13.                                     (8)
```

The guard is odd and the common gcd is one.  Direct exact evaluation of
(5) at `z_0` gives

| role | physical failure roots |
|---|---|
| guard `H=183` | `{10,11,12}` |
| `q_1=95` | `{3}` |
| `q_2=93` | `{6}` |
| `q_3=114` | `{5,9}` |
| `q_4=198` | `{8}` |
| `q_5=304` | `{2,7}` |

These six sets are pairwise disjoint.  Their union is

```text
{2,3,5,6,7,8,9,10,11,12},                                  (9)
```

so (4) has the exact root mask

```text
e=A_0={0,1,4}.                                               (10)
```

This is not merely a point-fibre realization.  Enumerating the wall set
(6) for all six roles and all thirteen branches, the nearest wall below
`z_0` is

```text
z=13/1281,             H iota_0(z)=1/7,                     (11)
```

and the nearest wall above it is

```text
z=29/2772,             198 iota_12(z)=183-1/14.             (12)
```

There is no other role wall between (11) and (12).  Therefore the entire
table and (10) hold for every

```text
z in I=(13/1281,29/2772).                                   (13)
```

The exact masses are

```text
|I|=29/2772-13/1281=53/169092>0,

|iota_r(I)|=53/2198196                 for every root r.     (14)
```

The companion referee enumerates all `1,974` distinct-per-role wall values
used in this comparison and performs every membership check with rational
arithmetic.

The hostile also sits inside a genuine exclusive-owner cell of an allowed
valuation profile.  Add the blocker coefficients

```text
(c_1,c_2,c_3)=(13,169,13^5)=(13,169,371293),

(nu_13(c_1),nu_13(c_2),nu_13(c_3))=(1,2,5).               (14a)
```

Every blocker status is root-independent in (2).  Around `z_0`, the nearest
blocker walls come from `c_3` and cut out the nested interval

```text
J=(4117/399854,4129/399854) subset I.                      (14b)
```

For every `z in J`, the exact blocker word is

```text
D_(c_1)(iota_r(z))=1,

D_(c_2)(iota_r(z))=D_(c_3)(iota_r(z))=0
                                      for every r.          (14c)
```

Indeed the two endpoint values of the reduced deepest blocker are
`13^4z=4117/14` and `4129/14`, the consecutive walls surrounding
`13^4z_0=28561/97`.  Hence THM-2296's exclusive-owner event obeys

```text
E_1=A_0 intersection D_(c_1)
       minus (D_(c_2) union D_(c_3))
   =A_0

on J,                    |J|=6/199927>0.                    (14d)
```

Thus the blind mask is the literal direct-root mask of a positive
exclusive-owner fibre in the strict `(1,2,5)` profile.  This still does not
assert that the three blockers cover `A_0` globally.

## 3. The target-active failure is isolated and every selector misses it

At the local owner-pivot interface, designate the two ordinary roles

```text
k_a=95,                    q_*=k_b=93.                       (15)
```

This is the THM-2350/2461 target-role convention: `k_a` is the target-active
role remaining after the deepest graft `q_*=k_b` is fixed.  We are not yet
asserting the global scalar-cover hypotheses.  On every fibre above (13),
root `3` satisfies

```text
3 in D_(k_a),

3 notin D_(q_*),

3 notin E_H union D_114 union D_198 union D_304.             (16)
```

Thus the old `k_a` failure at root `3` is singleton, `q_*` is safe there,
and none of the four target-neutral roles cofails there.  This is stronger
than the support-only hostile in THM-2558.

For the mask (10), direct evaluation of THM-2558's lexicographic necklace
selector gives, in slopes `tau=1,...,12`,

```text
(t_tau)=(2,2,7,8,9,7,7,9,9,11,2,12),                        (17)

H(e)={2,7,8,9,11,12}.                                       (18)
```

In particular

```text
3 notin H(e).                                                (19)
```

Equations (13)--(19) prove a positive-interval physical hostile: every
lexicographic all-slope selected-head packet misses the isolated old
target-active failure, even though its branch has the positive mass in
(14).  Physical six-comb incidence, the one/two-root ordinary capacity,
the three/four-root guard capacity, `q_*` safety, and absence of neutral
cofailure are therefore insufficient to eliminate THM-2558's blind layer.

## 4. A target-informed chord bypasses the hostile locally

The obstruction belongs to the **lexicographic choice rule**, not to the
complete directed boundary graph.  In (10), choose

```text
s=0 in e,                    t=3 in F_(k_a)(z),

tau=t-s=3.                                                   (20)
```

This is an occupied-to-empty chord aimed at the known target failure.  The
THM-2537 Cayley identity gives pointwise

```text
(C_3 e)_0+(C_3 e)_3=e_3-e_0=-1.                             (21)
```

After multiplication by any nonnegative root-invariant owner weight, the
negative two-endpoint sum again scalarizes its positive mass.  Hence a
target-informed chord selector can bypass this blind necklace locally.
THM-2559 proves the universal measurable version on all `165` rows; it is a
subsequent related theorem, not a dependency here.

Equation (21) remains entirely at the old predecessor horizon.  It neither
places `k_a` genuinely later nor identifies its old ancestry root with a
future semantic root.

## 5. Exact stopping boundary

The example includes three correctly typed blockers and a positive exclusive
owner fibre, but it does **not** prove their global cover of `A_0`, select a
terminal owner word, or produce a live scalar row.  It therefore does not
prove that the blind interval survives every global counterexample constraint.
What it proves is the sharp logical boundary

```text
local physical six-comb geometry
+ positive strict-profile exclusive-owner fibre
+ isolated k_a failure
+ q_* safety
+ no neutral cofailure
+ all twelve lexicographic slopes
  -/-> selected target-active head.                          (22)
```

Consequently the first branch of THM-2558's remaining program must use
information not present in a single six-comb fibre: the global scalar cover,
deep-blocker/owner coherence across base cells, or a target-informed selector
such as (20).  No later-root map, ancestry/future-digit intertwiner, scalar-row
exclusion, or LRC(14) conclusion follows.

## Exact referee

Run

```bash
python3 04-computation/lrc14_physical_blind_necklace_hostile_thm2561.py
python3 -O 04-computation/lrc14_physical_blind_necklace_hostile_thm2561.py
```

Both runs byte-match

```text
05-knowledge/results/lrc14_physical_blind_necklace_hostile_thm2561.out
```

and the recorded hashes.  The referee checks coefficient admissibility,
enumerates the exact wall cell, reconstructs all six physical masks and
their disjoint partition, verifies the `(1,2,5)` blocker word and nested
positive exclusive-owner interval, checks isolated `k_a` failure and `q_*`
safety, independently computes all twelve selected heads, and checks the
target-informed Cayley chord.

**QED.**
