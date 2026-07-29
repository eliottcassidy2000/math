---
id: THM-2903
title: "One-hard actual H3 link and bad-triangle closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  All 61 non-direct one-hard roots from
  THM-2901 close: 1,907 of 1,961 actual H3 pair children close directly,
  the 54 unresolved edges have only two triangles, and exact literal
  two-cover caps close both triangle children.
source: codex/lrc-j6-parity-seed-audit-2026-07-29
depends_on:
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2902-omission-six-ranked-h1-depth-one-closure
related:
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2898-unique-max-gate-five-parity-matching-closure
verification:
  - 04-computation/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.ledger.out
---

# THM-2903 -- one-hard actual H3 link and bad-triangle closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

THM-2901 leaves `61` roots having exactly one scalar-hard marked suffix
and failing its direct `q_5+2beta_2<h` certificate.  Every one of those
`61` unique hard suffixes admits no five-cover.  Hence all `61` bodies are
whole seven-body/six-slot root closures.

The set of `61` bodies is pinned by the canonical digest

```text
f58c4143f329d215ff9bb7ec594d172e831fbccbab713a2398dc6cb53c60b8b7.
```

Its full sorted list is in the verification note.

## 2. Actual H3 cores

Fix one of the `61` literal hard carriers `C`, with mass `h`, ordered
excluded gate prefix `P`, and exact global pair cap `beta_2` from
THM-2901.  Since

```text
beta_2<4h/7,
```

THM-2893 forces at least three labels of every hypothetical five-cover
into

```text
H_3={w allowed : c_C(w)>=(h-beta_2)/3}.                   (1)
```

The discrepancy cutoff is scanned exactly and gives

```text
roots                                                61
cutoff labels scanned                         at most 920
sum of actual |H_3|                              461
actual |H_3| range                              2..23
actual unordered H_3 pairs                     1,961.     (2)
```

Thus the theorem-bearing flag universe is four orders of magnitude below
the cutoff-binomial upper bound that motivated the census.

## 3. Literal singleton links

For an actual pair `L subset H_3`, put

```text
C_L=C minus D_L,       h_L=|C_L|.
```

Every three-subset of a hypothetical five-cover is parent-heavy: its
complementary pair has coverage at most `beta_2`.  Therefore each of the
three child labels `z` satisfies the exact THM-2893 link law

```text
c_(C_L)(z)
 =U_C(L union {z})-U_C(L)
 >=h-beta_2-U_C(L)
 =h_L-beta_2.                                            (3)
```

All child labels consequently lie in

```text
J_L={z allowed outside L:c_(C_L)(z)>=h_L-beta_2}.         (4)
```

If

```text
delta_L=6h_L/7-beta_2>0,
gamma_L=(99/70)r_L/7,
```

then strict discrepancy seals `(4)` at

```text
N_L=ceil(gamma_L/delta_L)-1.                              (5)
```

Exactly `1,908` pairs satisfy `delta_L>0`; the maximum cutoff is `60,333`,
the maximum actual link size is `1,260`, and the sum of the actual link
sizes is `10,104`.  The remaining `53` pairs are recorded as cutoff
failures, not as cover witnesses.

The exact restricted child tests give

```text
|J_L|<3                                      1,409
top-three singleton sum <h_L                   484
q_3(J_L)+B_2(J_L)<h_L                           14
finite-link survivor                              1
delta_L<=0                                       53
                                                -----
                                                1,961.   (6)
```

The middle pair cap ranges only over distinct vertices in `J_L`; every
winning union is literal.  Thus `1,907` pair flags are universally
resolved.  The unresolved graph has `54` edges on ten roots.  Fifty-one
roots already have no unresolved edge; their sorted-list digest is

```text
e2a7f23111c3dc8b41a6518fa83571c3a334b1f601e8216ee5eb0a701539c92a.
```

## 4. Unresolved-pair triangle obstruction

Let `B` be the undirected graph on `H_3` whose edges are the `54`
unresolved pairs.  A hypothetical cover supplies at least three vertices
of `H_3`.  Every pair among any chosen three must be unresolved, since a
resolved pair cannot belong to a cover.  Hence `B` must contain a
triangle, exactly the `ell=2,m=3` clique obstruction of THM-2893.

Nine of the ten nonempty bad graphs are triangle-free.  The last body is

```text
E=(1,3,6,9,11,12,13),       rank 1,       apex 21,       P=(21).
```

Its `21` unresolved edges have exactly two triangles:

```text
(15,24,33),       (24,30,33).                             (7)
```

There are no other unresolved triangles.  Their canonical digest is

```text
635acf7fbc9a23257c9a746fc5f5ed1eb4db4538baff7ffe42632cfbd678ece3.
```

This is the decisive compression: it is enough to close the two literal
two-label children behind `(7)`, not all `1,961` pair children and not all
triples of the original core.

## 5. The two final triangle children

Both triangles in `(7)` are parent-heavy.  On the literal residual `R`
behind a triangle, let `h_R=|R|`, let `q_1` be its globally sealed allowed
singleton maximum, and put `gamma_R=(99/70)r_R/7`.  Since

```text
Delta_R=6h_R/7-q_1>0,
W_2=floor(gamma_R/Delta_R)+1,                              (8)
```

every pair with an endpoint at least `W_2` has union strictly below

```text
q_1+h_R/7+gamma_R/W_2<h_R.                               (9)
```

The exact finite pair head below `W_2`, combined with `(9)`, gives:

```text
triangle (15,24,33)
  h_R=57247/630630, r_R=32
  parent-heavy margin=661/19110
  q_1=1231/28028 at 27, W_2=191
  finite H_2=9469/145530 at {27,39}, 10 paid unions
  global no-two-cover margin=4399/112420308

triangle (24,30,33)
  h_R=33323/315315, r_R=34
  parent-heavy margin=2069/105105
  q_1=31235/756756 at 27, W_2=140
  finite H_2=293767/3783780 at {15,27}, 2 paid unions
  global no-two-cover margin=8023/33108075.                (10)
```

All inequalities in `(10)` are strict.  The two remaining triangle
children close, so the unresolved graph supports no hypothetical cover.
This proves the statement.

## 6. Root accounting and boundary

THM-2899 closes every other marked suffix on each of the `61` bodies, and
THM-2892 closes the eight-body chamber.  The `61` bodies are disjoint from
the fifteen-root THM-2895/2898/2899/2901 union.  They overlap THM-2902 in
exactly

```text
(1,2,3,4,5,10,12),       (1,2,3,4,5,12,13).              (11)
```

Thus THM-2903 adds `61-2=59` roots to the current seventeen-root proved
union.  The proved union has size

```text
17+59=76,
```

and the official seven-body residual becomes

```text
3432-76=3356 roots.                                      (12)
```

The result closes the one-hard bank, not all `12,919` H3-routed suffixes.
It does not close the `52` pair-cap exceptions, close the seven-body rung,
or prove LRC(14).

## 7. Verification

The verifier hash-pins the complete THM-2901 pair ledger and both interval
engines, reconstructs every literal carrier, scans actual rather than
cutoff-universe cores, checks the link identity on every retained vertex,
and independently subtracts every winning pair.  All theorem guards are
explicit and survive optimization.  Locked ordinary and optimized
replays agree.

The link-row digest and complete proof digest are

```text
b27ca638baf6194be3ee436916233389fe5b7ace7f49020bacc1bd1221c92ebb
08dc4a539544a417ff884ef4631b10d6eb14fd63d7b62b06e76d92e3d4d9b162.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.py
SHA-256 d62ff934f445b247a478d6b58f83f43f331003960291bce11f4ee313e6312707

05-knowledge/results/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.ledger.out
SHA-256 56236d7f8e51f12f80126c359415844762ff1b855b5d3b98d9b215781334754d

05-knowledge/results/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.out
SHA-256 5719083a83b275206907f47141fee8da2ba489194e31ba7c119f5313e3dfe73d
```
