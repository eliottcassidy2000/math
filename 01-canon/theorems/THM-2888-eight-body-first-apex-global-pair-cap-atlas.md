---
id: THM-2888
title: Eight-body first-apex global pair-cap atlas
status: >
  PROVED + FINITE-EXACT + VERIFIED.
  The 30030 literal first-apex carriers selected by THM-2885 have globally
  sealed exact two-comb caps.  Exactly 30025 caps are below 5h/7; 13802
  carriers close by a 2+2 pair partition, and the five nonpositive carriers
  close by an independently verified refined residual graph.  A genuine
  terminal-apex weighted gate, using no merely finite rank3 branch, closes
  1064 of 2508 active eight-body roots.  Together with THM-885's 495 low
  bodies this proves 1559/3003 whole roots.  The remaining 1444 roots and
  10202 honest nonterminal apex branches stay open; the j=5 rung and
  LRC(14) are not proved.
source: root/lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-2885-eight-body-global-top-fifteen-and-top-ten-hitting-gate
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-741-near-AP-four-slot-closure-all-2002-bodies-in-1-14
  - THM-885-covering-case-decomposition-j56-sweep
related:
  - THM-2883-ranked-apex-hitting-closure-all-thm741-roots
  - THM-2881-rank-impossible-pair-residual-closure-of-thirty-eight-thm741-roots
verification:
  - 04-computation/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.py
  - 05-knowledge/results/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.out
  - 04-computation/lrc14_j5_pair_hostiles_recursive_audit_codex_20260729.py
  - 05-knowledge/results/lrc14_j5_pair_hostiles_recursive_audit_codex_20260729.out
  - 04-computation/lrc14_j5_five_hostile_refined_k4_scout_codex_20260729.py
  - 05-knowledge/results/lrc14_j5_five_hostile_refined_k4_scout_codex_20260729.out
---

# THM-2888 — eight-body first-apex global pair-cap atlas

## 1. Universe and statement

Let

```text
E in C({1,...,14},8),        G_E = the good set of E,
m_E=|G_E|,                   c_E(w)=|G_E intersect D_w|.      (1)
```

THM-2885 gives ten global apex speeds `A_E`, with deterministic speed
tie-breaking, such that any five-speed external set whose scalar coverages
can exhaust `G_E` meets `A_E`.  For `a in A_E`, form the literal carrier

```text
C_(E,a)=G_E minus D_a,       h=|C_(E,a)|,
r=#components of C_(E,a),
c(w)=|C_(E,a) intersect D_w|,             w>=15, w!=a.       (2)
```

For distinct remaining speeds define the exact pair union

```text
U(x,y)=|C_(E,a) intersect (D_x union D_y)|.                  (3)
```

The exact atlas contains all

```text
3003 bodies x 10 apices = 30030 literal carriers.            (4)
```

For `30025` of them it proves a global rational cap

```text
U(x,y)<=B_(E,a)<5h/7                                         (5)
```

for every distinct `x,y>=15` outside the prefix.  Exactly five carriers
fail the last strict inequality; Section 5 closes all five by literal
residual computations.

Moreover,

```text
2B_(E,a)<h                                                   (6)
```

on `13802` carriers.  Any four remaining speeds can be partitioned into two
pairs, so `(6)` and `(3)` close all four slots immediately.

## 2. Global pair-tail seal

THM-735 and its discrepancy dependencies give the strict estimate

```text
c(w)<u(w):=h/7+(99/70)r/(7w).                               (7)
```

For every carrier the verifier scans all `15<=w<=2500`, excluding the fixed
apex.  Let `q_1` and `q_4` be the first and fourth coverages in this finite
scan.  It checks

```text
q_4>h/7,                u(2501)<=q_4,                       (8)
q_1<4h/7.                                                     (9)
```

Thus the scanned first four values are global.  The largest rank-one tail
threshold is

```text
263999736/239719 < 2501                                    (10)
```

at body `(4,6,8,9,10,11,12,14)`, apex rank `9`, speed `182`.

Put

```text
T_2=(99/70)r/[7(4h/7-q_1)].                                (11)
```

The largest value is

```text
1559025468/3415457 < 2501                                  (12)
```

at body `(1,2,6,8,9,10,13,14)`, apex rank `6`, speed `46`.
Consequently, if at least one endpoint of a pair is at least `2501`, then

```text
U(x,y)<=c(x)+c(y)<q_1+u(2501)<5h/7.                        (13)
```

Inside `15,...,2500`, order speeds by decreasing `c(w)`.  Branch-and-bound
uses `U(x,y)<=c(x)+c(y)` and pays an exact literal subtraction only while a
pair can improve the incumbent.  The finite-head maximum `H_2` and

```text
B_(E,a)=max(H_2,q_1+u(2501))                               (14)
```

therefore cap every pair globally.  There are `191306` paid pair
subtractions; the maximum on one carrier is `2037`.  Only `353` pairs in
the entire atlas even have scalar sum at least `5h/7`.

All comparisons in `(8)`--`(14)` preserve equality correctly.  In
particular `(7)` is strict, while `u(2501)<=q_4` is allowed; a cap equality
`B=5h/7` would be classified nonpositive, not silently promoted.

## 3. Scalar classes versus pair geometry

The globally sealed individual top four reproduce the independent scalar
census

```text
direct:          h-(p_1+p_2+p_3+p_4)>0,              13969;
rank3 finite:    direct fails and
                 h-(p_1+p_2+p_3)>h/7,                10939;
failure:         neither strict inequality,            5122. (15)
```

The unique equality in the second test is

```text
E=(1,2,4,7,8,10,12,14), apex rank 9, a=77,
h-(p_1+p_2+p_3)=h/7=211/6468.                         (16)
```

It is correctly in the failure class.  A finite rank3 head is not a
positivity proof and is never counted as a terminal apex below.

Among the `13802` pair-partition terminals in `(6)` are

```text
11654 scalar-direct, 2095 rank3-finite, and 53 scalar-failure carriers.
                                                                    (17)
```

Only the last two groups add genuinely new terminals beyond the scalar
direct class.

For every positive cap `(5)`, put

```text
delta=5h/7-B>0,
W=floor(2(99/70)r/(7delta)).                            (18)
```

Then

```text
B+2u(W+1)<h.                                           (19)
```

Hence an obstruction has at most one speed above `W`: choose any two such
tails as single combs and the other two speeds as the globally capped pair.
The smallest positive `delta` is

```text
709/8828820
```

at body `(2,4,6,7,8,10,12,14)`, apex rank `3`, speed `26`; its maximizing
pair is `(18,22)` and its raw `W` is `130827`.

This is a finite reduction, not a closure.  To normalize its scale, the
verifier treats the maximizing pair as an edge apex: it globally closes the
literal residual behind that edge, deletes the edge from the allowed pair
cap, and recomputes `(18)`.  This is done for the `29` positive scalar
failures with `W>2500` and the five nonpositive carriers below.  All `34`
forced-edge residuals close and every deleted-edge cutoff is at most `2044`.

## 4. The heavy-graph reframe

For later use, `(5)` has a useful exact graph interpretation.  Put

```text
theta=h-B>2h/7
```

and call `xy` heavy when `U(x,y)>=theta`.  If four speeds cover the carrier,
then for every pair `xy` and its complementary pair `zw`,

```text
h<=U(x,y)+U(z,w)<=U(x,y)+B,
```

so all six pairs are heavy: a cover induces a `K_4`.  Since
`U(x,y)<=c(x)+c(y)`, the finite set

```text
H={x:c(x)>=theta/2}
```

is a vertex cover of the heavy graph.  Every putative `K_4` therefore
contains a heavy triangle in `H`.  Its discrepancy cutoff is exactly the
same scale as `(18)`:

```text
(99/70)r/[7(theta/2-h/7)]
  =2(99/70)r/(7delta).                                  (20)
```

Thus the remaining proof obligation is a finite literal residual problem
indexed by heavy triangles, not an unrestricted four-speed enumeration.
No triangle closure is asserted in this theorem.

## 5. The five nonpositive pair caps

The five carriers and their unique pairs with `U>=5h/7` are

```text
E=(1,2,3,7,8,10,11,13), a=19 rank 5: (17,18);
E=(1,2,3,7,8,10,11,13), a=17 rank 7: (18,19);
E=(1,4,6,7,9,10,11,13), a=23 rank 6: (16,17);
E=(2,3,4,8,10,12,13,14), a=16 rank 6: (18,22);
E=(2,4,6,8,10,12,13,14), a=16 rank 4: (18,22).          (21)
```

The exact margins `5h/7-U` are respectively

```text
-50887/285170886,
-3650923/2851708860,
-480313/920551632,
-23/1070160,
-3553/1070160.                                         (22)
```

First, subtract each exceptional pair literally.  A new global pair cap on
the residual is strictly below its residual mass in all five cases; the
minimum margin is

```text
389293/63488880.                                        (23)
```

This closes every four-set containing the exceptional pair.

For a four-set avoiding it, the independent refined-hostile verifier deletes
that pair, constructs a nonexceptional cap `B<5h/7`, and applies the heavy
graph of Section 4.  Its five finite heads contain `1089` coarse heavy
edges.  For every one, the exact global pair cap on the literal two-edge
residual is below the residual mass.  Hence the refined graph has zero
edges, and in particular no `K_4`.  This closes every avoiding branch.
Together with `(23)`, all five carriers in `(21)` are genuine terminals.

The refined ledger digest is

```text
1d682c78bb1b2b9a9b907b3d25d62b15f86f85044334d0a19ea3b7839b6426f5.
                                                                    (24)
```

## 6. Genuine terminal-apex weighted gate

An apex is called terminal here **only** if

```text
its scalar direct margin is positive, or
2B<h, or
it is one of the five literally closed carriers in (21).            (25)
```

No rank3-finite branch and no change of distinguished observer is included.
There are

```text
16122 genuine terminal apices,
13908 raw nonterminal apex branches.                                (26)
```

For each root, filter its global THM-2885 top fifteen to retain every
nonterminal apex among ranks `1,...,10` and the five sentinel ranks
`11,...,15`.  If the five largest retained root coverages sum to less than
`m_E`, every scalar-dangerous five-set meets a terminal apex; a set avoiding
them closes by the union bound.  The five sentinels make the filtering
global: no later speed can enter the retained top five.

The exact genuine-terminal gate closes

```text
1297/3003 roots in the full atlas,
1064/2508 active roots meeting {13,14}.                              (27)
```

The five refined terminals do not change the root count, although their
individual closure is still exact.  THM-885 independently closes all
`C(12,8)=495` low bodies.  Therefore the disjoint whole-root composition is

```text
495+1064=1559/3003 closed eight-body roots.                           (28)
```

If any added speed is at most `14`, the extension acquires a ninth
in-window speed and global THM-741 closes it.  Thus `(28)` is a whole-root
statement, not merely a pure-tail census.

The residual for the next exact computation is

```text
1444 active roots,
10202 honest nonterminal apex branches inside those roots.           (29)
```

The root-gate ledger digest is

```text
9177cd9bcc34e301ae0297b011e533ba8bddd6a172be599e46efe8437580b180.
                                                                    (30)
```

## 7. Verification

The principal exact verifier and stored output are

```text
04-computation/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.py
SHA-256 cba433bce508ca8cc1e90c813e1988bb73c765ffafe350b74c3bad240eeca10f

05-knowledge/results/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.out
SHA-256 b6bcbf90257942523e7d26d33a2c06bde67805e1e008bfe190ccca7df0d83669
```

It hash-pins the THM-2885 and THM-2883 exact kernels and the independent
five-hostile closure.  The rational scalar tooth primitive and guarded
integer-vector primitive agree on `60060` values.  Every one of the `30030`
first paid pair carriers agrees between both subtraction orders and a
simultaneous subtraction; one carrier per body also agrees with direct
full-family reconstruction.

The canonical digests are

```text
full profile
  bd99500652a1ac5aff95932cc748e395799b8d11e2db71bef58b37e924959f35
hierarchical paid pairs
  ef3113e4aac5f352f9acdb7f9dde2d01890a1e6edbc11e2b0be6620c1e451f23
threshold pairs
  8c674340ca4be92e8baa34cd436956ef911d9fc3c1f8dc337fc4f4fdee6e589b
forced residuals
  4eeee139c7b0d9ba2d36957bfb1808997d623b526f682c96439cfd5b853a7cf9
heavy-edge normalization
  4615fd5b30d83006cdb10c8f1f2d64ce8a6a38d1e9c540358a4ba89c7c02c3b9.
                                                                    (31)
```

Ordinary and optimized replays with six workers are byte-for-byte identical.

## 8. Exact scope and anchored guardrail

This theorem proves the pair atlas, `1559` whole roots, the five hostile
carriers, and the reduction `(29)`.  It does not close the remaining heavy
triangles, the full eight-body/five-slot rung, or LRC(14).

In particular, recentering at `max(E)` changes the distinguished runner.
THM-741 in that new frame would prove loneliness of the new observer, not
the original anchored runner whose good set is `G_E`.  No such observer
drift is used in `(25)`--`(29)`. ∎
