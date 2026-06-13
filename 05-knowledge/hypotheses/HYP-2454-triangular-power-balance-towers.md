# HYP-2454 - Triangular power-balance towers stop at squares and expose a 78/90 support shadow

**Status:** OPEN synthesis; exact p=1/p=2 identities and finite power-center/bracket scout.
**Source:** codex-2026-06-12.
**Companions:** HYP-2457, HYP-2453, HYP-2452, HYP-2451, HYP-2450,
HYP-2445, HYP-2444, HYP-2430, HYP-2425, HYP-2128.
**Computation:** `04-computation/triangular_power_balance_towers_codex.py`;
stored output `05-knowledge/results/triangular_power_balance_towers_codex.out`.
**Beatty/Pell side-word addendum:** HYP-2456;
`04-computation/triangular_tower_beatty_pell_decomposition_codex.py`.
**Faulhaber anchor addendum:** HYP-2457;
`04-computation/triangular_faulhaber_anchor_expansion_codex.py`.
**External anchors:** OEIS A059270 (`https://oeis.org/A059270`) and
OEIS A059255 (`https://oeis.org/A059255`).

## Statement

Let `T_n=n(n+1)/2` and define the power-balance defect

```text
D_p(C,n) = sum_{i=0..n} (C-i)^p - sum_{i=1..n} (C+i)^p.
```

The user's two triangular towers are the first two exact positive integer
solutions of `D_p(C,n)=0`:

```text
p=1:  C = 2T_n = n(n+1)
p=2:  C = 4T_n = 2n(n+1)
```

Thus the first tower is the ordinary equal-sum interval balance

```text
n^2 + ... + (n^2+n) = (n^2+n+1) + ... + (n^2+2n),
```

and the second tower is the consecutive-square balance

```text
(2n^2+n)^2 + ... + (2n^2+2n)^2
  = (2n^2+2n+1)^2 + ... + (2n^2+3n)^2.
```

The working hypothesis is that the integer-center phenomenon genuinely stops
after powers `1` and `2`: for every `p>=3`, the positive root of `D_p(C,n)=0`
is trapped between the consecutive integers `2pT_n` and `2pT_n+1`.  The script
verifies this bracket for `3<=p<=8` and `n<=40`; for `p=3,4` the displayed
formulas make the first proof target quite concrete.

## Exact Algebra Already In Hand

The script records:

```text
D_1(C,n) = C - n^2 - n.
D_2(C,n) = C(C - 2n^2 - 2n).

D_3(C,n) = C^3 - 3n(n+1)C^2 - n^2(n+1)^2/2
          = C^3 - 6T_n C^2 - 2T_n^2.

D_4(C,n) = C(C^3 - 4n(n+1)C^2 - 2n^2(n+1)^2)
          = C(C^3 - 8T_n C^2 - 8T_n^2).
```

For `p=3`,

```text
D_3(6T_n,n)   = -2T_n^2 < 0,
D_3(6T_n+1,n) = 34T_n^2 + 12T_n + 1 > 0.
```

For `p=4`, after removing the zero center,

```text
D_4(8T_n,n)/(8T_n)       < 0,
D_4(8T_n+1,n)/(8T_n+1)  > 0.
```

The same bracket pattern appears computationally for `p=5..8`.  This suggests
a Bernoulli/Faulhaber proof route: write `D_p(2pT_n+k,n)` and prove sign at
`k=0,1`.

HYP-2457 sharpens this route.  With `u=n(n+1)` and `c=a+n`, the exact
midpoint equation keeps only odd Faulhaber moments:

```text
D_p(c,n)=c^p - 2*sum_{r odd} binom(p,r)c^(p-r)S_r(n).
```

The real root has formal expansion

```text
c_p(n)=p*u
  + (p-1)(p-2)/(12p)
  - (p-1)(p-2)(2p^2-4p-1)/(180p^3*u)
  + O(u^-2),
```

with the next `u^-2` coefficient also carrying `(p-1)(p-2)`.  Thus the p=1
and p=2 exact towers are explained by the same odd-moment/Bernoulli address
rather than by separate low-degree coincidences.

## First Tower: Square-Shell Partition

The ordinary tower is not only an equality.  Its row `n` partitions the
integer shell between consecutive squares:

```text
F_L(n) = [n^2, ..., n^2+n]
F_R(n) = [n^2+n+1, ..., n^2+2n]
F_L(n) union F_R(n) = [n^2, ..., (n+1)^2-1].
```

The common sum is

```text
S_1(n)=n(n+1)(2n+1)/2 = 3 * (1^2+...+n^2).
```

This is OEIS A059270 in the indexing used by the prompt.

## Second Tower: Square Balance And Ordinary Defect

The second tower has center `4T_n` and square common sum

```text
S_2(n)=n(n+1)(2n+1)(12n^2+12n+1)/6,
```

matching OEIS A059255.  If the squares are forgotten, the two ordinary sides
do not agree; instead they have a controlled defect:

```text
L_2(n)=n(n+1)(4n+3)/2,
R_2(n)=n(n+1)(4n+1)/2,

L_2(n)-R_2(n)=n(n+1)=2T_n,
L_2(n)+R_2(n)=2n(n+1)(2n+1)=4S_1(n).
```

This is the cleanest addition/multiplication bridge in the packet.  Squaring
turns the defective ordinary comparison into an equality, and the defect is
itself triangular.

## Crossover Atlas

The script verifies the user's visible overlaps:

```text
Q_L(2) = [10,11,12] inside F_L(3)=[9,10,11,12],
Q_R(2) = [13,14]    inside F_R(3)=[13,14,15],
Q_L(3) = [21,22,23,24] = F_R(4),
Q_R(3) = [25,26,27] inside F_L(5),
Q_L(4) = [36,37,38,39,40] inside F_L(6).
```

Through `Q` rows `n<=100` and `F` rows `m<=150`, the only exact full-side
equality is

```text
Q_L(3) = F_R(4) = [21,22,23,24].
```

The length constraint already explains why this is exceptional: side equality
forces `m=n+1`; the endpoint equations then give `n=3` in the nontrivial
case.

Endpoint alignments continue in Pell-like families:

```text
2n^2+n       = m^2,
2n^2+2n     = m^2+m,
2n^2+2n+1   = m^2,
2n^2+3n     = m^2+m,
```

and their start/end variants.  The stored run finds `25` such boundary
coincidences through `Q` row `100`, including the conspicuous values
`36, 90, 144, 210, 420, 421, 840, 841, 4900, 14280, 14281`.

## The 78/90 Support Shadow

The row `Q(3)` is the best transfer beacon:

```text
Q_L(3) = [21,22,23,24], ordinary sum 90,
Q_R(3) = [25,26,27],    ordinary sum 78,
Q_L(3)^2 sum = Q_R(3)^2 sum = 2030.
```

The scout finds `L_2(3)=90=S_1(4)` as the unique checked `L_2(n)=S_1(m)` hit
for `n,m<=500`, while `R_2(3)=78=C(13,2)`.

This matters because `78` is already present in the repo as the `lambda_5`
of the hypothetical Type II minimum-word design `5-(72,16,78)`, and also as
the `D7` index in HYP-2445's product-quotient support-gate atlas.  HYP-2454
does not claim this constructs a self-dual `[72,36,16]` code.  It says the
same `78` is now linked to an adjacent additive shadow `90=S_1(4)` by an
exact square-balance equality, so the next support ledger should track the
pair `(78,90)` rather than either scalar alone.

## Transfer To LRC14 And Hidden Lifts

HYP-2128 already identifies the LRC worry modulus by the triangular identity

```text
8*C(n,2)+1 = (2n-1)^2.
```

HYP-2454 adds a sibling identity family: the row centers `2T_n` and `4T_n`
are the only observed integer centers where a left-heavy interval can balance
the immediately following right interval after taking powers.

For LRC14, this suggests a defect ledger:

```text
ordinary side defect     -> triangular resource 2T_n
square equality          -> multiplicative closure of that defect
endpoint Pell alignment  -> shell/clock boundary address
78/90 row                -> code-design/LRC-resource beacon
```

For `[72,36,16]`, the bridge should be tested at the support level.  Weight
enumerator coefficients are scalar totals; minimum-word supports and their
incidences are the hidden lift.  The `78/90` pair is therefore a candidate
address in the support-incidence lift, parallel to HYP-2452's convolution
lift of polynomial coefficients.

## Tournament Analysis

Tournament vertices were deliberately not the interval entries alone.
Alternatives considered: rows, intervals, endpoints, centers, sums, defects,
square-equality shadows, Pell boundary events, code design parameters, LRC
shell moduli, Fourier modes, and proof obligations.

The stored tournament uses vertices:

```text
power_center_rigidity,
78_90_code_shadow,
square_tower_multiplier_lift,
unsquared_defect_channel,
first_square_row_partition,
pell_endpoint_alignments,
overlap_atlas,
lrc14_worry_modulus_bridge,
convolution_boundary_lift.
```

The pairwise observable is the majority comparison of

```text
(exactness, novelty, transfer_to_lrc, transfer_to_72, computability, proof_potential),
```

with the listed order as the fixed Hamiltonian tie path.  The run gives a
nontransitive proof-carrier tournament:

```text
ranking leader: 78_90_code_shadow
score_hist = {0:1, 2:2, 3:1, 4:1, 5:1, 6:1, 7:2}
directed_3cycles = 6
scc_sizes = [8,1]
hamiltonian_paths = 53
edge_flips_vs_exactness = 12
```

This is a useful warning.  The exact identities are not automatically the
highest-leverage proof route; the best next work may be the support transfer
or hidden-lift interpretation.

## Open Tasks

1. Prove or refute the general bracket
   `D_p(2pT_n,n)<0<D_p(2pT_n+1,n)` for all `p>=3`.
2. Use HYP-2456's Beatty/Pell normal form as the side-word classifier, then
   extend the same endpoint-wall treatment to the remaining `F`/`Q` power
   balance families.
3. Turn the `Q(3)` `78/90` shadow into a concrete `[72,36,16]`
   minimum-design support-ledger constraint.
4. Attach the same defect ledger to LRC14: compare the `27` shell, `78`,
   `90`, and `91` resources.
5. Lift interval balances into convolution/support balances, following
   HYP-2452's hidden-boundary-total program.
