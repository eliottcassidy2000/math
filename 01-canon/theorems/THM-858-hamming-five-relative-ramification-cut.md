---
id: THM-858
title: Hamming-five relative-ramification cut and uniform effective-order bound
status: PROVED STRUCTURAL + FINITE-EXACT — uniform complement-lcm cut; no new common-sheet language through effective order 21; every common-sheet effective order is at most 42,336
source: codex-2026-07-15-S10 Hamming-five relative-order continuation
depends_on: [THM-810, THM-823, THM-844, THM-845, THM-847]
related: [THM-815, THM-840, THM-857, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_effective_order_21_common_sheet_census_codex_S10.cpp
  - 05-knowledge/results/lrc13_hamming_five_effective_order_21_common_sheet_census_codex_S10.out
---

# THM-858 — Hamming-five relative ramification

Put

```text
delta=1/13,                         [12]={1,...,12}.
```

Use THM-823's oriented five-colour sheet model.  Thus the replacement labels
are a five-set `R={r_1,...,r_5} subset F_13^*`; colour `i` has effective order
`D_i>=1`, with `13` not dividing `D_i`, and unit `e_i mod D_i`.  Its CRT lift
`u_i mod 13D_i` is defined by

```text
u_i=D_i r_i (mod 13),                 u_i=e_i (mod D_i). (1)
```

For `c=lcm(D_1,...,D_5)`, the sheets covered by colour `i` at owner `o in R`
are

```text
E_i(o)={ell in Z/cZ:
        -D_i < <u_i(o^(-1)+13ell)>_(13D_i) <= D_i}.     (2)
```

A **common-sheet presentation** satisfies

```text
union_i E_i(o)=Z/cZ                         for every o in R. (3)
```

The half-open orientation in (2) is the one fixed in THM-810/823.

## Theorem

### A. Uniform relative-ramification cut

Let `S` be any nonempty set of colours, put `T={1,...,5}\S`, and use the
convention that the lcm of the empty set is one.  Define

```text
L_T=lcm(D_j:j in T),
m_i=D_i/gcd(D_i,L_T),                         i in S,
rho(m)=ceil(2m/13)/m.                                      (4)
```

Every common-sheet presentation satisfies the simultaneous family of cuts

```text
sum_(i in S) rho(m_i) >= 1             for every nonempty S. (5)
```

This is an all-order statement.  It is independent of the chosen unit word
and is not a scalar cutoff on the individual `D_i`.

### B. Prime-power carrier rigidity

The cuts (5) have the following uniform consequences.

1. Every effective order is lcm-redundant:

   ```text
   D_i divides lcm(D_j:j!=i)                         for every i. (6)
   ```

   Equivalently, no maximal prime power is carried by only one colour.

2. Fix a prime `p!=13`, and let `S_p` be the colours on which `v_p(D_i)` is
   maximal and positive.  If `|S_p|=2`, then necessarily

   ```text
   p=2,                  m_i=2 for both i in S_p.       (7)
   ```

   If `p` is odd and `|S_p|=3`, then necessarily

   ```text
   p=3,                  m_i=3 for all i in S_p.        (8)
   ```

   In particular, a maximal seven-power is carried by at least four colours.

3. No prime `p>=11`, `p!=13`, divides a surviving effective order.  If five
   divides one of the orders, then all five colours carry the maximal
   five-adic exponent and (5) forces

   ```text
   D_i in {5,10,15,20}                         for every i. (9)
   ```

Part C and THM-823 reject (9).  Consequently every still-unclassified
common-sheet presentation has

```text
D_i is {2,3,7}-smooth for every i.                       (10)
```

Thus the live arithmetic object is not five unrelated orders: it is the
hypergraph of shared maximal prime powers, decorated by the complement-lcm
fibres appearing in (4).

4. The relative cuts at every adjacent `p`-adic valuation level give the
   uniform range bounds

   ```text
   max_i v_2(D_i)-min_i v_2(D_i) <= 5,
   max_i v_3(D_i)-min_i v_3(D_i) <= 2,
   max_i v_7(D_i)-min_i v_7(D_i) <= 1.                (10a)
   ```

   THM-823 supplies an index `i_*` with `D_(i_*)<=21`.  Combining (10) and
   (10a), every effective order in every common-sheet presentation satisfies

   ```text
   D_i <= D_(i_*)*2^5*3^2*7 <= 21*2016 = 42,336.      (10b)
   ```

This makes the remaining effective-order presentation bank uniformly finite.
It does not by itself decide the arbitrary metric lifts within each surviving
presentation language.

### C. The exact effective-order-at-most-21 boundary

Consider the no-order-one bank

```text
1 in R subset F_13^*,       |R|=5,
D_i in {2,...,12,14,...,21}.                              (11)
```

The normalization means that label `1` is present; rows are not quotiented by
the remaining multiplicative action.  Separate as the **new branch** the rows
with at least one `D_i>12`.  The exact census is

```text
conceptual normalized presentations                    817,112,670
scalar-capacity presentations                                4,245
scalar presentations with every D_i<=12                     2,190
scalar presentations in the new branch                       2,055
rejected by singleton cuts (6)                                1,870
lcm-redundant residual rows                                     185. (12)
```

The `185` residual rows have exactly eight unordered order patterns:

| order pattern | presentations |
|---|---:|
| `(2,2,5,20,20)` | 50 |
| `(2,2,9,9,18)` | 40 |
| `(2,2,10,20,20)` | 30 |
| `(2,3,3,10,15)` | 20 |
| `(2,3,3,18,18)` | 10 |
| `(3,3,3,5,15)` | 20 |
| `(3,3,3,18,18)` | 10 |
| `(3,10,10,15,15)` | 5 |

Top-prime instances of (5) reject all `185`, partitioned by the first
rejecting prime as

```text
p=2: 20,                 p=3: 40,                 p=5: 125. (13)
```

There are therefore zero common-sheet presentations in the new branch of
(11).  Independently, a literal replay of every unit word on all `185`
residual rows checks

```text
unit words                                                51,360
common-sheet unit words                                        0
direct/affine membership identities                     264,960
owner-capacity identities                                 18,288
best minimum sheet coverage                                  4/5. (14)
```

This source does **not** scan order one.  THM-823(E) already proves, without
an order bound, that a row containing `D_i=1` is either all order one or one
order-one colour plus an order-three THM-810 quartet.  THM-823(C) classifies
the no-order-one bank through twelve.  Combining those symbolic branches with
(11)--(14), every common-sheet presentation with all `D_i<=21` belongs to one
of the already known all-one, all-three, or mixed one-plus-three languages.
THM-845, THM-844, and THM-847 respectively close their proper metric lifts.
Thus no new sporadic common-sheet language occurs through effective order
`21`.

Together with THM-823's scalar pivot `min_i D_i<=21` and (10b), every
unclassified no-order-one common-sheet row now lies in the precise finite
strip

```text
2<=min_i D_i<=21 < max_i D_i<=42,336,
all D_i {2,3,7}-smooth,
D_i divides lcm(D_j:j!=i),
and every relative cut (5) holds.                         (15)
```

### D. Tournament audit and the faithful carrier

On each of the `185` lcm-redundant rows, take the five colours as vertices.
The raw pair observable is

```text
rho(D_i)-rho(D_j).                                       (16)
```

For the complement-conditioned gauge on the pair `{i,j}`, first quotient
`D_i,D_j` by their gcd with the lcm of the other three orders and then compare
the two resulting `rho` values.  The switch is raw to complement-conditioned;
ties are oriented by increasing replacement label.  Both gauges give a
transitive tournament on every row:

```text
score histogram                 {0:1,1:1,2:1,3:1,4:1}
directed triangles                                      0
SCC-size histogram                                  {1:5}
Hamiltonian paths                                        1
raw/conditioned edge flips                             574
flip histogram                    {0:8,1:17,2:38,3:49,
                                   4:38,5:28,6:7}
raw ties / conditioned ties                      610 / 1,850. (17)
```

This is a negative tournament result with a constructive replacement.  Pair
rankings remain transitive even while the subset cuts reject every row.  They
forget which complement misses which owner sheet, which set carries a top
prime power, and how the affine intervals meet the complement fibre.

The predicate-preserving object has three layers:

```text
prime-power carrier hyperedges
    -> complement-lcm fibres and relative orders m_i
    -> owner-labelled affine cyclic intervals.           (18)
```

The first layer tells which cuts can fire; the second gives their fractional
capacity; the third retains literal common-sheet phase.  Bare runners,
replacement colours, pairwise order comparisons, or a gain-free tournament
destroy at least one of these data.  In this branch, hyperedges and fibres,
not runners or arcs, are the useful vertices and incidence objects.

## Proof

### 1. Every complement misses an owner sheet

For a colour of order `D`, replacement label `r`, and owner `o`, write

```text
C_D(r,o)=|E_i(o)|/c.
```

This is THM-823's oriented scalar capacity.  We first show that every set
`T` of at most four colours has an owner `o in R` such that

```text
sum_(j in T) C_(D_j)(r_j,o)<1.                          (19)
```

Augment `T`, if necessary, to four colours `T'`.  If those four colours fail
scalar coverage at one of their own four owners, (19) follows immediately.
Otherwise THM-810's equality classification applies: either all four orders
are one, or all four orders are three and their labels form a multiplicative
coset of `{1,5,8,12}`.  At the fifth owner the order-one alternative
contributes zero.  In the order-three alternative, THM-823's oriented coset
matrix gives total contribution zero or `2/3`.  It is again strictly below
one.  Removing the colours added to form `T'` preserves the strict deficit,
proving (19).

For `j in T`, sheet membership is periodic modulo `D_j`.  Since `D_j` divides
`L_T`, (19) implies

```text
|union_(j in T) E_j(o) inside Z/L_T Z|<L_T.
```

Choose a residue `ell_0 mod L_T` missed by every colour in `T`.

### 2. Affine intervals on the missed fibre

Fix `i in S` and put `D=D_i`, `e=e_i`.  Represent
`a=r_i o^(-1) mod 13` in `{1,...,12}`.  Every residue in (2) can be written

```text
D a+13k,
k=e ell+e o^(-1)13^(-1)                    (mod D).     (20)
```

The condition in (2) becomes

```text
-D < D a+13k <= D.                                     (21)
```

Thus the eligible `k mod D` form one consecutive cyclic interval, obtained
from a real half-open interval of length exactly `2D/13`.

Now restrict `ell` to the full common-sheet fibre

```text
F={ell in Z/cZ:ell=ell_0 (mod L_T)}.                    (22)
```

As (22) is traversed, (20) advances by the unit multiple `eL_T mod D`.
It therefore samples one gcd coset of size

```text
m=D/gcd(D,L_T),                                        (23)
```

with every coset point repeated equally often if the global fibre is larger.
After dividing the lifted interval (21) by `gcd(D,L_T)`, its length is
`2m/13`.  A half-open interval of this length contains at most
`ceil(2m/13)` integer points.  Hence colour `i` covers at most the fraction

```text
rho(m)=ceil(2m/13)/m                                   (24)
```

of `F`.  Every colour in `T` misses all of `F`.  If (3) holds, the colours in
`S` must cover `F`; the union bound and (24) give (5).

### 3. Arithmetic consequences

For a singleton `S={i}`, (5) gives `rho(m_i)>=1`.  For every integer `m>=2`,
`ceil(2m/13)<m`, so `m_i=1`.  This is exactly (6).

Let `S_p` be a top-prime carrier set.  The complement lcm has smaller
`p`-adic valuation, so every relative order `m_i`, `i in S_p`, is divisible
by `p`.  The elementary bounds

```text
m>=3                 => rho(m)<=1/3, equality only at m=3;
m even               => rho(m)<=1/2, equality only at m=2;
5 divides m          => rho(m)<=1/5,
                         equality only for m in {5,10,15,20};
p divides m, p>=11,
p!=13                => rho(m)<1/5                         (25)
```

now give (7)--(9).  For the last inequality, `m>=22` gives
`rho(m)<2/13+1/22<1/5`; the only smaller prime-divisible cases are
`m=11,17,19`, checked directly.  Two odd carriers contribute at most `2/3`.
Three odd carriers reach one only at `m_i=3`.  Four five-carriers contribute
at most `4/5`; five carriers force equality in every term, yielding (9).
Five carriers of a prime at least eleven still contribute strictly less than
one.  Part C and THM-823 eliminate (9), leaving (10).

It remains to prove the uniform cap (10a).  Fix `p in {2,3,7}` and list the
distinct values of `v_p(D_i)` in increasing order.  At any adjacent-level
boundary, let `S` be the colours strictly above the lower level and `T` the
colours at or below it.  Then `T` is nonempty, so `|S|<=4`, and every relative
order in the cut (5) is divisible by `p^g`, where `g` is the gap between those
two adjacent valuation levels.

The required exact bounds are

```text
2|m  => rho(m)<=1/2,       4|m  => rho(m)<=1/4,
16|m => rho(m)<=3/16,
3|m  => rho(m)<=1/3,       9|m  => rho(m)<=2/9,
7|m  => rho(m)<=2/7,      49|m  => rho(m)<=8/49.       (25a)
```

For example, if `m=4k`, then `ceil(8k/13)<=k`, and if `m=16k`, then
`ceil(32k/13)<=3k`; the other lines are identical one-step ceiling bounds.

For `p=2`, a boundary can have only `|S|=2,3,4`.  At sizes two and three,
`g>=2` would make every `m_i` divisible by four and the sum in (5) at most
`1/2` or `3/4`; hence `g<=1`.  At size four, `g>=4` would make every `m_i`
divisible by sixteen and the sum at most `3/4`; hence `g<=3`.  As the boundary
moves upward, `|S|` strictly decreases, so the total valuation range is at
most

```text
3+1+1=5.                                                (25b)
```

For `p=3`, sizes one and two are impossible.  At sizes three or four, `g>=2`
makes the cut sum at most `2/3` or `8/9`; hence both possible gaps are at most
one and the total range is at most two.  For `p=7`, sizes at most three are
impossible, while at size four a gap at least two gives sum at most `32/49`.
Thus there is at most one gap and it is one.  This proves (10a).

Choose `i_*` using THM-823's `min_i D_i<=21`.  Since all prime divisors are
now in `{2,3,7}`, the three range bounds imply, prime by prime,

```text
v_2(D_i)<=v_2(D_(i_*))+5,
v_3(D_i)<=v_3(D_(i_*))+2,
v_7(D_i)<=v_7(D_(i_*))+1.
```

Multiplying the three inequalities gives (10b).

### 4. Exact census and replay

Multiplication in `F_13^*` preserves all owner ratios, so any label set can be
normalized to contain `1`.  There are `C(11,4)=330` such label sets and
nineteen allowed no-order-one values in (11), giving

```text
330*19^5=817,112,670                                  (26)
```

conceptual rows.  The source evaluates scalar capacities as integers after
multiplication by

```text
232792560=lcm(1,...,22),                               (27)
```

and uses only exact optimistic upper bounds to prune the order tree.  It then
applies the singleton cuts, groups the survivors into the eight displayed
patterns, and applies the top-prime cuts.  Every count in (12)--(14) is
asserted in the source.

The literal cross-check does not reuse the top-prime verdict.  For every unit
word on every lcm-redundant row it constructs the centred CRT masks from (2)
and asks whether their union is every sheet at every owner.  A separate affine
implementation of (20)--(21) agrees on all `264,960` membership decisions and
on all `18,288` owner capacities.  Optimized and unoptimized builds produce
byte-identical output, an AddressSanitizer/UndefinedBehaviorSanitizer build is
clean, and the stored hashes are

```text
source SHA-256  7d77552c0873a59532a509fe08f88b924c625ce5a308c7f5d00a0166c7dc98f1
output SHA-256  19ec318f781f883d2e1b25119fd3a726616559b894dc196317d19ac745172134
payload FNV-64  fb30afb76cb41691.                            (28)
```

The same exact row loop produces the tournament fingerprints (17).  Since
all `185` rows are rejected only after a subset/complement operation while
both pair gauges remain transitive, (18) is forced as the faithful carrier
for this proof rather than inferred from a visual analogy. ∎

## Scope guardrail

This theorem closes the common-sheet bank with `max_i D_i<=21` and proves that
the remaining effective-order bank is finite; it does not
prove that common-sheet presentations themselves are absent, since the known
all-one/all-three/mixed equality languages remain before their metric lifts
are tested.  It does not classify the finite strip (15), close the arbitrary
metric lifts of any new presentation there, handle other arbitrary-scale
Hamming-five deck ramifications, close arbitrary-scale/common-sheet Hamming six, or
settle the global `n=12` sporadic branch.
