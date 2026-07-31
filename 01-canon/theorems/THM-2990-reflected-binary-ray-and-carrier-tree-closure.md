---
id: THM-2990
title: "Reflected binary ray and carrier-tree closure"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.  The candidate proves a
  sufficient tree-Hunter cone for every nonconstant six-label reflected k=1
  packet, closes the single binary ray with body (1,2,4,9,11,12) and levels
  (2Q,Q,Q,Q,Q,2Q) for every Q>=1 by one addressed endpoint word, and gives
  that ray an all-Q carrier-averaged spanning-tree certificate.  It makes no
  arbitrary-k=1, denominator-ledger, first-drift-cap, or LRC(14) closure claim.
source: codex-reflected-ray-carrier-audit-2026-07-30
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCTreeHunter.lean
verification:
  - 04-computation/lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.out
  - 04-computation/lrc14_j7_reflected_binary_dilation_resonance_referee_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_binary_dilation_resonance_referee_thm2941.out
  - 04-computation/lrc14_j7_reflected_hostile_carrier_tree_allq_bridge_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_hostile_carrier_tree_allq_bridge_thm2941.out
---

# THM-2990 -- reflected binary ray and carrier-tree closure

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**

This theorem candidate has three deliberately separate conclusions.

1. A translation-free pair floor, summed on a spanning tree, gives a uniform
   sufficient cone for every nonconstant reflected level word.
2. One concrete ray outside that cone is closed for every height by retaining
   a body-safe cell and its endpoint-owner word.
3. Averaging actual pair overlaps over the whole body carrier gives the same
   ray a spanning-tree certificate for every height.

The third conclusion is a structural re-expression of the second, not a new
row exclusion.  Failure of the sufficient cone is only an address and never a
physical survivor statement.  Arbitrary reflected `k=1`, the zero-aligned
sector, the rung, and LRC(14) remain open.

## 1. Reflected cells and the exact singleton debt

Let `E=(e_0,...,e_5)` be the increasing enumeration of a six-element subset
of `{1,...,14}`, and put

```text
L=14 lcm(E).
```

For positive integral levels `q_i`, define the six reflected labels

```text
z_i=q_i L-e_i.                                           (1)
```

Let `J_E` be the integer body-safe cells from THM-2941.  On a cell `j in J_E`
use the normalized coordinate `t=(j+u)/L`, `0<=u<=1`, and write

```text
A_i(j)={u: ||z_i(j+u)/L||<1/14},
U_j=union_i A_i(j).
```

The reflected-cell identity gives, exactly in Lebesgue mass,

```text
mu(A_i(j))=q_i L/[7(q_i L-e_i)]
            =1/7+e_i/[7(q_i L-e_i)].                   (2)
```

Thus the six singleton masses have the same value on every body-safe cell,

```text
sum_i mu(A_i(j))=6/7+epsilon(E,q),
epsilon(E,q)=sum_i e_i/[7(q_iL-e_i)].                  (3)
```

By THM-2941, one cell with `mu(U_j)<6/7` makes the projected residual have
mass greater than `1/7`, so it cannot fit in the one remaining aligned danger
comb.  All closures below use this strict implication.

## 2. A cellwise pair floor

Fix a positive integer base `b` and two indices `i,k` with unequal levels.
Put

```text
lambda_i=q_i-b-e_i/L,       lambda_k=q_k-b-e_k/L,
S=|lambda_i|+|lambda_k|,    v=lambda_i-lambda_k.        (4)
```

For `ell>0`, define

```text
G(ell)= [ floor(ell)/49
          + max(0,{ell}-5/7)^2/4 ] / ell.              (5)
```

Then on every body-safe cell,

```text
mu(A_i(j) intersect A_k(j))
 >= G(|v|)-(4S+|v|/2)/b.                               (6)
```

Here is the mechanism.  Subdivide `u=(r+x)/b` and put `s_r=r/b`.  Modulo an
integer, the exact phase is

```text
z_i(j+u)/L
 =x+lambda_i(j+s_r+x/b).                               (7)
```

Deleting the last term gives a moving limiting arc.  The intersection of the
two limiting arcs has mass

```text
phi(v(j+s)),        phi(y)=max(0,1/7-||y||).            (8)
```

The least integral of this period-one tent over an interval of length `ell`
is

```text
floor(ell)/49+max(0,{ell}-5/7)^2/4.                    (9)
```

Indeed, each complete period contributes `1/49`; a remainder of length at
most `5/7` fits in the zero gap, and any excess is split equally between the
two linear tent flanks.  Dividing by `ell` proves the `G` term.

The two exact arcs differ from their sampled limits only within the two
endpoint neighborhoods, costing at most `4S/b`.  Since `phi` is one-
Lipschitz, left Riemann sampling costs at most `|v|/(2b)`.  This proves `(6)`.
The exact referee checks the direct multiplier pullback against `(7)`, the
tent primitive, and the error directions with rational arithmetic.

## 3. Summing pair floors on a tree

For any spanning tree `T` on the six indices, Hunter's inequality gives

```text
mu(U_j)+sum_(ik in T) mu(A_i(j) intersect A_k(j))
 <=6/7+epsilon(E,q).                                  (10)
```

Let

```text
m=min_i q_i,              D=max_i q_i-m.               (11)
```

Assume `D>=1`.  Join the minimum-level vertices to the higher-level vertices
by a bipartite spanning tree.  For an edge whose higher endpoint has level
gap `d>=1`, orient it from its minimum endpoint `i` to its higher endpoint
`k` and put

```text
delta=(e_i-e_k)/L.
```

The exact all-body bounds are

```text
-11/168<=delta<=11/168.                                (12)
```

With base `b=m`, the two offsets in `(4)` have opposite signs, so

```text
S=|v|=d+delta.                                         (13)
```

Also `(3)` and `sum(E)/L<=1/6` give

```text
epsilon(E,q)<=1/(39m).                                 (14)
```

Consequently a sufficient tree condition is

```text
m sum_(ik in T) G(d_ik+delta_ik)
 > (9/2)(sum_(ik in T)d_ik + sum_(ik in T)delta_ik)
   +1/39.                                              (15)
```

The inequality is strict.  Equality in a later rounded threshold is not
silently accepted.

### 3.1 Exact all-body tree census

For each of the `3003` bodies and each of the `62` nontrivial choices of the
minimum-label part, give an eligible bipartite edge the unit-gap floor

```text
G(1+(e_i-e_k)/L).                                      (16)
```

Kruskal's algorithm chooses a maximum-floor tree and, among ties, minimizes
the signed label correction.  The `186186` exact rows prove that if `H` is
the chosen floor sum and `C` its correction, then

```text
(45/2)/H < 313,
[(45/2)+(9/2)C+1/39]/H < 297.                          (17)
```

The unique simultaneous extremal body/part is

```text
E=(1,2,3,4,6,12),       minimum-label mask=31,
H=142918133/1987853616, C=-11/42.                      (18)
```

For fixed `delta`, `G(d+delta)` is nondecreasing in the integer gap `d>=1`.
Since a tree has five edges and every gap is at most `D`, `(17)` makes `(15)`
strict whenever

```text
m>=313D-16.                                            (19)
```

In particular the `D=1` branch closes for `m>=297`.

### 3.2 Gap-sensitive sharpening for `D>=2`

The exact worst label-corrected edge floor is

```text
g_*(d)=min_(|delta|<=11/168) G(d+delta)
       =G(d-11/168)
       =(2304d-935)/[672(168d-11)].                    (20)
```

On `[-11/168,0]`, the derivative numerator is

```text
d(98delta+24)+49delta^2>0,                             (21)
```

so the negative branch is minimized at the left endpoint.  On the positive
branch the minimum is at `11/168`, and direct subtraction from the value in
`(20)` gives

```text
11(9672d+935)/[672(168d-11)(168d+11)]>0.               (22)
```

Thus `(20)` is uniform.  Moreover both `g_*(d)` and `d/g_*(d)` increase with
`d`; their forward-difference numerators are respectively `5489` and

```text
387072d^2+72912d-146795>0       (d>=1).                (23)
```

Choose the label-optimized minimum/higher double-star.  Its signed label
correction satisfies `C<=11/42`.  Applying `(20)` edge by edge in `(15)` gives
the sufficient threshold

```text
R(D)=3024D(168D-11)/(2304D-935)+330328/17797.           (24)
```

The forward difference of `R(D)-221D` is the negative of

```text
[2654208D^2+499968D+161024765]
 /[(2304D-935)(2304D+1369)],                           (25)
```

and `R(2)<554`.  Hence every nonconstant reflected packet closes under the
clean sufficient cone

```text
D=1:   m>=297;
D>=2:  m>=221D+112.                                   (26)
```

Equivalently, failure of this certificate forces

```text
D=1:   m<=296;
D>=2:  m<=221D+111.                                   (27)
```

Statement `(27)` is an **uncertified wedge**, not a survivor classification.
The constant case `D=0` is the common-level diagonal already closed by
THM-2941.

## 4. The binary ray outside the cone

Now fix

```text
E=(1,2,4,9,11,12),          L=5544,
(q_0,...,q_5)=(2Q,Q,Q,Q,Q,2Q),       Q>=1.             (28)
```

This has `m=D=Q`, so `(26)` never reaches it.  It is therefore a genuine test
of whether an address can replace a bounded-spread hypothesis.

On a body-safe cell, write `u=(r+x)/Q`, `K=QL`, and `J=Qj+r`.  Label by label,

```text
(a_eQL-e)(j+(r+x)/Q)/L
 =a_ex-e(J+x)/K                 (mod 1),                (29)
```

where `a=(2,1,1,1,1,2)`.  Thus the coarse union mass is the average of `Q`
binary-kernel masses at the consecutive fine cells `Qj,...,Qj+Q-1`.

The body-safe cell `j=L/14=396` is a sharp hostile.  Its `8Q` arcs are
pairwise disjoint in the repeated owner order

```text
1_even,2,4,12_even,1_odd,9,11,12_odd,
```

and hence

```text
mu(U_396)=sum_e a_eQL/[7(a_eQL-e)]
          =6/7+sum_e e/[7(a_eQL-e)]>6/7.               (30)
```

So neither the first safe cell nor an unaddressed minimum is monotone under
heterogeneous dilation.

At the different body-safe cell

```text
j_*=2L/7=1584,                                         (31)
```

the low labels phase-lock in pairs: `(2,9)` at `4/7` and `(4,11)` at `1/7`.
The exact endpoint-owner word has, for every `0<=r<Q`, three components:

```text
({4,11,1,12}: left owner 4, right owner 12),
({2,9,1}:     left owner 2, right owner 1),
({12}:        left and right owner 12).                (32)
```

All consecutive endpoint differences are positive affine forms in
`Q>=1, 0<=r<Q`; cross multiplication supplies `27` exact sign certificates.
Summing the component lengths gives

```text
M_Q=mu(U_(j_*))
 =33Q(511253875440Q^3-994881888Q^2+494683Q-73)
  /[(924Q-1)(1386Q-1)(2772Q-1)(11088Q-1)].             (33)
```

Furthermore

```text
1/2-M_Q=P(Q)/[2(924Q-1)(1386Q-1)(2772Q-1)(11088Q-1)],

P(Q)=5619650962464Q^4-23087810592Q^3
     +31384122Q^2-11352Q+1.                            (34)
```

For `Q>=1`, group `(34)` as

```text
Q^3(5619650962464Q-23087810592)
 +Q(31384122Q-11352)+1>0.                              (35)
```

Thus

```text
M_Q<1/2<6/7                                             (36)
```

for every `Q>=1`.  By Section 1, the entire ray `(28)` is empty.  Its limiting
mass is `9505/22176=3/7+1/22176`.

This is one labelled ray.  Permuting labels together with their level types
is harmless, but changing the body, the type assignment, or the reflected
residue packet is not covered.

## 5. Carrier-averaged tree certificate on the same ray

For `(28)`, the body carrier has exactly `|J_E|=2260` safe unit cells.  Define
the actual carrier weights

```text
omega_ik(Q)=|J_E|^(-1) sum_(j in J_E)
              mu(A_i(j) intersect A_k(j)).             (37)
```

Averaging `(10)` shows that any spanning tree `T` with

```text
sum_(ik in T) omega_ik(Q)>epsilon(E,q)                  (38)
```

forces the average of `mu(U_j)` below `6/7`, and hence supplies a closing
cell.  This is the exact source-to-target map:

```text
cellwise overlap -> carrier average -> one low-union cell.
```

It forgets which cell closes.  The endpoint word `(32)` is the sidecar that
restores enough local mass for the infinite tail.

The unequal-level graph of `(28)` is `K_(2,4)`.  It is connected and has a
spanning tree, but it has no unequal-level perfect matching.  This explains
why the earlier matching topology stalls while a tree can still use physical
same-level overlaps.  It does **not** say that every translation-free pair
floor is zero.

At `j_*`, the fixed tree

```text
T_*={1-9,1-11,2-9,4-11,11-12}                          (39)
```

induces a tree on the active label set of every open endpoint segment in
`(32)`.  Its five overlaps therefore telescope exactly:

```text
sum_(ik in T_*) mu(A_i(j_*) intersect A_k(j_*))
 =6/7+epsilon(E,q)-M_Q>5/14.                           (40)
```

Every other cell contribution is nonnegative, so

```text
sum_(ik in T_*) omega_ik(Q)>5/(14*2260).                (41)
```

Also `epsilon_Q<=epsilon_1/Q`, and the exact value satisfies

```text
5 < (14*2260/5)epsilon_1 < 6.                          (42)
```

Equations `(41)`--`(42)` prove `(38)` for every `Q>=6`.

For `1<=Q<=5`, use the `1296=6^4` labelled spanning trees of `K_6`.  Each
edge occurs in `432` of them, so their average weight is one third of the
total pair weight.  Pointwise,

```text
binom(k,2)>=k-1_(k>0),                                 (43)
```

where `k` is the number of active clauses.  Hence some tree has carrier
weight at least

```text
[6/7+epsilon_Q-average_j mu(U_j)]/3.                   (44)
```

The exact five-row census makes `(44)` strictly larger than `epsilon_Q`; its
smallest margin is

```text
956812421768354854111/11531953382864939723790          (45)
```

at `Q=1`.  A second path reconstructs all fifteen literal carrier edge
weights, all `2260` singleton laws, and the maximum tree for each of these
five levels.  It agrees with the direct average-union route and verifies the
actual consequence `(38)`, rather than only the lower-bound model.

Therefore the carrier-averaged spanning-tree criterion closes `(28)` for
every `Q>=1`.  This conclusion is redundant with `(36)` but identifies the
faithful intermediate object: not an unaddressed pair table and not a forced
tournament, but a weighted graphic-matroid carrier plus one endpoint word.

## 6. Exact verification and formal boundary

All three Python artifacts use `Fraction` arithmetic and `RuntimeError`
guards, so assertions remain active under `python -O`.  The final audit must
record byte-identical ordinary and optimized transcripts, source/output
SHA-256 values, and the pinned hashes of the THM-2941 carrier engine and
common-level reflected engine here before this candidate is promoted.

`TournamentH7.LRCTreeHunter.tree_hunter_add_le` formalizes the measure-
theoretic inequality `(10)` for an arbitrary parent-pointer tree.  The module
is imported by the `TournamentH7` root and contains neither `sorry` nor
`native_decide`.  Its role is partial: the reflected coordinate `(7)`, the
tent floor, the `186186`-row census, the endpoint word, and the all-`Q`
polynomial identity remain exact prose/computation, not Lean theorems.

## 7. Failure boundary and frontier

What is proved after promotion:

- the sufficient cone `(26)` for every nonconstant six-label reflected word;
- the one ray `(28)` for every positive integral `Q`;
- the carrier-tree certificate `(38)` on that same ray.

What is not proved:

- that a packet in `(27)` is realizable or dangerous;
- closure of every two-level, binary, or heterogeneous reflected word;
- a finite bound for arbitrary `k=1`;
- removal of a row from the current `k=2` or `k=3` ledgers;
- any improvement of their first-drift caps;
- LRC(14).
