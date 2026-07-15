---
id: THM-802
title: Affine lifting of prefix-legal diagonal loops and noncentral packets on a fixed five-core-safe interval
status: PROVED (general affine phase-cell lifting lemma; existential realization of all five unequal one-fast-owner multiplicity classes; exact all-height wall order, prefix legality, reduced return, fixed metric/core incidence, and LRC lift) + VERIFIED (H=1 through 300; 406,350 walls and 406,650 chambers; finite-exact k=2 through 6 class audit; exhaustive 181,440-word legality census and tournament fingerprints)
source: codex-2026-07-14-S10 reduced-holonomy continuation
depends_on:
  - THM-779   # normalized collision transition and prefix-legality equation
  - THM-773   # prime-seven sheet tokenization and LRC lift
related: [THM-778, THM-783, THM-786, THM-788, THM-794, HYP-6840, MISTAKE-149]
verification:
  - 04-computation/lrc14_noncentral_diagonal_return_packets_codex_S10.py
  - 05-knowledge/results/lrc14_noncentral_diagonal_return_packets_codex_S10.out
---

# THM-802 — Affine lifting and noncentral diagonal-return packets

## 1. General affine phase-cell lifting lemma

Let `(d_a,beta_a)` be pairwise distinct pairs, with `d_a` positive and
`beta_a` nonzero modulo seven.  Put

```text
s_a=beta_a^(-1) (mod 7).
```

Suppose there are a real phase `alpha_0`, integers
`N_a=round(beta_a alpha_0)`, and a simple merge word `omega` of the limiting
wall times

```text
theta_(a,j)=(N_a+j+1/2-beta_a alpha_0)/d_a,
0<=j<d_a,                                               (L1)
```

such that:

1. no `beta_a alpha_0` is a half-integer and all times in (L1) are distinct
   in `(0,1)`;
2. the normalized token multiset `r_a=-s_a N_a`, modulo one common
   translation, is `{0,0,1,2,3,4,5,6}`;
3. `omega` is prefix-legal from `r` by THM-779's equation; and
4. there is one `c in F_7` with

   ```text
   d_a s_a=c       for every a.                         (L2)
   ```

Then there are an open phase interval `U` about `alpha_0` and `L_0` such that
for every multiple of seven `L>=L_0`, the positive speeds

```text
w_a=d_a L+beta_a                                  (L3)
```

are distinct and realize the word `omega` on every block
`[M/L,(M+1)/L]` whose two endpoints lie in `U`.  Every such block is fully
prime-seven blocked, and its normalized collision return is the identity.
Consequently the number of consecutive copies is `Theta(L)`.

If `U` is chosen inside a core-safe component for a fixed core `P`, the same
packet chain retains that metric/core incidence in the thirteen-speed LRC
row `7P union {w_a}`.

### Proof of the lifting lemma

At `x=M/L`, where `M/L` remains in a sufficiently small neighborhood on
which every nearest integer is the fixed `N_a`,

```text
round(w_a M/L)=d_a M+N_a.
```

The scaled position of its `j`th wall relative to the left boundary is

```text
L(z_(a,j)-M/L)
 = (N_a+j+1/2-beta_a M/L)/(d_a+beta_a/L).               (L4)
```

Formula (L4) converges uniformly to (L1) as `M/L -> alpha_0` and
`L -> infinity`.  Strictness of the finite set of endpoint and pairwise-order
inequalities therefore gives one open `U` and `L_0` on which the complete
wall word is still `omega`.  Shrinking `U` slightly ensures that both
endpoints of every retained block stay in the same nearest-integer cell, so
owner `a` has exactly `d_a` walls.  Since the pairs `(d_a,beta_a)` are
distinct, each speed difference is a nonzero affine function of `L`; enlarging
`L_0` makes every speed positive and avoids all of their finitely many
possible equality scales.

By (L2), the boundary tokens are

```text
-s_a(d_a M+N_a)=-cM-s_aN_a,
```

so their normalized state is the fixed `r`.  Prefix legality covers every
wall and intervening chamber.  After the packet, owner `a` has moved by
`-d_as_a=-c`, independent of `a`, hence the absolute return is diagonal and
the normalized return is the identity.  An interval of fixed positive length
contains `Theta(L)` consecutive mesh blocks.  Core incidence is unchanged
when `U` lies inside the core-safe set. ∎

This lemma is the structural result: a phase-realizable collision loop is an
inflatable object.  THM-794 is its equal-multiplicity instance.  The family
below proves that the inflatable language is already strictly larger.

## 2. Corollary — all five unequal one-fast-owner classes occur

Fix a five-owner core `P` and any nonempty open interval `J` contained in its
strict core-safe set.  For every

```text
k in {2,3,4,5,6},                                     (C1)
```

there are pairwise distinct integers `beta_a`, with

```text
d=(1,1,1,1,1,1,1,k),
beta_a=1 (mod 7) for 0<=a<=6,
beta_7=k (mod 7),                                      (C2)
```

such that the affine lifting lemma produces `Theta(L)` consecutive blocked
packets inside `J` for all sufficiently large multiples of seven `L`.  One
packet has owner word

```text
omega_k=(7,0,1,2,3,4,5,6, 7,7,...,7),                 (C3)
```

with `k-1` final occurrences of owner seven.  Its owner-count vector is
`(1,1,1,1,1,1,1,k)`, and its reduced return is the identity.  Thus every
unequal residue multiplicity `k=2,...,6` occurs in a phase-realizable
diagonal-loop language, already after restricting to one faster owner and
seven unit-rate owners.

This corollary is existential.  It asserts no effective upper bound on the
chosen `|beta_a|` or on the onset scale `L_0`.

### Proof

Choose an irrational `alpha_0 in J`, fix `k` as in (C1), and put

```text
u=k^(-1) (mod 7),        theta_0=1/(4k).
```

Prescribe the `k` limiting walls of owner seven at

```text
theta_0+h/k,       0<=h<k,                              (C4)
```

and prescribe the singleton wall of owner `j`, `0<=j<=6`, at

```text
theta_j=theta_0+(j+1)/(8k).                             (C5)
```

The first fast wall occurs at `2/(8k)`, the seven singleton walls at
`3/(8k),...,9/(8k)`, and the second fast wall at `10/(8k)`.  All later fast
walls follow in order and remain below one.  Hence the strict merge word is
exactly (C3).

Choose desired nearest-integer residues

```text
N_7=0 (mod 7),       N_j=u+j (mod 7), 0<=j<=6.          (C6)
```

The corresponding normalized initial offsets are

```text
r_7=0,               r_j=-(u+j) (mod 7).               (C7)
```

The seven singleton offsets in (C7) run through all of `F_7`, and owner seven
duplicates the zero, so the initial state is covered.  At the first fast wall
the prefix sum is zero.  Before singleton `j` it is `u+j`, cancelling
`r_j`.  After all seven singletons it is again `u`; the one previous fast
event contributes the same `u`, so every remaining consecutive fast event is
legal.  Finally

```text
1*1=1,            k*u=1 (mod 7),                        (C8)
```

which proves the diagonal return.

It remains only to realize the prescribed phases with integers satisfying
(C2).  For each fixed nonzero residue `b mod 7`, the forward orbit

```text
{(b+7n)alpha_0 mod 7 : n>=0}                            (C9)
```

is dense in `R/7Z`, because its rotation step is `7alpha_0` on a circle of
length seven and `alpha_0` is irrational.  For a singleton owner, choose a
small target neighborhood of

```text
N_j+1/2-theta_j;
```

for owner seven choose one about

```text
N_7+1/2-k theta_0=N_7+1/4.
```

The neighborhoods can be chosen disjoint from all relevant half-integer
boundaries and small enough to preserve the strict order (C4)--(C5).  Every
target neighborhood is hit infinitely often, so the choices below can also
be made pairwise distinct.  Density (C9), with
`b=1` for the singleton owners and `b=k` for owner seven, supplies distinct
positive choices `beta_a=b+7n` hitting those neighborhoods.  Their actual
nearest integers have the required residues (C6), even though their absolute
integer lifts are immaterial.  All hypotheses of the affine lifting lemma now
hold.  Shrink its phase interval `U` inside `J`; the resulting packet count is
`Theta(L)` and its metric lift stays core-safe. ∎

The five rows in (C1) are the nonzero unequal count residues left after a
same-owner block of seven events is treated as its own zero-return refinement.
The next subsection supplies an explicit quantitative representative of the
first row `k=2`.

## 3. Explicit unequal-multiplicity family

For every integer `H>=1`, put

```text
L=105H,
d=(1,1,1,1,1,1,1,2),
beta=(1,22,-20,-6,15,-13,8,2),
w_a=d_a L+beta_a,
W=(L+1,L+22,L-20,L-6,L+15,L-13,L+8,2L+2).       (1)
```

The eight speeds are positive and distinct, none is divisible by seven, and
their inverse residues are

```text
s=(w_a^(-1) mod 7)=(1,1,1,1,1,1,1,4).           (2)
```

Let

```text
I=[2/15,1/7],
N=(0,3,-3,-1,2,-2,1,0).
```

For every integer block index

```text
14H <= M <= 15H-1,                                (3)
```

the complete global midpoint-wall word in `[M/L,(M+1)/L]` is simple and has
owner word

```text
omega=(7,2,5,3,0,6,4,1,7).                        (4)
```

More precisely, the walls of owner `a` in the block are

```text
z_(a,j)(M)=(d_a M+N_a+j+1/2)/(d_a L+beta_a),
0<=j<d_a,                                          (5)
```

and (4) means

```text
M/L
 < z_(7,0) < z_(2,0) < z_(5,0) < z_(3,0)
 < z_(0,0) < z_(6,0) < z_(4,0) < z_(1,0)
 < z_(7,1) < (M+1)/L.                              (6)
```

Every point of the **fixed interval** `I`, including all its walls, is
prime-seven-sheet blocked by `W`.  Thus `I` contains `H` consecutive packets,
`9H` covered walls, and `8H` owner switches.  At every block boundary the
normalized collision state is

```text
r=(0,4,3,1,5,2,6,0),                               (7)
```

and one packet returns it exactly.  Indeed its count vector is `d`, and

```text
d_a s_a = 1 (mod 7)       for every owner a.         (8)
```

So the raw token return is diagonal translation by `-1`, hence zero in
`F_7^8/Delta`.  This packet is **noncentral in word/multiplicity shape**—owner
seven occurs twice while the others occur once—even though its return
holonomy is central.  It is not a repetition or relabelling of THM-794's
once-per-owner packet.

The metric sidecar is genuine LRC14 core incidence.  For

```text
P={1,2,11,12,13},                                    (9)
```

the whole interval `I` is strictly core-safe at threshold `1/14`.  Hence the
thirteen nonzero speeds

```text
A_H=7P union W                                       (10)
```

are distinct, and for every `x in I` each of the seven lifted times
`t=(x+k)/7`, `k in F_7`, is obstructed: the five core speeds are safe while
some member of `W` is strictly dangerous.  This is a fixed obstructed region,
not a counterexample to LRC14; witnesses may occur outside these lifts.

Finally the fastest and second-fastest owners are

```text
f=2L+2,       g=L+22,       ceil(f/g)=2.              (11)
```

The span from the first to the last displayed wall is

```text
E_H=(2H-1)/(210H+2),                                  (12)
```

and for every `H>=3`,

```text
E_H > 1/g+2/f.                                        (13)
```

Thus THM-794's failure of the universal extent target is not confined to
uniform once-per-owner packets.

## 4. Proof of the explicit family

### 4.1. The boundary phase cell and the complete event list

Every block boundary in (3) lies in `I`, because

```text
14H/L=2/15,          15H/L=1/7.
```

For `x in I`, direct interval arithmetic gives

| owner `a` | `beta_a` | range of `beta_a x` | nearest integer `N_a` |
|---:|---:|---:|---:|
| 0 | 1 | `[2/15,1/7]` | 0 |
| 1 | 22 | `[44/15,22/7]` | 3 |
| 2 | -20 | `[-20/7,-8/3]` | -3 |
| 3 | -6 | `[-6/7,-4/5]` | -1 |
| 4 | 15 | `[2,15/7]` | 2 |
| 5 | -13 | `[-13/7,-26/15]` | -2 |
| 6 | 8 | `[16/15,8/7]` | 1 |
| 7 | 2 | `[4/15,2/7]` | 0 |

Every displayed range lies strictly inside `(N_a-1/2,N_a+1/2)`.  At
`x=M/L`, therefore,

```text
round(w_a x)=d_a M+N_a.                              (14)
```

The same statement at `(M+1)/L` differs by exactly `d_a`.  Hence owner `a`
has exactly `d_a` walls in the open block, and they are precisely (5).

Cross-multiplication reduces all comparisons in (6) to four elementary
conditions:

```text
M/L < z_(7,0)                         iff 4M<L;
z_(7,0) < z_(2,0)                    iff 84M>11L-10;
z_(a,0) < z_(b,0),
 (a,b)=(2,5),(5,3),(3,0),(0,6),(6,4),(4,1)
                                         iff 14M<2L-5;
z_(1,0) < z_(7,1)                    iff 84M>11L-52.
```

The final comparison `z_(7,1)<(M+1)/L` is automatic for positive `L,M`.
For `L=105H` and (3),

```text
4M < 60H < 105H,
84M >= 1176H > 1155H-10,
14M <= 210H-14 < 210H-5.
```

The weaker fourth inequality follows as well.  This proves (6), simplicity,
and completeness of the global event list.

### 4.2. Prefix legality and the nonuniform reduced return

At `x=M/L`, combine (2), (8), and (14):

```text
k_a=-s_a(d_a M+N_a)=-M-s_a N_a (mod 7).              (15)
```

After subtracting the common duplicate root `-M`, (15) is exactly (7).  Its
multiset is `{0,0,1,2,3,4,5,6}`, so the starting chamber is covered.

For a word prefix, let `S_j` be the sum of inverse steps already crossed and
`C_a(j)` the number of earlier occurrences of owner `a`.  THM-779's exact
prefix condition is

```text
r_a+S_j-C_a(j)s_a=0                                  (16)
```

when the next owner is `a`.  Along (4), the right required initial residues
`C_a(j)s_a-S_j` are

```text
(0,3,2,1,0,6,5,4,0),                                 (17)
```

in owner order `(7,2,5,3,0,6,4,1,7)`.  These are exactly the corresponding
entries of (7).  Thus every prefix is legal: the walling owner holds one copy
of the duplicate root, the other seven tokens form `F_7` at the wall, and the
next chamber remains covered.

At the end of the packet the total inverse step is

```text
7*1+2*4=15=1 (mod 7),
```

while (8) says that every owner has also accumulated step one.  The exact
state update

```text
r'_a=r_a+S_9-d_a s_a
```

therefore gives `r'=r`; in absolute coordinates every token changes by `-1`.
Induction across the `H` blocks proves coverage of every wall and chamber in
`I`, including its two generic endpoints.

### 4.3. Metric/core incidence and the LRC lift

On `I`, the exact minimum clearances of the five core speeds are

```text
p             1       2       11      12      13
min ||px||    2/15    4/15    3/7     2/7     1/7.      (18)
```

For example `11I=[22/15,11/7]` straddles the half-integer but no integer, so
its minimum occurs at an endpoint; the other four rows are monotone inside a
single nearest-integer cell.  Every value in (18) exceeds `1/14`.

For a lifted time `t=(x+k)/7` and a core speed `7p`,

```text
||7p t||=||p(x+k)||=||px||.
```

The core is therefore safe.  Prime-seven blocking by `W` says exactly that
for each `k` some exception satisfies `||w_a(x+k)/7||<1/14`.  This proves the
obstruction assertion for (10).  Since every member of `W` is nonzero modulo
seven, it is distinct from `7P`; direct inspection of (1) gives thirteen
distinct positive speeds.

### 4.4. Ratio and extent

The largest `beta` among the first seven owners is `22`, so (11) follows, and

```text
1 < (2L+2)/(L+22) < 2.
```

The first wall is `z_(7,0)(14H)` and the last is
`z_(7,1)(15H-1)`.  Subtracting their formulas gives (12).  Substitution in
(13), followed by positive-denominator cross-multiplication, verifies it for
`H=3`; its numerator is

```text
210H^2-481H-68 = 210u^2+779u+379,       u=H-3,
```

which is positive for every `u>=0`.  This finishes the proof. ∎

## 5. The computable reduced state that survives this audit

The construction isolates a minimal useful packet certificate:

```text
Xi_M=(s, r, omega, d;
      exact wall coordinates z_(a,j)(M),
      [M/L,(M+1)/L] subset I subset G_P).              (19)
```

- `(s,r,omega)` decides every prefix through (16);
- `(d_a s_a)_a` decides the reduced return;
- the exact wall coordinates decide centered-Beatty realizability; and
- the two inclusions retain metric/core incidence.

The packet-level collision vertex in (19) never changes: all `H` packets are
self-loops in one normalized state, while the metric base advances by `1/L`.
Thus merely counting SCC exits is also insufficient unless **all** legal
diagonal-return loops—not only THM-794's once-per-owner generator—are
contracted and their metric lift remains attached.

The owner multiplicity vector is a particularly sharp lossy quotient.  All

```text
9!/2 = 181,440
```

words with multiplicity `d` have the same reduced-return congruence (8), but
an exhaustive application of (16) finds only three legal from (7):

```text
(0,6,4,1,2,5,3,7,7),
(7,2,5,3,0,6,4,1,7),
(7,7,6,4,1,2,5,3,0).                                 (20)
```

Even retaining the first-occurrence order
`(7,2,5,3,0,6,4,1)` leaves eight possible insertion slots for the second
owner-seven event; only the displayed slot is legal from (7).  Thus
multiplicity preserves return but destroys prefix legality, and the
first-owner tournament still destroys the decisive repeated-owner insertion.

## 6. Tournament Analysis

Three vertex choices expose the quotient boundary.

1. **Wall occurrences.**  The nine events in one packet are vertices, with
   chronological precedence as pairwise observable and time reversal as the
   switch.  The tournament is `T_9`: score histogram `{0:1,...,8:1}`, no
   directed triangles, nine singleton SCCs, and one Hamiltonian path—the raw
   word (4).  It is lossless for order but forgets that the endpoint state
   returns.
2. **First owner occurrences.**  Chronology gives `T_8`; the count vector is a
   sidecar.  This retains the first-owner path but loses the insertion slot of
   the second owner-seven event, as the `8`-versus-`1` census above proves.
3. **Owners in collision-sheet gauge.**  Use the eight owners as vertices.
   The marked observable compares `(r_a, first-occurrence position)`
   lexicographically, resolving the equal-sheet pair `(7,0)` by chronology.
   Switch to circular orientation

   ```text
   a -> b iff r_b-r_a in {1,2,3} (mod 7),
   ```

   retaining the same tie rule.  The marked tournament is transitive with
   score histogram `{0:1,...,7:1}`, no directed triangles, singleton SCCs,
   one Hamiltonian path `(7,0,3,5,2,1,4,6)`.  The circular switch has score
   histogram `{3:4,4:4}`, `20` directed triangles, one eight-owner SCC, `629`
   Hamiltonian paths, and lexicographic tie path
   `(0,2,1,4,6,7,3,5)`.  The gauges differ on `9` edges.

The circular tournament sees the collision geometry, and the occurrence
tournament sees exact chronology; neither alone preserves both prefix
legality and reduced return.  In particular, its SCC is **telemetry**, not the
LRC predicate: the predicate is exact owner-word prefix legality in the
centered mechanical order, coupled to the reduced return and metric/core
incidence.  The decorated carrier (19) retains precisely those fields.
