---
id: THM-2174
title: "Phase-danger clutters, finite-window lifts, and the endpoint-scale obstruction"
status: >
  PROVED + VERIFIED-EXACT. Exact-denominator deletion success is generated
  by the antichain of minimal numerator danger sets; for Q<=14 this collapses
  to the divisor mask, while genuine metric residue geometry begins at 15.
  Matching h-digit boundary windows makes a THM-2171 pump preserve every
  denominator dividing q^h, but digitwise zero masks already fail at 25.
  At unbounded scale, a fixed core's endpoint phase word is periodic in W
  while its signed cocycle is C_r/W. An explicit thirteen-speed pump
  preserves the full endpoint phase word modulo 1680 and three independent
  crossing relations while changing exact safe measure. Thus both aligned
  finite residue words and scale-sensitive endpoint current are real phase
  coordinates; no fixed window proves LRC(14).
source: codex-2026-07-24-LRC-post-2168-closure
depends_on:
  - THM-2162
  - THM-2163
  - THM-2171
related:
  - THM-2161
  - THM-2166
  - THM-2167
  - THM-2168
  - THM-2169
script: 04-computation/lrc14_endpoint_phase_scale_obstruction_thm2174.py
output: 05-knowledge/results/lrc14_endpoint_phase_scale_obstruction_thm2174.out
script_sha256: 561050e451f5d60526381e3c750da9fd45904cac9ce4af92a9263682bba552ce
output_sha256: ba4c1eba837c247b64ce6d6c407dbf1f49b36bb2edd4e46612fc7a319c3c2e49
independent_script: 04-computation/lrc14_phase_danger_clutter_referee_thm2174.py
independent_output: 05-knowledge/results/lrc14_phase_danger_clutter_referee_thm2174.out
independent_script_sha256: 96b96fa6bd4e1d29e3cfba801562d70a065740944c7ca21566e4ce24e0fc4081
independent_output_sha256: ed37d444f6690b2a2fbc2c7d39d476dc49690e23e5ed11f71be12c7b4e7e5942
hash_basis: working-tree bytes (LF)
---

# THM-2174 -- endpoint-phase scale obstruction

At radius `1/14`, write

```text
G_E={t in R/Z:||et||>=1/14 for every e in E}.         (1)
```

THM-2171 makes bounded relation systems pumpable without losing positivity,
order, distinctness, or primitivity. This theorem proves that even the
natural endpoint-phase sidecar from THM-2162 is not enough to preserve the
LRC target: its phase labels are finite, but its signed magnitude remembers
the unbounded speed scale.

## 1. Endpoint phases versus endpoint current

Let the positive-length components of `G_E` be

```text
I_s=[l_s,r_s],                    s=1,...,K,           (2)
```

with harmless choices at their measure-zero endpoints. Choose an integer
`L_E` divisible by every denominator of every `l_s,r_s`. Let `H` be the
continuous periodic primitive from THM-2162:

```text
H'=1_(-1/14,1/14)-1/7,             H(0)=0.            (3)
```

For a residue `r mod L_E`, define the endpoint numerator

```text
C_E(r)=sum_(s=1)^K [H(r r_s)-H(r l_s)].              (4)
```

THM-2162's exact endpoint identity immediately gives:

> **Phase-scale lemma.** If `W` is a positive integer and
> `W=r mod L_E`, then
>
> ```text
> epsilon_W(G_E)
> :=|G_E intersect {t:||Wt||<1/14}|-|G_E|/7
> =C_E(r)/W.                                         (5)
> ```

Indeed, `Wl_s` and `rl_s` agree modulo one, as do `Wr_s` and `rr_s`,
while the Jacobian factor in the endpoint identity remains `1/W`.

This separates two coordinates that a finite phase word conflates:

```text
endpoint phase label:       r=W mod L_E;
signed endpoint magnitude:  C_E(r)/W.                (6)
```

If `C_E(r)!=0`, then the speeds `r+nL_E` have one common endpoint phase
word but infinitely many distinct endpoint currents. More sharply, inside
that residue class,

```text
C_E(r)/W=C_E(r)/W'  iff  W=W'.                       (7)
```

Thus preserving the exact current on a nonzero class retains the speed
scale itself.

## 2. A seven-speed core with nonzero phase currents

Take

```text
E=(1,2,3,4,5,6,8).                                   (8)
```

Every boundary of `G_E` is a danger-comb endpoint

```text
(14m+/-1)/(14e),                  e in E,             (9)
```

so one may take

```text
L_E=lcm(14e:e in E)=1680.                            (10)
```

Exact interval arithmetic gives

```text
|G_E|=27/70,             K=16,                       (11)
C_E(1)=-27/490,          C_E(2)=-27/245.             (12)
```

Both residue classes therefore have the infinite scale sensitivity in
(7).

## 3. The relation- and phase-preserving pump

For `a=1,2,3`, put

```text
f_(a,-)=3360a+1,             f_(a,+)=3360a+2,
f'_(a,-)=1680a+1,            f'_(a,+)=1680a+2.       (13)
```

The two ordered thirteen-speed rows are

```text
V =E union (3361,3362,6721,6722,10081,10082),
V'=E union (1681,1682,3361,3362,5041,5042).          (14)
```

Both are primitive because they contain `1`. Apply the base-two THM-2171
pump which deletes radix level `j=4` and splices in level `k=5`. The core
coordinates are below `2^4` and remain fixed, while

```text
3360a+r
 =r+2^5(105a)
 |-> r+2^4(105a)
 =1680a+r,                       r=1,2.              (15)
```

Thus the pump sends `V` exactly to `V'`. Its two boundary quotients are

```text
Z_4=(0,0,0,0,0,0,0,210,210,420,420,630,630),
Z_5=(0,0,0,0,0,0,0,105,105,210,210,315,315).        (16)
```

They have the same owner suffix `{8,...,13}` and the same one-based
quotient-tie cut mask `{7,9,11}`.

There are three independent height-one crossing relations. With `e_i`
denoting the coordinate vectors of the ordered row, set

```text
rho_a=-e_1-e_(6+2a)+e_(7+2a),        a=1,2,3.        (17)
```

For each pair,

```text
-f_(a,-)+f_(a,+)=1,
-f'_(a,-)+f'_(a,+)=1,                              (18)
```

and coordinate `1` pays `-1`. Hence every `rho_a` annihilates both rows.
Their private far-coordinate pairs make them independent over every field.
Equation (16) also shows that their carries at levels four and five are
both zero, so (17) is preserved by the exact THM-2171 splice.

These relations have a stronger version of the coefficient shape forced
by THM-2166:

```text
cut carry=1,
far height=1 with every nonzero far coefficient a 7-unit,
core support=1 and core height=1.                    (19)
```

The rows in (14) are not zero-safe, so (19) is a carrier stress test, not
an invocation or converse of THM-2166.

Finally,

```text
f_(a,+/-)-f'_(a,+/-)=1680a.                          (20)
```

Thus the pump preserves every coordinate modulo `L_E=1680`. In particular,
it preserves the complete endpoint phase word of every far comb against
every component of the fixed core `G_E`; this is exactly THM-2171's optional
finite-residue sidecar with `N=1680`.

## 4. The exact current and target still move

The far speeds in each minus column have residue one modulo `1680`, and
those in each plus column have residue two. Equations (5) and (12) give

```text
epsilon_W(G_E)=-27/(490W),       W=1 mod 1680,
epsilon_W(G_E)=-27/(245W),       W=2 mod 1680.        (21)
```

Therefore every paired old/new far comb has identical endpoint phases but
different exact current.

The failure persists for the full thirteen-speed target, not only for a
one-comb marginal. Two independent exact interval calculations give

```text
|G_V|
 =1561405750435498559/10390707539702618590,

|G_(V')|
 =317645844187362436/2113439446871390435,             (22)
```

and

```text
|G_(V')|-|G_V|
 =126112336463776271349/4405994577616688706478598
 >0.                                                  (23)
```

One evaluator merges the complete family of rational danger intervals.
The other forms the complete rational boundary arrangement and sums exactly
the safe open cells. They agree on both fractions in (22). Normal and
optimized Python produce the stored transcript.

## 5. Consequence and exact residual

The example preserves more than the algebraic state of THM-2171:

```text
fixed core;
positive strict order and primitivity;
owner and quotient-tie state;
three independent height-one crossing relations and their carries;
the full coordinatewise residue word modulo the core endpoint modulus;
every resulting endpoint phase label.               (24)
```

Nevertheless the exact safe measure changes. Therefore:

1. endpoint phase labels, even at the full core denominator, are not a
   target-preserving pump sidecar;
2. the missing datum already appears in the one-comb THM-2162 current as
   the scale factor `1/W`;
3. retaining the exact nonzero current on a fixed phase class is not finite
   state compression--by (7), it retains `W` itself.

This sharpens the post-THM-2171 residual. Algebraic feasibility has a finite
terminal, and any chosen finite residue bank can be preserved, but a
counterexample descent still needs a scale-sensitive phase/current theorem
that preserves the **zero-safe predicate**, not merely its endpoint labels.
The theorem does not say that no other finite target invariant can exist,
does not close THM-2168's `(0,0)`, `(0,1)`, `4+3`, or `5+3` profiles, and
does not prove or refute LRC(14).

## 6. Exact-denominator danger clutters

The endpoint witness proves that no fixed residue word is the whole target.
Within one fixed denominator, however, there is an exact finite carrier.
For an integer residue put

```text
|x|_Q=min(x mod Q,-x mod Q),
B_(Q,a)(v)={i:14|a v_i|_Q<Q},      a in (Z/QZ)^*.    (25)
```

If `S` is a set of deleted labels, the remaining row is safe at the reduced
phase `a/Q` exactly when

```text
B_(Q,a)(v) subset S.                                  (26)
```

It admits some exact-denominator-`Q` phase exactly when `S` contains one
edge of the **phase-danger clutter**

```text
A_Q(v)=Min_subset{B_(Q,a)(v):a in (Z/QZ)^*}.          (27)
```

This is the unique antichain of minimal generators of the corresponding
upper ideal in the Boolean deletion lattice. It is the canonical
existential carrier. To remember which numerators work, retain the labelled
profile `a |-> B_(Q,a)` instead.

The object is a hypergraph, not a tournament: edges can have every
cardinality and containment, not pairwise orientation, is the target
operation.

## 7. Divisor masks are complete exactly through denominator 14

The danger inequality in (25) is

```text
|a v_i|_Q<ceil(Q/14).                                 (28)
```

For `Q<=14`, only residue zero is dangerous, so every unit numerator has
the same danger set

```text
E_Q(v)={i:Q divides v_i},             A_Q(v)={E_Q(v)}. (29)
```

For the original thirteen-speed row, the target is the bit `E_Q=empty`.
For it together with all single deletions, the exact quotient has fifteen
values:

```text
empty;          singleton {i}, 1<=i<=13;          plural. (30)
```

Arbitrary deletion recursion needs the full mask, since deleting `S`
succeeds exactly when `E_Q subset S`.

This mask is not contained in the rank-two carry-owner state. At base five
let

```text
A=(6,12,18,13,23,21,26,41,51,22,32,37,48),
B=(6,12,18,13,23,20,27,41,51,22,32,37,48),           (31)

r=(1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1),
s=(1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1).               (32)
```

Both rows are primitive and distinct, both relations vanish and have rank
two modulo five, and their complete carry pairs and owner masks agree:

```text
K_j=(0,0),(-1,1),(-1,1),(0,0),
O_j=[13],[13],{7,8,9,11,12,13},empty.                (33)
```

But `E_5(A)=empty`, `E_5(B)={6}`, and their integer margins are
`m_5(A)=1`, `m_5(B)=0`.

This loss is sharp over one carrier fibre. Fix

```text
Z_1=(1,2,3,2,4,4,5,8,10,4,6,7,9)
```

and enumerate primitive distinct `v=5Z_1+D` with
`(r;s)D=(-5,5)^T`. Every row has the same complete carry-owner path.
Exact dynamic programming gives

```text
E_5=empty:                  76787,
|E_5|=1:                   346271,
|E_5|>=2:                 1783379,
total:                    2206437,
distinct masks:             1836.                    (34)
```

Every singleton occurs, so all fifteen values in (30) occur over one
carry-owner fibre.

Along a base-`q` radix, the divisor-mask regime persists while `q^h<=14`.
The complete sidecar is the nested prefix-zero mask

```text
F_h={i:D_(0,i)=...=D_(h-1,i)=0}={i:q^h divides v_i}. (35)
```

The first metric depths for `q=2,3,5,7,11,13` are respectively
`4,3,2,2,2,2`. At those depths nonzero residues enter the danger interval.

## 8. Fixed windows are sufficient, separate zero masks are not

Suppose a THM-2171 carry-owner-tie state repeats at boundaries `j<k`.
Fix `h>=1`. If `j<h`, also require equality of aligned windows

```text
(D_j,...,D_(j+h-1))=(D_k,...,D_(k+h-1)).             (36)
```

If `j>=h`, no window is needed. After deleting levels `j,...,k-1`, all
digits below `j` are unchanged; (36) replaces each remaining residue digit
by an equal one. Hence the unnormalized pumped row satisfies

```text
v'=v mod q^h                                         (37)
```

coordinatewise. It preserves every labelled profile, clutter, numerator
margin, and deletion phase predicate for every `Q|q^h`. This is an exact
target-preserving lift at fixed depth. Its raw window alphabet is
`q^(13h)`, and normalization should be postponed until after any labelled
residue certificate.

Digitwise zero masks do not fake the aligned window. With (32), put

```text
U=(6,11,16,12,21,22,26,42,51,23,32,37,49),
W=(6,11,16,12,21,22,26,41,51,23,31,37,49).           (38)
```

They share their maximum, both relations, the entire carry-owner path, and
every base-five zero-digit mask. For

```text
m_Q^*(v)=max_((a,Q)=1) min_i |a v_i|_Q,
```

exact computation nevertheless gives

```text
m_25^*(U)=1,
m_25^*(W)=2,       successful W numerators=3,7,18,22. (39)
```

Their unrestricted margins are both `m_25=5`, because those also allow a
numerator divisible by five. Their clutters expose the missing angular word:

```text
A_25(U)={
 {1},{2},{3},{5},{6},{8},{10},{11},{4,12},{7,9,13}
},
A_25(W)={empty}.                                      (40)
```

Thus adjacent residue words, not coordinatewise zero events, are the
finite-depth phase coordinate.

## 9. A positive scale consequence

The scale obstruction also yields a useful bound. Let a deletion core
`C` have safe-set mass `mu>0` and `K` positive-length BV components
`[l_t,r_t]`. If adjoining speed `W` makes the full row zero-safe, then
THM-2162 and equation (5) give

```text
sum_t[H(Wr_t)-H(Wl_t)]=6Wmu/7>0.                     (41)
```

Since `max H-min H=6/49`,

```text
W<=K/(7mu).                                          (42)
```

So a fixed positive-mass deletion core cannot support arbitrarily large
zero-safe extensions. The current (41) is the archimedean complement of
the rational clutters (27): finite windows preserve bounded denominator
geometry, while nonzero endpoint current retains scale itself.

The two companion computations independently verify (8)--(24) and
(31)--(40), normally and under `-O`. What remains open is to force from
zero-safety either bounded phase depth or a repeated state carrying enough
endpoint current. QED.
