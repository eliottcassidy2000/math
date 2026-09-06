---
id: THM-4434
title: "LRC14 universal scale-three network projection bound"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + FINITE-EXACT +
  INDEPENDENTLY AUDITED. For every primitive sorted distinct positive odd
  ternary-unit triple, the minimum of the three complete degree-zero raw
  network projections is at most 6/77, with equality only at (1,5,11).
  Outside the signed norm-four speed identities, all three projections are
  strictly below 6/77. The associated scale-three triple-comb completion is
  valid for bodies of Haar measure at least 6/77. This does not prove a
  quantitative ten-body Haar floor, arbitrary entry, synchronization, or
  LRC(14).
source: empty-core-sep06 + dense_geometry_dichotomy + projection_inequality + root, 2026-09-06
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
related:
  - THM-4425-lrc14-all-height-rank-one-carrier-closure
  - THM-4431-colored-lattice-basis-and-three-direction-lrc-network-closure
coefficient_script: 04-computation/lrc14_coefficient_box_empty_core_three_ray_sep06.py
coefficient_output: 05-knowledge/results/lrc14_coefficient_box_empty_core_three_ray_sep06.out
coefficient_script_sha256: 1552d098878b069c4f6e7e00737a982ef6b019105303bdb658a778cbe5a68ef7
coefficient_output_sha256: 09ecc3728b37540bbbde566e8ede59926afcdc8e394085aa15b5bd6082da82be
independent_coefficient_script: 04-computation/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.py
independent_coefficient_output: 05-knowledge/results/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.out
independent_coefficient_script_sha256: b326cf2ec6ae35aa4fd4c5434bab9028d17c7b21249047e5e34c01d76174c694
independent_coefficient_output_sha256: 7df7679bc41ae62b206b40da774d68ac09595e947cd59b66742a935e747effaa
literal_script: 04-computation/lrc14_universal_literal_empty_core_sep06.cpp
literal_output: 05-knowledge/results/lrc14_universal_literal_empty_core_sep06.out
literal_script_sha256: 64f2d209f5d93cb8f92250b8565c51a2a5437a61d7294f74bec9a712a502e5d0
literal_output_sha256: 49e334421cc50bc8768bbcc7d6457e63e86b89ddb106b27be4689690a99033b9
hash_basis: raw LF repository bytes
audit: >
  PASS. The analytic chain was independently checked from the zonotope
  determinants through section discrepancy, peak width, coefficient split,
  intercept, planar Minkowski bound, and the exact height cutoff. A separate
  import-free coefficient compiler reproduces all 308 patterns and the same
  semantic digest in 437,776 optimization-live gates. The native C++ head
  evaluates all 1,317,935 eligible triples through height 601 directly from
  six sheet permutations, without carrier congruences or roof formulas.
  A full H79 dump agrees with the older literal and raw implementations;
  root also compiled ordinary and optimized binaries and replayed that head
  byte-for-byte. All comparisons are exact integer or rational comparisons.
---

# THM-4434 -- LRC14 universal scale-three network projection bound

**PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + FINITE-EXACT +
INDEPENDENTLY AUDITED.** The universal degree-zero local network inequality is
closed. Chart entry, synchronization, and `LRC(14)` remain **OPEN**.

## 1. Statement

Let

```text
w=(a,b,c),       1<=a<b<c,                              (1)
```

be primitive, odd, and nonzero modulo three. Retain THM-4414's complete live
carrier dictionary

```text
Lambda(w)={C in Z^3:C dot w=0,
  every C_i nonzero mod 3,
  14|C_i|<3(w_j+w_k) for {i,j,k}={1,2,3}}.             (2)
```

For each coordinate, put

```text
E_i(w)=sum_(C in Lambda(w)) min(3/(7c),
  [3(w_j+w_k)-14|C_i|]/(14w_jw_k)).                    (3)
```

Then

```text
min_i E_i(w)<=6/77,                                    (4)
```

with equality if and only if `w=(1,5,11)`. If `w` has no signed primitive
relation with coefficient magnitudes `(1,1,2)`, then the stronger assertion
holds:

```text
E_i(w)<6/77 for every i.                               (5)
```

The excluded speed identities are exactly

```text
c=2a+b,       c=a+2b,       or       2b=a+c.           (6)
```

They are all covered by THM-4422, including the unique equality in (4).

THM-4414 also gives the sharp physical conclusion
`mu(F_w)<=min_i E_i(w)<=6/77`. Equality holds exactly at `(1,5,11)`,
where all three projections and the physical mass equal `6/77`. Common
odd ternary-unit dilation preserves physical Haar mass and permutes sheet
labels, so physical equality for nonprimitive tails is equivalent to
primitive reduction `{1,5,11}`. The norm-four exception is necessary:
`(1,5,7)` has projections `(8/245,6/49,4/35)`.

## 2. A relation-independent slice integral

The complete address construction is retained from the
[affine-slice proof](../../05-knowledge/results/lrc14_empty_core_certificate_sep06.md)
and its [actual-zero-coordinate addendum](../../05-knowledge/results/lrc14_pair_relation_empty_core_certificate_sep06.md).
Put `r=3/14` and let `v.w=0` be a primitive relation. If `u.w=1`,
any integer carrier has the lift `n=-u cross C`, so
`C=w cross n`. The strict roofs mean the three intervals
`((n_i-r)/w_i,(n_i+r)/w_i)` overlap pairwise and hence jointly. Thus
`e=n-yw` lies in the open error cube for some `y`. The integer defect
`delta=v.n=v.e` is lift-independent, and `v cross C=delta w` identifies
one affine carrier line of full primitive step `v`. The owner-residue
proof in the linked notes retains every defect and multiplier; no
primitive-direction cutoff is assumed.

Put `r=3/14` and let `v dot w=0` be a primitive nonzero integer relation.
Write `S=||v||_1` and `M=max_i|v_i|`. Choose `i` with `v_i!=0` and map the
cube `[-r,r]^3` linearly by

```text
e -> (v dot e,(w cross e)_i/v_i).                      (7)
```

Its image `Z` is a centrally symmetric planar zonotope. Let `bar f(delta)`
be the length of its closed section at first coordinate `delta`. On the
interior of its support put `f=bar f`; set `f=0` at the two support
endpoints and outside. This endpoint reset matches the strict physical defect
list and changes no integral. It matters only when an actual zero coefficient
gives a vertical endpoint edge.

The section length is even and concave on the interior of its support, hence
`f` is nonincreasing on the positive half-line. The three generator-pair
determinants of `Z` have absolute values `a,b,c`. Therefore the exact
zonotope area and Fubini give

```text
I=integral_R f(t)dt=4r^2(a+b+c)
                   =9(a+b+c)/49 <=27c/49.              (8)
```

For every even nonnegative function decreasing on the positive half-line,
rectangle comparison gives

```text
|sum_(k in Z) f(hk)-I/h|<=f(0),       h>0.             (9)
```

Indeed, the positive-half integral lies between the right and left
step-`h` sums; doubling leaves an error of at most the central value.

The exact owner-residue decomposition gives the slope load

```text
F_v(w)=(2/3)sum_k f(3k),                         if every v_i is a unit mod 3;
F_v(w)=(1/3)(sum_k f(k)-sum_k f(3k)),            otherwise.               (10)
```

The second line includes an actual zero integer coefficient. A primitive
relation cannot have two zero residues modulo three: the third would then
also vanish. Applying (9) with `h=1,3` in the correct upper/lower directions
gives in both cases

```text
F_v(w)<=(2/9)I+(2/3)f(0).                              (11)
```

At defect zero, `w cross e=t v`. In a coordinate where `|v_i|=M`,

```text
|t|M<=r(w_j+w_k),
f(0)<=2r(w_j+w_k)/M<=6c/(7M).                          (12)
```

Combining (8), (11), and (12) yields

```text
F_v(w)/c<=6/49+4/(7M).                                 (13)
```

For `M>=19`, this is at most `142/931<15/98`.

## 3. The complete small-coefficient box

It remains to prove `F_v(w)<=15c/98` when `M<=18`, outside norm four.
Up to signs and permutations, the complete coefficient universe is

```text
0<=p_1<=p_2<=p_3<=18,
at least two p_i nonzero, gcd(p_1,p_2,p_3)=1,
p_1+p_2+p_3 even,
at most one p_i zero mod 3,
p!=(0,1,1),(1,1,2).                                    (14)
```

The parity condition follows by reducing `v dot w=0` modulo two. The first
exclusion forces equal speeds; the second is the separately closed norm-four
case. There are exactly `308` patterns: `293` of full support and `15`
with an actual zero coordinate.

For each signed permutation, exact rational polygon clipping constructs every
cube section and the normalized speed polygon

```text
W_v={w in [0,1]^3:v dot w=0}.                          (15)
```

A section width is a maximum of linear functions minus a minimum, hence is
convex in `w`; its maximum on `W_v` occurs at a vertex. The compiler retains
the complete defect list, all sign sectors, and all vertices. An independent
cube-edge implementation checks full-support patterns, while the support-two
rectangle formula checks every actual-zero pattern.

Every one of the `308` exact maxima is at most `15/98`. Equality occurs
only for coefficient magnitudes `(1,7,8)` and is a boundary maximum of the
closed normalized speed polygon, not equality in (4). The excluded norm-four
pattern is a genuine slope hostile with value `2/7`. Together with (13),
this proves for every non-norm-four primitive relation

```text
F_v(w)<=15c/98.                                        (16)
```

## 4. Intercept and a short relation

Each unit-relation affine line contains two live longitudinal classes modulo
three and each one-zero line contains one. Exact open-interval counting gives

```text
N:=|Lambda(w)|<F_v(w)+B_v,
B_v<=2S/7+4/3.                                         (17)
```

For unit relations, the defects are multiples of three in
`|delta|<3S/14`, giving
`B_v=4|D|/3<4S/21+4/3`. For one-zero relations, deleting the multiples of
three gives `B_v=|D|<2S/7+4/3`. By (16),

```text
N<15c/98+2S/7+4/3.                                     (18)
```

The allowed defect list is nonempty: zero occurs in the unit case, while
an eligible nonunit relation has even norm at least six and permits `+/-1`.
Thus the strict summed count has no empty-list exception. The bound (18)
actually gives `N<2c/11` at the following threshold, which pays every
projection strictly.

Consequently THM-4422's automatic gate `N<=2c/11` holds whenever

```text
c>=(308/31)S+4312/93.                                  (19)
```

It remains to supply a relation meeting (19). Project the relation lattice to

```text
L={(x,y) in Z^2:ax+by=0 (mod c)},       det(L)=c.       (20)
```

The determinant follows because `gcd(a,b,c)=1`. The projected `l1` ball

```text
K_L0={(x,y):|x|+|y|+|(ax+by)/c|<=L0}                  (21)
```

is a centrally symmetric hexagon of area

```text
area(K_L0)=
 [2c(ab+ac+bc)/((a+b)(a+c)(b+c))]L0^2>(3/4)L0^2.      (22)
```

After scaling `c=1`, the positive numerator difference is

```text
8(ab+a+b)-3(a+b)(a+1)(b+1)
=3a(1-a)(b+1)+3b(1-b)(a+1)+2a(1-b)+2b(1-a).           (23)
```

Apply the elementary planar Minkowski argument to the interior of (21) with
`L0=4sqrt(c/3)`. Its area exceeds `4 det(L)`, so it contains a nonzero
lattice point. Dividing the corresponding relation by its content gives a
primitive relation of even norm

```text
S<4sqrt(c/3).                                          (24)
```

If it is norm four, THM-4422 applies. Otherwise (18)-(19) apply. For `c>=603`:

- if `S<=56`, the right side of (19) is at most
  `56056/93<603`;
- if `S>=58`, then
  `3S^2/16-(308/31)S-4312/93` is positive at `58` and increasing,
  while (24) gives `c>3S^2/16`.

Even parity exhausts the cases. Thus every eligible row with `c>=603`
satisfies (4), and every row outside (6) satisfies (5).

## 5. Complete height-601 head

The remaining exact universe is

```text
1<=a<b<c<=601,
a,b,c odd and nonzero mod 3,
gcd(a,b,c)=1.                                         (25)
```

An independently written native engine evaluates all `1,317,935` rows before
any support or relation classification. It constructs the three interval
sheets for each speed on exact denominator `42abc`, scans all six distinct
sheet assignments, and accumulates the three THM-4414 projections by
simultaneous interval cursors. It uses neither carrier congruences nor the
roof formula.

The exact result is

```text
rows                                             1,317,935
signed norm-four speed identities                    9,201
rows outside norm four                           1,308,734
native positive contacts                        147,048,282
min-projection failures                                  0
nonnorm-four rows with any E_i>=6/77                    0.  (26)
```

The unique equality is `(1,5,11)`. Outside norm four, the largest individual
projection is `12/161` at `(1,19,23)`; the largest minimum is `12/301`
at `(1,11,43)`. Eight hostile and positive controls check all three
projections, physical mass, and the factor-three correspondence between
native contacts and raw carriers.

A separate raw-congruence census checks the complete multi-direction head.
At height 79, every native row was also compared with the older pair-then-
third interval engine and the raw carrier formula. Root recompiled the native
source at ordinary and optimized settings; both height-79 outputs were
byte-identical. This completes (4)-(5) at every height.

## 6. Exact scale-three body consumer

For a finite positive-speed body `C`, set

```text
G_C={y in R/Z:||cy||>=1/14 for every c in C}.           (27)
```

Let `T={a,b,c}` be any three distinct positive odd integers prime to three.
If `G_C` is nonempty, then the universal local bound and THM-4414 give

```text
mu(G_C)>=6/77       implies       G_(3C union T) nonempty. (28)
```

Indeed, the preimage of `G_C` under multiplication by three is compact and
has the same Haar measure. Failure would place it inside the proper open
three-tail danger network, whose measure is at most (4). Strict measure
inequality is impossible; equality is also impossible because a nonempty
compact subset of a proper open subset of the connected circle cannot have
the same Haar measure. If the tails have a common divisor, multiplication by
that odd ternary unit reduces to the primitive case without changing Haar
measure or sheet labels up to permutation.

More sharply, it suffices that `mu(G_C)>=min_i E_i(T_0)`, where `T_0`
is the primitive tail reduction. Any failure in this exact typed domain
would force `mu(G_C)<mu(F_T)<=min_i E_i(T_0)<=6/77`. The
[consumer audit](../../05-knowledge/results/lrc14_universal_haar_consumer_empty_core_certificate_sep06.md)
retains the actual core, weak safety/strict danger endpoints, and the
THM-4032 hostile to lifting an arbitrary core witness.

The unqualified ten-body Haar floor is **REFUTED**: the recovered body
`{1,2,3,5,7,8,9,11,12,13}` has exact measure `14249/252252<6/77`.
A second inherited body defeats the necessary small-clock sieve while
remaining below that floor; its completed row has a safe phase. The
[body audit](../../05-knowledge/results/lrc14_haar_body_empty_core_sep06.md)
keeps the full closed geometry, including isolated safe points, and credits
the earlier numerical evidence. Completion therefore needs a stronger
actual-entry restriction or joint phase information when this gate fails.

Oddness also cannot be dropped: `(2,5,7)` has physical mass and minimum
projection `22/245>6/77`. The
[parity audit](../../05-knowledge/results/lrc14_parity_empty_core_sep06.md)
recovers the new norm-three obstruction and a separate mixed-parity
norm-four exception. The separate
[additive-family theorem](../../05-knowledge/results/lrc14_additive_parity_empty_core_sep06.md)
now proves the sharp replacement `6/55` for every primitive ternary-unit
`a+b=c`, with equality only at `(1,10,11)`. It does not close all mixed
parity. No arbitrary decomposition or synchronization theorem follows.

## 7. Reproduction and scope

```powershell
python -B 04-computation/lrc14_coefficient_box_empty_core_three_ray_sep06.py
python -B -O 04-computation/lrc14_coefficient_box_empty_core_three_ray_sep06.py
python -B 04-computation/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.py
python -B -O 04-computation/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.py
g++ -std=c++17 -O3 -DNDEBUG 04-computation/lrc14_universal_literal_empty_core_sep06.cpp -o lrc14_literal_sep06
./lrc14_literal_sep06 601
```

The coefficient-box producer performs `216,521` explicit gates plus
`1,172` comparison gates. The import-free coefficient referee adds `437,776`
exact gates. The native output is the complete exact transcript
of its bounded head. Full implementation details, controls, and independent
audit routes are recorded in
`05-knowledge/results/lrc14_global_slope_empty_core_certificate_sep06.md`,
`lrc14_coefficient_box_empty_core_three_ray_sep06.md`, and
`lrc14_universal_literal_empty_core_sep06.md`.

The universal degree-zero projection and typed triple-comb targets are proved.
Actual entry-produced body/phase control, arbitrary entry, synchronization,
and `LRC(14)` remain **OPEN**. The unqualified body floor is **REFUTED**.

## 8. Concurrent proof and reservation lineage

The empty-core session reserved this namespace honestly before the native
head was complete. Its root and two subagent referees audited the whole
analytic argument and compiler sources; the independently concurrent
[hexagon-session referee](../../05-knowledge/results/lrc14_universal_referee_overnight_hexagon_sep05.md)
reviewed the assembled candidate at `2e1ef24c4a`. That referee read the
frozen native output rather than rerunning the full census. Dense geometry
then supplied the separate import-free coefficient audit recorded above,
and the projection session promoted the same theorem concurrently with the
empty-core root. These are audits of one shared result, not separate novelty.
Both integrations retain the stronger physical equality and exact consumer.
No Lean formalization is claimed. The norm-54 wording in the earlier slice
note was corrected to an upper bound; the complete primitive box is unchanged.

An independent concurrent [universal audit](../../05-knowledge/results/overnight3_20260906_lrc_universal_audit.md) reproduces the coefficient box and literal head. Its [consumer reconstruction](../../05-knowledge/results/overnight3_20260906_lrc_consumers.md) recovers the THM-530 Section A body, retains all components and six isolated safe points, and proves exact adjacent-family oscillations. These are audits and refinements of the shared mechanisms.
