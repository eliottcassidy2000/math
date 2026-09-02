---
id: THM-4335
title: "LRC(14) owner permutations, component addresses, and minority renewal"
status: >
  PROVED ELEMENTARY OWNER/COMPONENT/RENEWAL LEMMAS + PROVED RELATIVE TO
  CITED LRC THROUGH THIRTEEN TOTAL RUNNERS FOR THE SCALE-FREE COROLLARY +
  VERIFIED-EXACT CONTROLS + INDEPENDENT INTERVAL AUDIT; LRC(14) OPEN. At an
  equality-capacity row dC union T, 2<=d<=7, with exactly d distinct d-unit
  tails, strict failure is equivalent to an owner permutation at every
  core-safe quotient phase. On each positive component, owners have constant
  affine tooth addresses; endpoint denominators and pairwise determinants
  give sharper necessary width inequalities and sufficient escape tests.
  Initial-segment selected-tail thresholds become (17,17,17,16,13,12). A
  p=5 h=420 row
  passes the inherited marginal gates but is closed by owner disappearance.
  Separately, every physical minority-anchor component has an exact shortest
  tooth-cover renewal budget with odd endpoint tolls and gcd-quantized
  overlap tolls. This is an obstruction/certificate compiler, not a closure
  of the minority branch or LRC(14).
source: root + minority_owner_completion + owner_walk_wildcard / LRC14 continuation session, 2026-09-02
depends_on:
  - LRCUpTo13
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
related:
  - THM-592-radius-derivative-structure
  - THM-771-sheet-endpoint-defect-and-reduced-winding-pierce
  - THM-773-prime-seven-sheet-monodromy-and-tournament-fibre
  - THM-1156-tooth-seam-chi7-bipartition
  - THM-1244-slowest-spoke-component-handoff-debt
  - THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy
  - THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape
component_script: 04-computation/lrc14_unit_component_address_escape_probe_20260902.py
component_output: 05-knowledge/results/lrc14_unit_component_address_escape_probe_20260902.out
component_independent_script: 04-computation/lrc14_unit_component_address_escape_probe_20260902_independent_audit.py
component_independent_output: 05-knowledge/results/lrc14_unit_component_address_escape_probe_20260902_independent_audit.out
component_script_sha256: 8b8e92c9bd0647cccdb6a78729d3b88a4b9447053d14893e28de49a519cf2111
component_output_sha256: 77d6f33a336cddab5b426fd1912343fef3245db1d52fdd1b0051cd251d3224fd
component_independent_script_sha256: 5b36be9aea4ac0d6d04b07ac16ec1e530cf0619cc39a57f0510bacc1e0e4f9c1
component_independent_output_sha256: 29359f4cd3e04bcd2e1aabb0ea799c25aaf1e88dbc5888894ad46c2c1ef2a575
owner_script: 04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py
owner_output: 05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out
owner_script_sha256: cb4f7118398c0b8f2049669742fe18b99688dc679bd72f01b59e444cea77ae79
owner_output_sha256: ff3df7ff3d6412a9597e1acf63bc4d32e24ffb2e53d6d603c8510d3228917383
renewal_script: 04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py
renewal_output: 05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out
renewal_script_sha256: 18f42dfc2ed67706cc61e84c1a03fab937d357a41afc7d8fa0c419402e3fd083
renewal_output_sha256: 57821c6521a081943a5be377197f819b69570fdb611caacc007d12f45621240e
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The owner probe checks 710,910 exact requirements, including
  72,044 rational owner formulas, every strict equality boundary and d=7
  half-integer tie, 47,716 gauge shifts, 40,690 pair-tooth intersections,
  two independent divisor constructors, and the complete p=5 control. The
  component probe and independent successive-intersection implementation
  verify 159,642 lift/activity cases, 1,418 actual effective-tooth
  containments, all six initial-segment rows, the endpoint strictness
  hostile, and the p=5 component control. The renewal probe independently
  cell-audits 46,124 small component instances and exactly classifies both
  h=420 hostiles. Normal/optimized/hash-seeded replays match frozen outputs.
---

# THM-4335 -- owner permutations, component addresses, and minority renewal

**PROVED ELEMENTARY OWNER/COMPONENT/RENEWAL LEMMAS + PROVED RELATIVE TO
CITED LOWER-DIMENSIONAL LRC FOR ONE COROLLARY + VERIFIED-EXACT CONTROLS.
LRC(14) REMAINS OPEN.**

## 1. Unit-tail owner theorem

For a finite positive integer set `A`, put

```text
G_A={theta in R/Z:min_(a in A)||a theta||>=1/14}.      (1)
```

Let `2<=d<=7`, let `C` be a nonempty finite set of positive integers, and
let `T` be a set of exactly `d` distinct positive integers satisfying
`gcd(w,d)=1` for every `w in T`. Put

```text
S=dC union T.                                          (2)
```

For a chosen real representative of `theta`, define the strict bad-sheet
mask

```text
B_w(theta)={j in Z/dZ:||w(theta+j)/d||<1/14}.          (3)
```

When `||w theta||<d/14`, let `n_w(theta)` be the unique nearest integer to
`w theta`. Then

```text
B_w(theta)=empty                         if ||w theta||>=d/14,
B_w(theta)={kappa_w(theta)}              if ||w theta||< d/14,
kappa_w(theta)=-n_w(theta)w^(-1) mod d.                (4)
```

Equality in the first line gives an empty **strict** mask. At `d=7`, an
active phase still has distance strictly below `1/2`, so its nearest integer
is unique; a half-integer is an equality-safe inactive boundary and receives
no owner.

Changing the representative from `theta` to `theta+1` sends

```text
n_w -> n_w+w,                    kappa_w -> kappa_w-1. (5)
```

Individual sheet labels are gauge-dependent, while coverage and bijectivity
are invariant.

> **Owner-permutation characterization.** The following are equivalent:
>
> 1. `M(S)<1/14`, equivalently `G_S` is empty;
> 2. for every `theta in G_C`, all tails are active and
>
>    ```text
>    w -> kappa_w(theta)
>    ```
>
>    is a bijection from `T` to `Z/dZ`.

Thus equality in the scalar sheet-capacity bound becomes an exact labelled
permutation condition. Activity without bijectivity is not enough: owner
collisions leave a free sheet.

## 2. Constant component addresses and pair separation

Assume the equivalent strict-failure conditions above. Let `I` be a
positive-length connected component of `G_C`, with its unique nonwrapping
closed lift

```text
I=[L,R] subset (0,1),                   W=R-L>0.        (6)
```

For each `w in T`, there is one integer `n_w`, constant on `I`, such that

```text
I subset (n_w/w-d/(14w), n_w/w+d/(14w)).              (7)
```

The owner labels `-n_w w^(-1)` form one fixed permutation on `I`. For
distinct `u,v in T`, define

```text
Delta_(u,v)=n_u v-n_v u,       D_(u,v)=|Delta_(u,v)|. (8)
```

Then

```text
d does not divide Delta_(u,v),
gcd(u,v) divides Delta_(u,v),                              (9)

W < min(
 d/(7u),
 d/(7v),
 d/(14u)+d/(14v)-D_(u,v)/(uv)
).                                                       (10)
```

The last expression is positive. The minimum in `(10)` is the exact length
of the intersection of the two addressed open teeth in `(7)`, not a union
bound.

## 3. Endpoint-denominator component escape

Write the endpoints in lowest terms as

```text
L=A_L/Q_L,                     R=A_R/Q_R.              (11)
```

Every positive component endpoint is a wall of a member of `C`, with
unreduced numerator `14k+-1`. Therefore

```text
14|Q_L,                         14|Q_R.                (12)
```

For every `w in T`, strict failure forces

```text
wW <  d/7-1/Q_L,
wW <  d/7-1/Q_R,
wW <= d/7-1/Q_L-1/Q_R.                                (13)
```

Consequently the row `(2)` is safe whenever there are a positive component
`I` and a tail `w in T` satisfying either

```text
wW >= d/7-1/min(Q_L,Q_R),                             (14)
```

or

```text
wW > d/7-1/Q_L-1/Q_R.                                (15)
```

The inequality in `(14)` is nonstrict, while `(15)` is strict. This is the
`d`-sheet component-address extension of THM-4331's two-sheet theorem.

The strictness in `(15)` is sharp for this proof operation. At `d=3`, take

```text
C={9},          I=[19/42,23/42],          w=4.         (16)
```

The component `I` lies strictly inside the open effective tooth

```text
(25/56,31/56),                                      (17)
```

with both endpoint gaps equal to `1/168`, yet

```text
wW=8/21=3/7-1/42-1/42.                               (18)
```

This is a geometric equality hostile, not an unsafe full row.

## 4. Apex-intercept obstruction

Put

```text
R_C=max C,       m_C=M(C)>1/14,       rho=(m_C-1/14)/R_C, (19)
```

and choose a core maximizer `theta_0`. Strict failure necessarily gives,
for every `w in T`,

```text
||w theta_0||+w rho<d/14.                             (20)
```

For every pair `u!=v`, with the tooth integers on the radius-`rho`
core-safe interval,

```text
D_(u,v)/(uv)+2rho<d/(14u)+d/(14v).                   (21)
```

Either reversed nonstrict inequality is therefore a safety certificate.
When `|C|=13-d`, cited lower-dimensional LRC gives

```text
m_C>=1/(14-d),              rho>=d/[14(14-d)R_C].     (22)
```

Hence every strict counterexample on this equality branch obeys

```text
w<(14-d)R_C                                             (23)
```

for every tail and

```text
2uv<(14-d)R_C(u+v-14D_(u,v)/d)                        (24)
```

for every pair. Equation `(23)` recovers THM-4330's scale-free boundary;
`(20)` retains the missing intercept, and `(24)` retains the nonzero
owner-address determinant.

For `d=5`, at any phase where all five owners exist, they cover all sheets
if and only if

```text
product_(w in T)(X-kappa_w)=X^5-X       in F_5[X],     (25)
```

equivalently

```text
(sum kappa, sum kappa^2, sum kappa^3, sum kappa^4)
 =(0,0,0,-1)                             in F_5^4.     (26)
```

## 5. Proof of Sections 1--4

If `j in B_w(theta)`, there is an integer `m` with

```text
|w theta+wj-dm|<d/14<=1/2.                            (27)
```

Thus `n=dm-wj` is the unique nearest integer to `w theta` and
`n+wj=0 mod d`. The unit hypothesis gives the unique solution in `(4)`.
The converse follows by choosing that sheet and dividing
`|w theta-n|<d/14` by `d`. This proves `(4)`, including the strict equality
boundary. Equation `(5)` is immediate.

For `theta in G_C`, all `d` lifts preserve the core:

```text
||dc(theta+j)/d||=||c theta||.                         (28)
```

Exactly `d` singleton masks cover the `d` sheets if and only if their owners
are distinct. Conversely, every physical phase maps under multiplication by
`d` either outside `G_C`, where the core is strictly dangerous, or into this
permutation cover. This proves both directions of the owner characterization.

Under strict failure, a connected closed component `I` lies in one connected
component of every tail's open active locus. Those components are exactly
the teeth `(7)`, so `n_w` and the owner are constant. Moreover

```text
Delta_(u,v)=uv(kappa_v-kappa_u) mod d.                 (29)
```

Distinct owners prove the first part of `(9)`; the second is immediate.
The tooth centres are separated by `D/(uv)` and have radii
`d/(14u),d/(14v)`. Direct interval intersection gives the minimum in `(10)`.
The closed `I` is strictly inside that open intersection, proving `(10)`.

For `(13)`, let the effective `w`-tooth containing `I` be

```text
(a,b)=((14k-d)/(14w),(14k+d)/(14w)).                  (30)
```

Strict containment and `(12)` give quantized positive gaps

```text
L-a=[wA_L-(Q_L/14)(14k-d)]/(wQ_L)>=1/(wQ_L),
b-R>=1/(wQ_R).                                        (31)
```

The tooth width is `d/(7w)`. Keeping one gap and the strict positivity of
the other gives the two strict bounds in `(13)`; keeping both quantized gaps
gives its nonstrict additive bound. The contrapositives are `(14)--(15)`.
At `d=7`, adjacent teeth meet only at excluded half-integer boundaries, so
the same connected-component argument applies. A wrapping tooth is handled
by its compatible real lift.

The core envelope is `R_C`-Lipschitz. Thus the closed circular interval of
radius `rho` around `theta_0` lies in `G_C`; strict failure puts it strictly
inside every tail tooth. Comparing its centre and radius gives `(20)`, and
comparing its length with a pair-tooth intersection gives `(21)`. Substituting
`(22)` proves `(23)--(24)`.

Finally, a permutation of `F_5` plainly gives `(25)--(26)`. Conversely,
Newton identities in degrees one through four turn `(26)` into
`e_1=e_2=e_3=0,e_4=-1`. The owner polynomial is `X^5-X-e_5`; evaluating it
at any one of its five listed roots and using `a^5=a` gives `e_5=0`. This
proves the converse. **QED.**

## 6. Initial-segment and minority-wall controls

For `m=13-d`, the first positive component of

```text
C_m={1,...,m}
```

is

```text
I_m=[1/14,13/(14m)],        W=d/(14m),
(Q_L,Q_R)=(14,14m).                                  (32)
```

Indeed, throughout this interval every speed `c=1,...,m` stays in its
closed safe band `[1/14,13/14]`; immediately to the left speed `1` fails,
and immediately to the right speed `m` fails. Thus `(32)` is the whole
connected component, including its two safe boundary points.

For one selected `d`-unit tail, the least admissible thresholds from
`(14)--(15)`, compared with the width-only threshold `wW>=d/7`, are

| `d` | 2 | 3 | 4 | 5 | 6 | 7 |
|---:|---:|---:|---:|---:|---:|---:|
| width-only selected-tail threshold | 23 | 20 | 19 | 16 | 17 | 12 |
| addressed selected-tail threshold | 17 | 17 | 17 | 16 | 13 | 12 |

Thus the endpoint address numerically improves the selected-tail threshold
for `d=2,3,4,6`. The `d=6` entry is not a new row-level gain: six distinct
positive `6`-units already force `max T>=17`, so some tail always triggers
the raw-width certificate; the exact control at selected tail `13` also
contains tail `17`. The exact probe gives a literal safe lift at every
displayed addressed threshold, while isolating a genuine endpoint gain for
`d=2,3,4`.

There is also a `p=5` minority-wall component control. Take

```text
C=(1,3,5,7,9,11,13,168),
T=(1,3,7,9,141).                                      (33)
```

The physical row `5C union T` has thirteen distinct speeds and is of the
minority form `{840} union W` with twelve positive odd speeds. It satisfies
all THM-366 denominator gates. The core has `58` positive components,
including

```text
I=[281/2352,293/2352],       W=1/196,
(Q_L,Q_R)=(2352,2352).                                (34)
```

For this component the raw, one-endpoint, and additive selected-tail ratios
are

```text
140,                 1679/12,                 839/6.  (35)
```

The least admissible five-unit is `141`, far below the old scale-free value
`9 max(C)=1512`. At quotient phase `81/658`, tails `3` and `141` are
inactive, and sheet `3` gives the physical witness

```text
t=411/658,                      clearance=1/14.        (36)
```

### 6.1 A marginal hostile closed by owner disappearance

A sharper `h=420`, `r_5=5` control is

```text
C=(13,168,349,375,711,737,1073,1099),
T=(1,21,23,327,689),
S=5C union T.                                         (37)
```

It is primitive, contains the anchor `840`, has twelve other positive odd
speeds, and is complete for every THM-366 denominator. It represents every
needed sign class in each useful doubled-denominator unit bank, and every
represented adaptive capacity is at least one, with equality uniquely at
`d=5`. Both THM-4330 half-turn clocks are killed, and the old height-only
condition does not fire.

The core maximum is

```text
M(C)=90/181,                       theta_0=7/181.      (38)
```

The pair `{13,168}` gives the upper bound: if both distances exceeded
`90/181=1/2-1/362`, writing their phases as a half-integer plus errors of
absolute value below `1/362` would make a half-integer cancel an error of
absolute value below `181/362=1/2`, impossible. Every core label is
congruent to `13` or `168` modulo `181`, so `theta_0` attains the bound.

At `theta_0`, the tail owners are exactly

```text
(0,4,3,1,2),                                           (39)
```

a perfect permutation. But tails `327` and `689` violate `(20)`. Moving
within the guaranteed core-safe radius to

```text
theta_1=59/1526
```

makes both masks empty. The remaining owner masks are `{0},{4},{3}`, so
sheet `2` gives

```text
t=(theta_1+2)/5=3111/7630,
min_(s in S)||st||=551/7630=1/14+3/3815.              (40)
```

This is a strict new certificate for one safe row, not a counterexample or
a family-level closure.

## 7. Deterministic minority-anchor renewal theorem

There is a complementary physical-component formulation that does not assume
capacity equality. Let

```text
S_h={2h} union W,
```

where `h>=1` and `W` is any finite set of positive odd integers. The
`2h` anchor-safe components are

```text
I_k=[L_k,R_k]
   =[(14k+1)/(28h),(14k+13)/(28h)],     0<=k<2h,       (41)
```

each of width `3/(7h)`. Write the open danger tooth of tail `w` with
integer address `n` as

```text
T(w,n)=((14n-1)/(14w),(14n+1)/(14w)).                 (42)
```

> **Shortest-cover-chain trichotomy.** For every component address `k`,
> exactly one of the following occurs.
>
> 1. Some point of `I_k` belongs to no tail tooth. It is a `1/14` witness
>    for `S_h`.
> 2. One tooth `T(w,n)` strictly contains `I_k`. Necessarily
>
>    ```text
>    3w<h.                                             (43)
>    ```
>
> 3. A farthest-reach greedy cover gives a minimum-cardinality chain
>
>    ```text
>    T(w_1,n_1),...,T(w_m,n_m),             m>=2,      (44)
>    ```
>
>    with strictly increasing left and right walls and proper consecutive
>    crossings. If
>
>    ```text
>    D_i=n_(i+1)w_i-n_iw_(i+1),
>    q_i=w_i+w_(i+1)-14D_i,                            (45)
>    ```
>
>    then
>
>    ```text
>    D_i>0,
>    q_i in 2gcd(w_i,w_(i+1)) Z_(>0),                 (46)
>
>    |T(w_i,n_i) intersect T(w_(i+1),n_(i+1))|
>      =q_i/(14w_iw_(i+1))
>      >=1/[7lcm(w_i,w_(i+1))].                       (47)
>    ```

The complete chain has an exact budget. Define

```text
A=w_1(14k+1)-2h(14n_1-1)>0,
B=2h(14n_m+1)-w_m(14k+13)>0.                          (48)
```

Both `A` and `B` are positive odd integers, and

```text
3/(7h)
 =sum_(i=1)^m 1/(7w_i)
  -A/(28hw_1)-B/(28hw_m)
  -sum_(i=1)^(m-1)q_i/(14w_iw_(i+1)).                 (49)
```

Repeated speed labels may occur in `(44)`: the states are addressed teeth
`(w,n)`, not merely runners. Every successful renewal pays both endpoint
tolls and every gcd-quantized collision toll.

### Proof

If the open teeth do not cover `I_k`, compact one-dimensional geometry gives
case 1. Otherwise choose, among all teeth containing `L_k`, one with farthest
right wall. If it does not pass `R_k`, choose among all teeth strictly
containing the current right frontier one with farthest right wall, and
iterate. Coverage supplies each choice, and only finitely many teeth meet
`I_k`, so the construction terminates.

The greedy rule excludes containment; consecutive left and right walls
increase strictly and cross properly. A tooth two steps later cannot overlap
the tooth two steps earlier, since it would have been eligible and extended
farther at the preceding choice. Thus only adjacent chain teeth overlap.
The usual frontier exchange proves minimum cardinality: after `r` choices no
`r`-tooth chain starting at `L_k` can reach beyond the greedy frontier,
because an interval doing so would have been an eligible farther-reaching
choice at the first index of disagreement.

In the one-tooth case, strict containment makes `1/(7w)>3/(7h)`, proving
`(43)`. For a proper crossing, increasing tooth centres give `D_i>0`, and
direct subtraction in the order “exiting right wall minus entering left
wall” gives

```text
right(T(w_i,n_i))-left(T(w_(i+1),n_(i+1)))
 =q_i/(14w_iw_(i+1))>0.                               (50)
```

Let `g_i=gcd(w_i,w_(i+1))`. Both `D_i` and the pair sum are divisible by
`g_i`; after division by the odd `g_i`, the pair sum is even and `14D_i/g_i`
is even. Thus `q_i` is a positive multiple of `2g_i`, proving `(46)--(47)`.

At the left endpoint, `A` is odd minus even; at the right, `B` is even minus
odd. Strict containment gives the positive signs in `(48)`, hence `A,B>=1`.
The union of the chain is one interval and only adjacent overlaps have
positive length. Subtracting its two endpoint overhangs and all consecutive
overlaps from the sum of tooth widths proves `(49)`. **QED.**

## 8. Owner-bit and stochastic-process connections

Write

```text
k=j+h epsilon,          0<=j<h,          epsilon in {0,1}. (51)
```

Then `I_k=(J_j+epsilon)/2`, where `J_j` is the `j`th component of
THM-4330's quotient set `G_h`. If `T(w,n)` is active on this physical
component, its quotient nearest integer is

```text
N_w=2n-epsilon w,                    N_w mod 2=epsilon. (52)
```

Thus every state of a renewal chain has the same MC2 owner bit. At the right
wall of an exiting tooth, that owner is equality-safe and the incoming tooth
is strictly active. The transition is exactly a collision-backed replacement,
not an independent random renewal.

The link to THM-592 is also exact but limited. Each edge in `(45)` has

```text
rho_i=D_i/(w_i+w_(i+1))<1/14,
q_i=14(w_i+w_(i+1))(1/14-rho_i)>0.                    (53)
```

Interior overlaps need not be exposed kinks of the global lonely-measure
function. The addressed renewal graph retains them because they preserve the
pointwise cover predicate. The faithful state is
`(j,epsilon;w,n,q_i)`, not a Bernoulli owner bit.

There is nevertheless a useful exact occupation identity. Let

```text
D_v={t:||vt||<1/14},             E=G_{2h},
m(t)=#{w in W:t in D_w},
F=mu({t in E:m(t)=0}),
Omega=int_E(m(t)-1)_+ dt,
A_0=sum_(w in W)mu(D_(2h) intersect D_w).             (54)
```

Since `mu(E)=6/7` and each tail-danger set has mass `1/7`, integrating the
pointwise identity `m=1_(m>=1)+(m-1)_+` gives

```text
F=(6-|W|)/7+A_0+Omega.                                (55)
```

For the twelve-tail minority branch,

```text
F=A_0+Omega-6/7.                                      (56)
```

This is deterministic defect bookkeeping, not a probability or independence
claim. A strict counterexample would have `F=0`; the identity alone does not
force `F>0`.

If `P_h(W)` denotes the component addresses spanned by a single tooth, then
a strict counterexample must support collision-backed renewal on at least

```text
2h-|P_h(W)|                                            (57)
```

addresses. The unresolved aggregation problem is reuse: one tail pair may
pay transition obligations at many different component addresses.

## 9. Exact minority-renewal controls

The renewal probe independently agrees with literal cell reachability on
`46,124` small component instances. It finds a sharp overlap-quantization
control

```text
h=10,               W={3,13},               k=6,
chain=((13,4),(3,1)),                     q_1=2.       (58)
```

For each THM-4330 joint hostile at `h=420`, with `P=1287` or `P=9009`, the
exact `840`-component classification is

```text
missing tail cover       726
one-tooth span           110   (all owned by w=11)
collision renewal          4
maximum chain length       3
minimum q_i                 4.                         (59)
```

The occupation masses differ and remain separately frozen in the output;
only the component counts coincide. Both rows return the same first exact
witness

```text
t=939/141274,
clearance=1/14,
unique binding speed=10091.                            (60)
```

This supplies a new exact certificate format for two already-known safe
rows. It is not extrapolated to the full `420|h` wall.

## 10. Connection contract and scope

```text
source:       equality-capacity rows dC union T, and minority rows
              {2h} union W with W odd
target:       sheet-owner permutations; addressed tail-tooth intersections;
              directed component-cover chains
map:          theta=dt; n_w=round(wtheta); kappa=-n_ww^(-1); or lift each
              quotient component to I_(j,epsilon) and orient tooth crossings
preserved:    exact strict-cover/witness predicate, equality, sheet and
              component addresses, endpoint and overlap slack
destroyed:    owner-only quotients forget distance to the next endpoint;
              greedy chains forget alternative covers and cross-component
              reuse; component width alone forgets endpoint denominators
sidecar:      Q_L,Q_R; n_w; Delta; core apex and rho; physical (j,epsilon);
              q_i, gcd/lcm, and cyclic reset order
decisive test: one failed owner permutation or one failed directed cover
              returns an exact physical safe lift
```

The unit-owner theorem requires exactly `d` distinct unit tails and
`2<=d<=7`. With more tails, full cover is not a bijection; without units or
above seven sheets, one tail may kill more than one sheet. The scale-free
corollary alone uses lower-dimensional LRC. The component criteria are
sufficient, not necessary, and their failure is not danger. The renewal
theorem is a local obstruction compiler; no theorem here forces a missing
component address in every row. No arbitrary `r_5>5` branch, affine entry
theorem, minority-branch closure, or proof of LRC(14) follows.

## 11. Reproduction

From the repository root:

```bash
python3 -B 04-computation/lrc14_unit_component_address_escape_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_unit_component_address_escape_probe_20260902.out -
python3 -B 04-computation/lrc14_unit_component_address_escape_probe_20260902_independent_audit.py \
  | diff -u 05-knowledge/results/lrc14_unit_component_address_escape_probe_20260902_independent_audit.out -
python3 -B 04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out -
python3 -B 04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py \
  | diff -u 05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out -
```
