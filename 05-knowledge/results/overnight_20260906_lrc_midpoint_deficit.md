# Midpoint curvature payments and the affine deficit they cannot see

**Status: PROVED ANALYTICALLY; FINITE-EXACT named controls; INDEPENDENTLY
AUDITED.** Compulsory midpoint progressions give an exact, optimally paired
lower bound on the THM-4422 boundary deficits. They do not recover those
deficits: a complete two-line carrier family has arbitrarily many compulsory
progressions and identically zero midpoint curvature, while its first
boundary deficit is positive. The resulting sufficient payment condition
is not necessary; whether it holds for every remaining count-dense full-rank
dictionary is **OPEN**. No universal projection inequality, physical entry, or LRC(14)
claim follows.

Root independently checked completeness of the two-line family, every
residue and strict-roof boundary, the hinge identity, the convex matching
exchange proof (including unmatched endpoints), and the count-safe scope.
The normal and optimized verifier outputs were replayed and both equal
the frozen transcript. The audit passed without a mathematical correction.

## Inheritance and the bounded question

The closest proved mechanism is
[THM-4422 / projection-deficit-and-beatty-row-reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md),
which retains the full support and asks whether one weighted boundary
deficit pays the carrier-count surplus. The new input is the audited
[full-cap carrier theorem](overnight_20260906_lrc_cap_carriers.md): two
carriers of the same owner color and lattice parity force their midpoint
to be live. The canonical count-dense `A_2` hostile is `(19,23,29)`; the
corrected near miss is the speed-independent convex average refuted in
THM-4422. The least-used sidecar here is the **affine baseline of a deficit
along an owner-colored line**, which survives when its second differences
vanish. THM-4423 is not used.

This is a bounded continuation of the exact incidence transfer from
Guy--Kelly recorded in the cap report. It makes no further literature claim.
Before deriving the family below, targeted searches for `2m+3`, `m-3k`,
`5m-9`, `3m-1`, and Jensen/deficit language in the recent LRC theorem and
results surfaces found no existing statement of this mechanism. The earlier
norm-four one-ray families are inherited positive controls, not new results.

The concrete proposal was: **can compulsory AP curvature pay the exact
boundary invoice?** The concept board has the following final state.

| Concept | Decisive observation | Surviving information |
|---|---|---|
| Full convex dictionary | Same color and parity force the midpoint | Every payment uses actual live carriers |
| AP incidence | Many progressions can have zero curvature | Counts need a boundary-location sidecar |
| Deficit hinges | Curvature is supported at two explicit kinks | Exact weighted payment is possible |
| Pairwise budget | Reusing endpoints can overspend their deficit | Disjoint endpoint pairs give a rigorous bound |
| Affine carrier lines | Positive deficit can be affine everywhere | Baseline or anchored values must be retained |
| Projection selector | Different coordinates may use different pairings | One successful coordinate suffices |

## 1. Exact local payment

Use the complete dictionary and filters of THM-4422:

```text
w=(a,b,c), 1<=a<b<c, gcd(a,b,c)=1,
a,b,c odd and nonzero modulo 3,
L_w={C in Z^3 : C dot w=0},
p_i(C)=3(w_j+w_k)-14|C_i|,
Omega_w={C in L_w : all p_i(C)>0, all C_i!=0 mod 3}.
```

Write `N=|Omega_w|`, `q=3/(7c)` and

```text
d_i(C)=[1-c p_i(C)/(6w_jw_k)]_+,
D_i=sum_(C in Omega_w) d_i(C),
I=N-2c/11.
```

The inherited exact target is `max_i D_i>=I`, equivalent to
`min_i S_i<=6/77`. Put

```text
alpha_i=7c/(3w_jw_k),
T_i=3(w_j+w_k)/14-3w_jw_k/(7c)>0.
```

Then `d_i(C)=alpha_i (|C_i|-T_i)_+`. Positivity of the thresholds follows
from `T_1=3(c-b)/14`, `T_2=3(c-a)/14`, and
`c(a+b)>2ab` for the third one.

Partition the dictionary into its two owner colors
`tau(C)=aC_1=bC_2=cC_3 mod 3`, and its four classes in `L_w/2L_w`.
The lattice has rank two; because all three speeds are odd, these parity
classes can also be represented by `000,011,101,110` in ambient coordinates.
For two distinct members `C=M-H`, `E=M+H` of one bucket, the midpoint `M`
is integral, remains in the kernel, satisfies all three strict roofs by
convexity, and retains the nonzero owner color. It belongs to the **full**
dictionary.

The exact local curvature payment is

```text
kappa_i(C,E)=d_i(C)+d_i(E)-2d_i(M)
 =alpha_i { [|H_i|-|M_i-T_i|]_+
           +[|H_i|-|M_i+T_i|]_+ } >=0.                 (1)
```

To prove (1), expand `(|x|-T)_+=(x-T)_++(-x-T)_+` and take the second
difference of each hinge. This also identifies the equality boundary:
curvature vanishes exactly when the open coordinate interval between the
endpoints crosses neither kink `T_i` nor `-T_i`. Touching a kink only at an
endpoint does not pay. In particular, AP existence alone gives no positive
payment.

For each coordinate `i`, sort each bucket by `C_i`, and pair its smallest
coordinate with its largest, then its second smallest with its second
largest, and so on. Leave the median unmatched if the bucket has odd size.
Different coordinates are allowed different pairings. Let `P_i` be these
endpoint-disjoint pairs and `U_i` the unmatched carriers. Define

```text
B_i=sum_((C,E) in P_i) kappa_i(C,E),
R_i=2 sum_((C,E) in P_i) d_i((C+E)/2)
      +sum_(C in U_i) d_i(C).
```

Every original endpoint is used exactly once or is unmatched, giving

```text
D_i=B_i+R_i,       B_i>=0, R_i>=0.                    (2)
max_i B_i>=N-2c/11  ==>  min_i S_i<=6/77.              (3)
```

Repeated midpoints are counted with their pair multiplicity in `R_i`;
this does not affect the identity. Arbitrarily summing over all APs has no
such endpoint budget and is not justified by (2).

## 2. Extreme pairing is optimal for this payment model

For any convex real function `f`, let

```text
W(x,y)=f(x)+f(y)-2f((x+y)/2).
```

Among all matchings of a finite sorted multiset of real coordinates,
pairing the extremes successively maximizes the sum of `W`. For odd
cardinality it leaves a median unmatched. This includes ties and matchings
that initially leave more than one point unmatched.

Indeed `W>=0`, so a maximum can be extended to cardinality `floor(n/2)`.
For `a<=b<=c<=d`, the two midpoint coordinates of the outer pairing
`(a,d),(b,c)` are majorized by those of `(a,c),(b,d)`, which are majorized
by those of `(a,b),(c,d)`. All have the same sum. Convexity therefore gives

```text
W(a,d)+W(b,c) >= W(a,c)+W(b,d)
              >= W(a,b)+W(c,d).                      (4)
```

If both extremes are paired elsewhere, (4) exchanges those pairs for the
extreme pair without reducing payment. If one extreme is unmatched,
replace the pair containing the other extreme by the full extreme pair:
`W(x,y)` is nondecreasing when the interval `[x,y]` is enlarged. For
completeness, this last fact follows from convex secant slopes, or first
for a hinge from
`W(x,y)=min(t-x,y-t)_+` for `f(z)=(z-t)_+`, then from the positive hinge
representation of a convex function on the finite interval in question.
Remove the chosen extreme pair and induct.

Thus `B_i` is the **largest endpoint-disjoint same-bucket curvature
payment**, not just one convenient matching. This says nothing about
fractional matching relaxations or a model that also budgets baseline
values. Its failure can identify missing information even after the
endpoint matching has been optimized.

## 3. A complete family invisible to every compulsory curvature

**Proposition.** For every odd integer `m>=5` with `3` not dividing `m`,
set `w=(1,m,2m+3)` and `c=2m+3`. Its complete carrier dictionary is

```text
Omega_w = +- { (m-3k,-2k-1,k) :
               k=2 mod 3,
               (5m-9)/42 < k < (3m-1)/14 }.            (5)
```

On every positive branch carrier in (5),

```text
d_1=(11m-9-42k)/(6m)>0,     d_2=d_3=0.                (6)
```

Consequently every same-owner AP, and hence every compulsory same-owner
same-parity AP, has `kappa_i=0` in all three coordinates. Nevertheless the
number of carriers and compulsory APs is unbounded. When the support is
nonempty, `D_1>0=B_1=B_2=B_3`.

**Completeness proof.** For any live carrier put `ell=C_2+2C_3`. The kernel
equation gives `C_1=-m ell-3C_3`, and the strict roofs imply

```text
|m ell| <= |C_1|+3|C_3| < 9(m+1)/7 < 2m.
```

Thus `ell` is `-1,0,1`. The middle value forces `C_1=0 mod 3` and is
forbidden. Negation reduces the two remaining cases to `ell=-1`. Write
`k=C_3`, giving the displayed vector in (5). The first roof implies
`k>(5m-9)/42>0`; the third gives `k<3(m+1)/14` and hence `C_1>0`.
The second roof then sharpens the upper bound to `k<(3m-1)/14`.
The two nonzero residue conditions for `C_2,C_3` force `k=2 mod 3`.
Conversely these bounds and this residue ensure all three strict roofs and
all three nonzero residues, proving equality with the full dictionary.

For the deficit formula, the first expression in (6) follows directly
from THM-4422. Its numerator is greater than `2m-6>0` by the upper bound
on `k`. A positive second deficit would require

```text
3m-4 < 14k < 3m-1.
```

The integer `14k-3m` would have to be `-3` or `-2`. Odd `m` excludes
`-2`, and `k=2 mod 3` makes this difference `1 mod 3`, excluding `-3`.
Finally

```text
T_3=3m/14+9/[14(2m+3)] > 3m/14 > k,
```

so the third deficit vanishes. All carriers on the positive branch have
owner color `m mod 3`; those on the negative branch have the other color.
A same-owner AP therefore stays on one affine line, where each of the
three deficits is affine. Its second difference is zero.

The open interval for `k` has length `c/21`, and one residue class out of
three is used. Hence

```text
N=2c/63+O(1).
```

The parity-bucket AP lower bound from the cap report proves unbounded
compulsory incidence as `N` grows. This family is always **count-safe**:
for `c>=17`, its positive branch count is less than `c/63+1<=c/11`;
the remaining case `m=5,c=13` has empty support. Therefore it is not a
counterexample to (3) in the hard regime `I>0`.

## 4. Named positive and hostile controls

The verifier rebuilds carriers from the raw integer kernel and strict
roofs. It does not import the family formula into this construction.
Its controls are fixed named rows, not a new height census.

| Triple | N | D | B | Invoice I | Interpretation |
|---|---:|---|---|---|---|
| `(19,23,29)` | 6 | `(86/23,50/19,1342/437)` | `(0,0,0)` | `8/11` | Canonical `A_2`: criterion is not necessary |
| `(23,29,37)` | 6 | `(82/29,70/23,1410/667)` | `(0,0,0)` | `-8/11` | Sharp cap control; count already pays |
| `(29,35,41)` | 10 | `(12/7,120/29,3838/3045)` | `D` | `28/11` | Genuine positive-curvature payment; inherited norm-four ray |
| `(49,59,61)` | 10 | `(1006/177,646/147,43420/8673)` | `(0,52/21,2684/1239)` | `-12/11` | Full-rank control; substantial affine residual |
| `(1,137,277)` | 10 | `(2660/411,0,0)` | `(0,0,0)` | `-444/11` | Complete two-line flat hostile |
| `(1,499,1001)` | 30 | `(9970/499,0,0)` | `(0,0,0)` | `-152` | Larger flat hostile |

At `(49,59,61)`, the compulsory triples are

```text
(8,13,-19), (11,-5,-4), (14,-23,11)
```

and its negative. Each contributes curvature
`(0,26/21,1342/1239)`. The first deficit is invisible to both progressions,
and the residual in (2) is
`(1006/177,94/49,24632/8673)`.

At `(1,137,277)`, the positive branch has `k=(17,20,23,26,29)` and vectors

```text
(86,-35,17), (77,-41,20), (68,-47,23),
(59,-53,26), (50,-59,29).
```

Their negatives complete the dictionary. This is a full-rank union of two
parallel affine lines, not a single relation ray. It is a named ten-carrier
witness, with no minimal-height claim.

The first failed implication in the proposed unweighted route is
**compulsory AP implies positive boundary curvature**. Formula (1) pinpoints
the failure: no kink lies inside the endpoint interval. The stronger
attempt to reconstruct all boundary deficit from such curvature also fails
on (5), where the entire positive deficit is an affine baseline. The
strongest survivor proved here is the exact optimal matching lower bound
(2)--(3).

## 5. Connection contract and next starting point

| Required field | Exact content |
|---|---|
| Source | Full owner-colored carrier lattice with compulsory midpoint APs |
| Target | Maximum normalized boundary deficit of THM-4422 |
| Map | Assign each same-bucket endpoint pair its three hinge second differences |
| Preserved predicate | Nonnegative endpoint-disjoint payments are bounded by the actual deficit |
| Destroyed information | Affine restrictions of each deficit to the owner-colored AP components |
| Needed sidecar | Anchored deficit values or affine baselines, plus midpoint multiplicities |
| Cheapest decisive test | Apply any proposed baseline recovery to the complete family (5) and the canonical `A_2` row |

The next precise question is whether a **baseline-aware** certificate can
pay `I>0` after the cap class and inherited ray classes have been removed.
One safe form is to prove a lower bound on the explicit nonnegative
residual `R_i`, using anchored boundary values or the exact Beatty-row
endpoints from THM-4422, and add it to `B_i`. Curvature-only reasoning has a
known kernel and must not silently discard this residual. No claim that
`B_i` alone suffices for every remaining count-dense full-rank dictionary
has been established or refuted here.

## Reproduction and scope

Run from the repository root:

```text
python 04-computation/overnight_20260906_lrc_midpoint_deficit.py
python -O 04-computation/overnight_20260906_lrc_midpoint_deficit.py
```

The standalone standard-library verifier uses exact integers and rational
numbers, with explicit non-optimizable gates. It checks the six named
carrier rows above, all eligible AP hinge identities in those rows,
the complete family formula at
`m=(5,7,11,13,17,31,43,61,79,101,137,199,499,1001)`, and extreme-pairing
optimality against an independent exhaustive matching recursion on 483
scalar cases, including coordinate ties. Family samples are finite controls
for the independent all-height proof, not its basis. The output records
5,220 checks. Normal and optimized output byte-match.

```text
script SHA-256 (LF bytes): 95c9fc35d6ca6ff5905d0550365f3abed1b8e8c71e62ccb7aa97e511c1f2fbbc
output SHA-256: c25fae5edeba53d3c303a32f3f5120d66716b9b2e7df42e5a04641efb4541f36
```

Companions:
[verifier](../../04-computation/overnight_20260906_lrc_midpoint_deficit.py),
[output](overnight_20260906_lrc_midpoint_deficit.out).
The finite scalar controls test the actual hinge functions used here;
the general convex matching statement is proved analytically in section 2.

The report makes no new physical-network, universal projection, or LRC(14)
claim. It closes one proposed incidence-to-payment mechanism with an exact
survivor, a complete infinite hostile family, and the missing affine
coordinate.
