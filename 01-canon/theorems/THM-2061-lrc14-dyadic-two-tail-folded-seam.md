---
id: THM-2061
title: "LRC(14) dyadic two-tail folded seam, divisor pins, and metric box"
status: >
  PROVED REDUCTION; NOT LRC(14). A strict counterexample in the sole
  imprimitive eleven-core seam 2C union {x,y} is equivalent to containment of
  the closed weak-safe set of C in one strict folded diamond at height 6/7.
  Nearest-integer ownership is constant on each component. Every surviving
  core is divisor-complete through 14; the odd tails satisfy x,y<12 max(C)
  and a gcd-sharpened centre-separation inequality; and the folded target has
  sharp measure at most 4/63. An exact bounded census closes every normalized
  quotient core contained in {1,...,19}, uniformly over the previously
  unbounded odd tails. The folding identity itself is THM-615/774; this
  theorem is its strict 1/14 specialization and atlas reduction.
source: codex-2026-07-21-LRC14-dyadic-seam
script: 04-computation/lrc14_dyadic_two_tail_folded_seam_codex_20260721.py
result: 05-knowledge/results/lrc14_dyadic_two_tail_folded_seam_codex_20260721.out
script_sha256: af1ad8ef16830a09ac030551774958894b0c9c36415dc2d0c2d1e6c8042793c9
result_sha256: 22fb09fa81f67418d8deaf62a5a330bc7aff3928189a81e5b3586b40203370da
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2060
  - THM-774
  - LRCUpTo13
related:
  - THM-615
  - THM-761
  - THM-769
  - THM-775
  - THM-2057
  - THM-2064
  - HYP-8846
---

# THM-2061 -- The LRC(14) dyadic two-tail seam

Put `delta=1/14`. Let `C` be a finite nonempty set of positive integers and
let `x>y` be distinct positive odd integers. Define

```text
phi_C(theta)=min_(c in C)||c theta||,
G_C={theta in R/Z:phi_C(theta)>=delta},
S=2C union {x,y}.                                         (1)
```

For the thirteen-speed application, `|C|=11` and `gcd(C)=1`: THM-2060 says
that every primitive row having an imprimitive eleven-speed deletion is safe
unless it has exactly this normalized form.

## 1. Exact strict folded obstruction

Put

```text
a=(x+y)/2,             b=(x-y)/2,
Q_(a,b)(theta)=||a theta||+||b theta||.                    (2)
```

Then

```text
M(S)<1/14
  iff G_C subset H^o_(a,b),
H^o_(a,b)={theta:Q_(a,b)(theta)>6/7}.                      (3)
```

Equivalently, for every `theta in G_C`,

```text
||x theta||<1/7,       ||y theta||<1/7,                    (4)
```

and the unique nearest integers `N_x(theta),N_y(theta)` have opposite
parity. The individual parities depend on the chosen real representative of
`theta`, but their difference does not: replacing `theta` by `theta+1` flips
both because `x,y` are odd.

### Proof

The two lifts of `theta` under doubling are

```text
t_j=(theta+j)/2,       j in Z/2Z.                          (5)
```

Both have even-core clearance `phi_C(theta)`. For an odd `z`, the tail is
strictly `1/14`-dangerous at a lift precisely when

```text
||z theta||<1/7,
j=N_z(theta) mod 2.                                       (6)
```

It can kill at most one lift. Thus both lifts over a point of `G_C` are
killed exactly under (4) and opposite parity. Points outside `G_C` are
already strictly killed by the core. This proves the parity criterion.

For a scalar version, THM-774's threshold-free folding identity says

```text
B_(x,y)(theta)
 =max_(j in Z/2Z) min(||x t_j||,||y t_j||)
 =(1-Q_(a,b)(theta))/2.                                   (7)
```

Therefore

```text
M(S)=max_theta min(phi_C(theta),B_(x,y)(theta)).           (8)
```

On `G_C`, the inequality `B<1/14` is exactly `Q>6/7`, proving (3).
The symbols are load-bearing: this theorem uses the **closed** core set
`phi_C>=1/14` and **open** tail/diamond inequalities because it describes a
strict counterexample. THM-769/774 use the dual convention—an open loose core
and closed teeth—to describe a tight equality row at `1/13`. QED.

## 2. Component ownership and the sharp metric box

Assume (3), and lift one connected component `I` of `G_C` to a closed real
interval of length `ell`. The integers `N_x,N_y` in (4) are individually
constant on this lift. Put

```text
D=|N_x y-N_y x|,       g=gcd(x,y).                         (9)
```

Then `D` is a positive odd multiple of `g`, and

```text
ell < min(
  2/(7x),
  2/(7y),
  (x+y-7D)/(7xy)
).                                                        (10)
```

In particular `x+y>7D>=7g`.

Let `R=max(C)`, `mu=M(C)`, and suppose `mu>1/14`. Around a maximizer of
`C`, the closed interval of radius

```text
rho=(mu-1/14)/R                                           (11)
```

lies in `G_C`. Applying (10) gives

```text
x,y < R/(7(mu-1/14)),
D/(xy)+2rho < 1/(7x)+1/(7y).                              (12)
```

If `|C|=s<=11`, settled lower-dimensional LRC gives

```text
x,y < 2(s+1)R/(13-s).                                     (13)
```

For the eleven-speed seam this becomes the finite strict box

```text
x,y<12R,       hence x,y<=12R-1,                          (14)
xy<6R(x+y-7D)<=6R(x+y-7g).                                (15)
```

### Proof

Condition (4) puts the closed interval `I` strictly inside one open `x`-tooth
of radius `1/(7x)` and one open `y`-tooth of radius `1/(7y)`. Connectedness
and uniqueness of nearest integers make `N_x,N_y` constant. Their tooth
centres are `N_x/x` and `N_y/y`, at distance `D/(xy)`. Opposite parity and
odd `x,y` make `N_x y-N_y x` odd and nonzero; divisibility by `g` is
immediate. The intersection of intervals of radii `1/(7x),1/(7y)` and centre
distance `D/(xy)` has length at most the three quantities in (10). Strict
containment gives strict inequalities.

The function `phi_C` is `R`-Lipschitz, which proves (11). Substitute
`ell=2rho` in (10). The LRC bound `mu>=1/(s+1)` gives

```text
7(mu-1/14)>=(13-s)/(2(s+1)),
```

proving (13). At `s=11`, (12) and `rho>=1/(84R)` give (14)--(15). QED.

## 3. Every surviving quotient core is divisor-complete through 14

In the eleven-speed seam, a strict counterexample forces

```text
for every N=2,...,14, some c in C satisfies N|c.           (16)
```

### Proof

Suppose a modulus `N` in that range divides no member of `C`. For every unit
`r mod N`, no `cr` vanishes modulo `N`, so

```text
||c r/N||>=1/N>=1/14.
```

Hence every unit fraction `r/N` belongs to `G_C`. An odd residue `z mod 2N`
is eligible at all those fractions exactly when

```text
7|zr|_N<N       for every r in U_N.                        (17)
```

The complete finite residue table is

| `N` | odd `z mod 2N` satisfying (17) |
|---|---|
| `2,4,6,8,10,12,14` | none |
| `3,5,7,9,11,13` | `z=N mod 2N` |

For even `N`, either tail already fails universal eligibility. For odd `N`,
both tails would be odd multiples of `N`; at `theta=1/N` their nearest
integers are both odd, contradicting opposite parity. Thus (3) is impossible
unless (16) holds. The table is a finite proof by the displayed residue
test; the companion script checks every row independently. QED.

This is the exact `1/14` counterpart of THM-772's divisor transfer, and it is
also a two-tail strengthening of THM-2057's missing-clock idea: instead of
taxing one tail by the missing modulus, the two sheet owners collide in
parity.

## 4. Exact folded measure and sharp cap

Put

```text
g=gcd(a,b)=gcd(x,y),       alpha=a/g,       beta=b/g,
r=alpha-beta=y/g,          s=alpha+beta=x/g,
n=floor(alpha/7+1/2).                                     (18)
```

Then `alpha>beta` are coprime and have opposite parity, and normalized
Lebesgue measure satisfies

```text
measure(H^o_(a,b))
 =2/(rs) sum_(j=1)^n min(
      2r/7,
      2alpha/7-(2j-1)
   ).                                                     (19)
```

Uniformly in the distinct odd tails,

```text
measure(H^o_(a,b))<=4/63.                                 (20)
```

The bound is sharp exactly at the reduced folded pair
`(alpha,beta)=(5,4)`, equivalently at the ratio `x:y=9:1`, and its common odd
dilates. Consequently every strict dyadic-seam counterexample must satisfy

```text
measure(G_C)<=4/63.                                       (21)
```

### Proof

Formula (19) is THM-774's half-grid component calculation with the half-grid
distance budget changed from `2/13` to `1/7`. Multiplication by `g` preserves
measure. The qualifying positive odd half-grid residues are
`1,3,...,2n-1`; their two component lengths sum to (19). Open versus closed
endpoints does not change measure.

The summands give

```text
measure(H^o)<=min(
  4n/(7s),
  2n(2alpha/7-n)/(rs)
).                                                        (22)
```

If `s>=9n`, the first term is at most `4/63`. If `s<9n`, then

```text
4rs>36n(2alpha-9n).
```

When `2alpha>=11n`, this dominates
`126n(2alpha/7-n)`, so the second term in (22) is at most `4/63`. In the
remaining case `2alpha<11n`, the definition of `n` forces `alpha<=12`; the
only nonempty possibilities are `alpha=4,5`. Direct substitution gives
respective maxima `2/49` and `4/63`, the latter only at `beta=4`. This proves
(20). Containment (3) gives (21). QED.

The measure cap is only a necessary scalar test. Pointwise component
containment, ownership words, and the determinant `D` remain load-bearing.

## 5. Finite-exact low-core closure

There is no strict LRC(14) counterexample of the form

```text
2C union {x,y},       C subset {1,...,19},       |C|=11,
gcd(C)=1,             x,y distinct positive odd.           (23)
```

This is uniform over the odd tails; no external speed cutoff is imposed.

### Exact certificate

Part 3 first leaves only `1,365` divisor-complete cores among the
`binom(19,11)=75,582` candidates. For each one, the companion constructs all
closed components of `G_C` with rational endpoints. If `ell_max(C)` is their
maximum length, (10) bounds every universally eligible odd tail by the
intrinsic strict cap

```text
w < 2/(7 ell_max(C)).                                     (22)
```

The script checks every permitted odd `w` exactly, requiring positive slack
at both endpoints of every closed component. It finds **zero** odd tails that
are eligible over all of `G_C`; therefore an opposite-parity pair is already
impossible. The largest intrinsic integer cap is `62`, at

```text
C={1,2,3,4,8,9,10,11,12,13,14},
ell_max(C)=1/220,
```

whose weak-safe set has `24` components. This proves (23).

## 6. Prior-art boundary and residual

- THM-615 first proved the two-odd-tail fold at LRC(14) and closed the AP
  even core; it explicitly left the general eleven-core confinement open.
- THM-769/774 proved the analogous equality-packet and folded diamond at the
  `1/13` threshold; THM-775 develops its dyadic deletion seam.
- THM-2060 supplies the sharp sheet-capacity reduction from an arbitrary
  imprimitive eleven-core to (1).

The residual is now exact: exclude (3) for primitive eleven-sets satisfying
the divisor pins (16), with odd `x,y` in the finite box (14) and determinant
constraint (15). This is finite for each fixed quotient core but not yet
uniform in `R`; THM-2052/2058's finite template intervals are the natural
place to attach these new congruence and component sidecars.

## Assumption challenge and tournament analysis

The proof vertices are the two lift sheets over each component of `G_C`.
Nearest-integer parity is the switch/gauge, and the interval order supplies a
tie Hamiltonian path. In a counterexample the two odd tails form a perfect
opposite-colour ownership at every component, so the local tournament is
transitive and contains no cycle. That tournament forgets tooth widths,
component locations, and the determinant separation `D`; those data are
exactly what (10), (15), and (3) use. The smallest faithful object is therefore
the ordered component list with its two parity words and rational tooth
slacks, not a tournament on runners.

## Computational audit

The companion uses rational and integer arithmetic only. It checks the global
folding identity, the strict parity/diamond equivalence, all unit-residue rows
in (17), the exact measure formula against independent piecewise-linear
integration, the sharp cap over a large reduced box, the determinant parity
and gcd law, and a finite low-core slice. It uses explicit runtime checks in
ordinary and optimized Python; the frozen output ends in `PASS`.
