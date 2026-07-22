---
id: THM-2077
title: "Terminal Kakeya needle, all-depth relative height, and recursive quarter escape"
status: >
  PROVED. THM-2075 lifts a terminal maximizer interval through the exact
  safe-child homeomorphism. Its length gives uniform relative-height bounds
  for every dyadic guard and both original odd tails, including the folded
  owner-determinant sidecar. Independently, every quotient level must fail
  the THM-2072 quarter-fan certificate, yielding an explicit residue-class
  ratio escape at every depth. This is a structural reduction, not LRC(14).
source: codex-2026-07-21-LRC-tangent-referee
depends_on:
  - THM-1234
  - THM-2072
  - THM-2073
  - THM-2075
  - THM-2076
related:
  - THM-2061
  - THM-2078
  - THM-841
  - THM-848
  - THM-1196
  - THM-1260
---

# THM-2077 -- terminal needle and recursive quarter escape

Put `delta=1/14`.  Retain a hypothetical strict dyadic seam and its tower
from THM-2073:

```text
S=2C union {x,y},             M(S)<delta,
C=Q_0,
Q_i=2Q_(i+1) union {h_i},     0<=i<r,                 (1)
```

where `x,y` and every `h_i` are odd, every `Q_i` is primitive and contains a
multiple of each integer `2,...,14`, and `Q_r` is hereditarily primitive.
The tower has `|Q_i|=11-i`.  THM-2073 gives `r<=8`, and the guard-capacity
rank floor of THM-2076 sharpens this to

```text
0<=r<=5.                                                (1a)
```

Write

```text
mu=M(Q_r),             B=max(Q_r),
rho=(mu-delta)/B.                                      (2)
```

The quantity `rho` is positive.  This theorem proves two independent
necessary structures: a metric needle with an all-depth height box, and a
scale escape forced at every quarter anchor.

## 1. The terminal maximizer supplies one literal Kakeya needle

Let `sigma_*` maximize

```text
phi_(Q_r)(sigma)=min_(q in Q_r)||q sigma||.
```

Since this function is `B`-Lipschitz, the closed circle interval

```text
I_r={sigma:dist(sigma,sigma_*)<=rho}                  (3)
```

lies in `G_(Q_r)` and has length `2rho`.  THM-2075 says that the inverse of
doubling is affine with one constant sheet address on every connected safe
component.  Iteratively lifting (3) therefore gives closed intervals

```text
I_i subset G_(Q_i),
|I_i|=2^(-(r-i))|I_r|=2^(i-r+1)rho,    0<=i<=r.       (4)
```

In particular the original core safe set contains the interval

```text
I_0 subset G_C,                |I_0|=2^(1-r)rho.      (5)
```

This is an actual one-dimensional Kakeya needle: it retains its affine
offset, its binary sheet address, and its inherited endpoint owners.  It is
not merely the measure identity `measure(G_C)=2^(-r)measure(G_(Q_r))`.

### Proof

For `sigma` within circle distance `rho` of `sigma_*` and every `q in Q_r`,

```text
||q sigma|| >= ||q sigma_*||-q rho
             >= mu-B rho=delta.                       (6)
```

Thus (3) holds.  Since `rho<1/2`, it is a proper circle interval.  It lies in
one connected component of `G_(Q_r)`.  The componentwise affine inverse in
THM-2075 halves its length at each of the `r-i` lifts and keeps the whole
interval on one sheet, proving (4)--(5). QED.

### Exact address and tail H-drift formula

The binary address in the preceding proof has an explicit nearest-integer
form.  On one terminal component choose compatible real lifts and put

```text
sigma_r=sigma,
epsilon_i=1-N_(h_i)(sigma_(i+1)) mod 2,
sigma_i=(sigma_(i+1)+epsilon_i)/2,       0<=i<r.       (6a)
```

The blocked child has bit `N_(h_i)(sigma_(i+1)) mod 2`, so `epsilon_i` is
the safe bit.  Indeed, for a candidate child with bit `k`,

```text
h_i(sigma_(i+1)+k)/2
```

is strictly `1/14`-dangerous exactly when
`N_(h_i)(sigma_(i+1))+h_i k` is even.  Since `h_i` is odd, that is exactly
`k=N_(h_i)(sigma_(i+1)) mod 2`.  Eligibility keeps each nearest integer,
hence each bit, constant on the relevant component.

Consequently the terminal address is

```text
a=sum_(i=0)^(r-1) 2^(r-1-i)epsilon_i,
sigma_0=(sigma+a)/2^r.                                 (6b)
```

For either original tail `z in {x,y}`, let `N_z` be its constant nearest
integer on the lifted top component and define

```text
K_z=2^r N_z-za.                                       (6c)
```

The strict top-level tail inequality is then exactly

```text
|z sigma-K_z|<2^r/7.                                  (6d)
```

Moreover complementary top-sheet ownership gives

```text
K_x y-K_y x
 =2^r(N_x y-N_y x)
 ==2^r mod 2^(r+1).                                   (6e)
```

The congruence is signed and unambiguous: `x,y` are odd and `N_x,N_y` have
opposite parity, so the parenthesized determinant is odd.  Thus the literal
H-drift state is not only the real interval displacement.  It is the packet

```text
(a,K_x,K_y),                                          (6f)
```

whose determinant carries exact dyadic valuation `r`.  Forgetting `a`
destroys the individual `K_z`, even though it cancels from their determinant.

## 2. All-depth relative-height bounds

Every tower guard and each original tail is forced to contain one of the
intervals (4)--(5) in a single strict danger tooth.  Consequently

```text
h_i < 2^(r-i-1)/(7rho),                0<=i<r,         (7)
2^i h_i < 2^(r-1)/(7rho),              0<=i<r,         (8)
x,y < 2^r/(7rho).                                      (9)
```

The settled lower-dimensional lonely-runner theorem applied to the
`11-r` speeds of `Q_r` gives

```text
mu>=1/(12-r),
rho>=(r+2)/[14(12-r)B].                               (10)
```

Hence the explicit terminal-scale box is

```text
2^i h_i < 2^r(12-r)B/(r+2),                           (11)
x,y < 2^(r+1)(12-r)B/(r+2).                           (12)
```

All thirteen speeds of `S` therefore obey

```text
max(S)
 <=2^(r+1)B max(1,(12-r)/(r+2)).                      (13)
```

The weak inequality in (13) allows the terminal speed `2^(r+1)B`; the guard
and tail bounds themselves are strict.

There is also a folded-owner sidecar.  Use (6c) on the terminal component
containing `I_r` and put

```text
Delta=|K_x y-K_y x|,            D=Delta/2^r.
```

Then `D` is a positive odd multiple of `gcd(x,y)`, equation (6e) fixes the
exact dyadic valuation of `Delta`, and

```text
D/(xy)+2^(1-r)rho < 1/(7x)+1/(7y).                   (14)
```

### Proof

THM-2073 makes `h_i` strictly `1/7`-eligible at every point of
`G_(Q_(i+1))`.  Thus the connected compact interval `I_(i+1)` is contained
strictly in one component of

```text
{sigma:||h_i sigma||<1/7},
```

whose length is `2/(7h_i)`.  Equations (4) and strict containment give

```text
2^(i-r+2)rho=|I_(i+1)|<2/(7h_i),
```

which is (7); multiplying by `2^i` gives (8).

THM-2061 makes each odd tail strictly `1/7`-eligible on all of `G_C`.
Applying the same tooth-length argument to (5) proves (9).  Since
`|Q_r|=11-r`, the known LRC bound is `mu>=1/(|Q_r|+1)=1/(12-r)`.  Subtracting
`1/14` proves (10), and substitution in (8)--(9) proves (11)--(12).
The speeds contributed by `2C` are `2^(r+1)q` with `q in Q_r` and
`2^(i+1)h_i`; taking their maximum together with the tails proves (13).

Finally, (6d) puts the terminal interval `I_r` strictly inside two open
intervals of radii `2^r/(7x)` and `2^r/(7y)`.  Their centres are separated by
`Delta/(xy)`, so their intersection has length at most

```text
2^r/(7x)+2^r/(7y)-Delta/(xy).
```

The compact interval `I_r` has length `2rho` and is strictly contained in
that intersection.  Divide the resulting inequality by `2^r` to obtain
(14).  Equation (6e) makes `D` odd and positive after taking the absolute
value; divisibility by `gcd(x,y)` is immediate. QED.

The open/closed convention is load-bearing and favourable.  The intervals
`I_i` are **closed weak-safe** intervals.  Guards and tails are **strictly
dangerous** on them.  A compact interval inside one open tooth has strictly
smaller length than that tooth, which is exactly why (7), (9), and (14) are
strict.  Replacing either side silently by a measure statement would lose
this endpoint gain.

### Deepest-lane pair-capacity saturation

If `r=5`, then `Q_5` has six speeds.  The sharp five-comb pair floor in
THM-1234 gives its six-danger-comb union bound

```text
measure(union_(q in Q_5){||q sigma||<1/14})<=212/273.
```

Consequently

```text
measure(G_(Q_5))>=61/273.                              (14a)
```

The final guard law puts this closed safe set inside

```text
E_(h_4)={sigma:||h_4 sigma||<1/7},
measure(E_(h_4))=2/7=78/273.                           (14b)
```

Compact/open separation and (14a)--(14b) therefore give the exact
near-saturation ledger

```text
61/273<=measure(G_(Q_5))<78/273,
0<measure(E_(h_4) minus G_(Q_5))<=17/273.              (14c)
```

Equivalently, the terminal safe set occupies at least `61/78` of the final
guard's entire eligibility region.  Since that region consists of `h_4`
equal teeth, at least one guard tooth contains terminal safe measure

```text
at least 61/(273h_4).                                  (14d)
```

THM-2075 transports (14a) to the original core as

```text
measure(G_C)>=61/(273*2^5)=61/8736.                    (14e)
```

This is the precise surviving contribution of the global five-comb/Fano
pair functional to the tower.  It supplies a high-occupancy guard tooth, not
a single long component, a located four-prefix survivor, or a `j=4` flood
chronology.  Turning (14d) into one of those objects requires an endpoint or
component-count sidecar.

## 3. Every quotient level has a quarter-anchor escape

No safe set `G_(Q_i)` contains an antipodal pair:

```text
there is no theta with
theta,theta+1/2 in G_(Q_i),             0<=i<=r.       (15)
```

For `i=0`, either original odd tail would be strictly `1/7`-dangerous at both
points, contradicting

```text
||z(theta+1/2)||=1/2-||z theta||.                     (16)
```

For `i>=1`, the preceding odd guard `h_(i-1)` gives the same contradiction
because it is strictly `1/7`-eligible everywhere on `G_(Q_i)`.

Let `c_i` be the smallest member of `Q_i` divisible by four; it exists by
divisor completeness.  The quarter-fan corollary of THM-2072 and (15) imply
that at least one of the following holds at every level:

```text
some c in Q_i with 4|c satisfies          c>7c_i;
some odd c in Q_i satisfies               c>5c_i/2;
some c in Q_i with c=2 mod 4 satisfies    c>6c_i.       (17)
```

In particular

```text
max(Q_i)>5c_i/2                                           (18)
```

at every level.  This is a recursive obstruction, not only a terminal-core
condition.

For `i<r`, let

```text
e_i=min{q in Q_(i+1):q is even}.
```

Since `Q_i=2Q_(i+1) union {h_i}` and `h_i` is its unique odd member,
`c_i=2e_i`.  Thus (17) has the exact child-level form

```text
some even q in Q_(i+1) has q>7e_i,
or h_i>5e_i,
or some odd q in Q_(i+1) has q>6e_i.                  (19)
```

### Proof

Identity (16) holds for every odd integer `z`.  Two numbers summing to
`1/2` cannot both be less than `1/7`, proving (15) in the two cases just
described.  If none of the alternatives (17) held, all three bounded-fan
hypotheses of THM-2072 would hold with anchor `c_i`; that theorem would put
an antipodal safe pair in `G_(Q_i)`, contradicting (15).  This proves
(17)--(18).  The residue classes in `Q_i=2Q_(i+1) union {h_i}` translate
exactly as stated above, proving (19). QED.

## 4. Frontier effect

The combination of (11)--(12) and (17) is the useful new restriction.  A
putative seam tower cannot be both residue-balanced and arbitrarily loose:

```text
every level must make a named mod-4 scale escape,
while every guard and both tails remain in one explicit box over
the terminal maximum B.                               (20)
```

This converts the non-hereditarily-primitive lane into a terminal-core
problem with two finite sidecars: the component address word from THM-2075
and the relative-height box above.  THM-2076 has already removed depths six
through eight, so (17) and (19) need be iterated at most five times.  This
does not upper-bound `B`, classify the hereditarily primitive terminal core,
or prove LRC(14).  THM-2078 separately uses the present height box to exclude
the whole finite slice `B<=24`; the remaining nontrivial terminal lane has
`B>=25`.

The old viewpoints now have precise scopes.

* **Toothpick/ladder.** THM-841 refutes homogeneous self-similarity of the
  Farey count.  Here the exact self-similar object is instead the legal
  operation `(I,a)->(I+a)/2`: lengths halve only because the component
  address is retained.  Discarding that bit recovers the old false-count
  failure mode.
* **H-drift.** THM-848 distinguishes closure of an averaged observable from
  closure of the transition state.  Likewise the measure-halving scalar is
  not enough for (7)--(14); the componentwise sheet section is the required
  state coordinate.
* **Fano/`chi_7`.** THM-1260 proves local colour-surjectivity even after a
  placed fork.  Those colours cannot determine the nearest-integer sheet bit
  or any alternative in (17).  A Fano quotient therefore supplies no missing
  implication here without address/owner incidence.
* **`j=4` flood tail.** THM-1196 and its positioned successors control how
  two combs can flood a four-prefix survivor.  The interval `I_i` is a
  candidate consumer for that technology, but (1) alone supplies no
  six-comb cover on the same located interval.  Importing the flood bound
  without such an incidence sidecar would change the predicate.
* **Kakeya.** The needle in (3)--(5) is literal one-dimensional interval
  transport.  No planar Kakeya dimension or maximal-function statement is
  used.

For Tournament Analysis, orienting quotient levels by depth gives a
transitive tournament with one Hamiltonian path and no nontrivial SCC.  It
preserves only recursion order.  It destroys the sheet address, interval
offset, mod-4 residue class, endpoint owner, strict tooth, and folded
determinant.  The faithful carrier is therefore

```text
(terminal component and owners; dyadic address word; lifted interval;
 guard/tail tooth centres; quarter-escape witness).                    (21)
```

There is no intrinsic binary orientation among these obligations, so the
tournament shadow does not advance the proof. QED.
