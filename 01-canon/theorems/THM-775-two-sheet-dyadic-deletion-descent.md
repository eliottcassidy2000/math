---
id: THM-775
title: Two-sheet dyadic deletion descent and the exact Z/4 seam
status: PROVED (elementary deck capacities and residue shells; independently audited)
source: codex-2026-07-14-S9
depends_on:
  - THM-765
  - THM-769
  - THM-772
related:
  - THM-774
  - THM-776
  - HYP-6820
---

# THM-775 — Two-sheet dyadic deletion descent

Put `delta=1/13` and, for a finite set `Q`, write

```text
phi_Q(tau)=min_(q in Q)||q tau||,
G_Q={tau in R/Z:phi_Q(tau)>delta}.
```

Let `A` be a primitive tight twelve-speed set in THM-769's two-sheet
equality packet:

```text
A=2U union {x,y},       |U|=10,       x,y odd.             (1)
```

Thus THM-772 says that `U` is primitive, contains a multiple of every
`m=2,...,12`, and has no multiple of 13.

## 1. Every imprimitive deletion is a dyadic seam

For `u in U`, put `P=U\{u}`.  If `P` is imprimitive, then

```text
gcd(P)=2,       u is odd,       P=2V with gcd(V)=1.         (2)
```

Consequently exactly one of the following alternatives holds.

1. Every one-element deletion of `U` is primitive.
2. `U` has a unique odd member `u` and

   ```text
   U=2V union {u},       |V|=9,       gcd(V)=1.             (3)
   ```

   Deleting `u` is the only imprimitive deletion of `U`.

### Proof

Let `d=gcd(P)>1`.  Since `U` is primitive,

```text
gcd(d,u)=1.                                                (4)
```

The set `A\{2u}=2P union {x,y}` is primitive by THM-765.  Therefore `d`
cannot divide both `x` and `y`; choose

```text
w in {x,y} with d not dividing w,
D=d/gcd(d,w)>=2.                                          (5)
```

Choose a maximizer `tau_0` of `P`.  The lower-dimensional LRC bound gives
`M(P)>=1/10>delta`.  All points

```text
tau_j=tau_0+j/d,       j in Z/dZ,
```

are therefore strictly `P`-safe.  At each `tau_j`, either `u` is
`delta`-dangerous or `tau_j in G_U`; in the latter case THM-769 makes the
chosen odd exception `w` `2delta`-eligible.  Hence the `d` deck points are
covered by the two closed conditions

```text
||u tau_j||<=delta          or          ||w tau_j||<=2delta. (6)
```

By (4), the first phases form a complete `d`-grid.  Its capacity is

```text
alpha(d)=(floor(2d/13)+1)/d <= 1/2,                       (7)
```

with equality only for `d=2`.  The second phases form a `D`-grid with
repetitions and have capacity

```text
beta(D)=(floor(4D/13)+1)/D <= 1/2.                        (8)
```

Here equality in (8) can occur at `D=2` or `D=4`; that extra equality case
does not affect the argument.  Covering the whole deck in (6) forces
`alpha(d)+beta(D)>=1`, so both capacities are one half.  The equality case
in (7) gives `d=2`.  Equation (4) makes `u` odd, every member of `P` is even,
and

```text
2=gcd(P)=gcd(2V)=2gcd(V),
```

so `V` is primitive.  Applying the same argument to any imprimitive
deletion shows that its deleted member must be the unique odd member.  This
proves the dichotomy and uniqueness.  Notice that only one of `x,y` was
chosen in (5); the proof does not claim both are nondivisible by `d`.

## 2. The first seam is an exact `Z/4` ownership tiling

Assume (3).  Then

```text
A=4V union {2u,x,y}.                                      (9)
```

For every `sigma in G_V`, let `N_z(sigma)` denote the unique nearest integer
to `z sigma` whenever the following eligibility bound holds.  In fact

```text
||u sigma||<=2/13,
||x sigma||,||y sigma||<=4/13.                            (10)
```

On the four lifts

```text
t_j=(sigma+j)/4,       j in Z/4Z,
```

the ownership sets are

```text
2u: {j:j=N_u(sigma) (mod 2)},
x : {j_x},       j_x=-x^(-1)N_x(sigma) (mod 4),
y : {j_y},       j_y=-y^(-1)N_y(sigma) (mod 4).           (11)
```

They form a literal disjoint partition of `Z/4Z` of sizes `2+1+1`.
In particular, `j_x` and `j_y` are the two distinct classes of the parity
opposite to the class owned by `2u`.

### Proof

Every speed in `4V` is strictly safe at all four lifts.  The speed `2u` can
kill at most the two lifts of one parity, while each odd speed `x,y` can kill
at most one lift.  These capacities sum to four.  Tightness of `A` requires
all four lifts to be killed, so all three capacities saturate and no two
ownership sets overlap.  Solving the respective congruences gives (10)--(11).

This also says that exactly one of the two lifts

```text
(sigma+k)/2,       k in Z/2Z,
```

lies in `G_U`: `u` blocks the other one.  Equivalently, `u` is a persistent
odd `2/13`-eligible guard over all of `G_V`.

## 3. Divisor transfer and the seam geometry

The quotient `V` again contains a multiple of every

```text
m=2,...,12.                                               (12)
```

Moreover the singleton owners in (11) obey the exact cocycle

```text
N_x(sigma)y-N_y(sigma)x = 2 (mod 4).                      (13)
```

If `B=max(V)`, `mu=M(V)`, and `rho=(mu-delta)/B`, then

```text
2/(xy)+2rho <= 4/(13x)+4/(13y),                           (14)
u <= 20B/3,       x,y <= 40B/3.                           (15)
```

### Proof of divisor transfer

Suppose that `V` has no multiple of some `2<=m<=12`.  Every unit fraction
`a/m` then lies in `G_V`.  The first eligibility condition in (10), imposed
at all unit fractions, has the same `s=2` residue shell as THM-772:

- even `m` admit no odd eligible residue;
- odd `m` force `u=m r` with `r` odd.

Thus `m` is odd.  For `m=3,5,7,9,11`, all-unit `4/13`-eligibility forces an
odd speed to be divisible by `m`.  Indeed the thresholds

```text
floor(4m/13)=0,1,2,2,3
```

are smaller than the maximum least residue in every nonzero gcd orbit
modulo `m` (respectively `1`; `2`; `3`; `4` or `3`; `5`).  Hence

```text
x=m r_x,       y=m r_y,
```

with `r_x,r_y` odd.  At `sigma=1/m`, the guard `2u` owns the odd parity
because `N_u=u/m` is odd, whereas

```text
j_x=-x^(-1)N_x=-m^(-1) (mod 4),
j_y=-y^(-1)N_y=-m^(-1) (mod 4).
```

The two odd exceptions therefore own the same singleton, contradicting the
disjoint partition (11).  This proves (12).

### Proof of the metric statements

The two singleton classes differ by `2 mod 4`; multiplying their difference
by the odd unit `xy` gives (13).  Around a maximizer of `V`, the
`B`-Lipschitz property supplies an interval of radius `rho` contained in
`G_V`.  On it the `x`- and `y`-eligibility teeth have respective radii
`4/(13x)` and `4/(13y)`.  Their centres differ by at least `2/(xy)` by (13),
which proves (14).

Since `V` has nine speeds, `mu>=1/10`, and hence

```text
rho >= (1/10-1/13)/B=3/(130B).
```

Containment of that interval in the `u` tooth of radius `2/(13u)` and in
the `x,y` teeth of radii `4/(13x),4/(13y)` gives (15).

## 4. The full dyadic deletion tower

The preceding seam iterates.  There is a finite chain

```text
U=Q_0,
Q_i=2Q_(i+1) union {h_i},       h_i odd,                  (16)
```

ending at a quotient `Q_r` whose one-element deletions are all primitive.
Every `Q_i` is primitive and contains a multiple of every `m=2,...,12`.
At each displayed seam:

1. deleting `h_i` is the sole imprimitive deletion of `Q_i`;
2. `h_i` is `2/13`-eligible at every point of `G_(Q_(i+1))`;
3. of the two lifts of each point of `G_(Q_(i+1))`, exactly one lies in
   `G_(Q_i)` and the other is blocked by `h_i`.

If `Q_0` is already hereditarily primitive the chain has length zero.  If it
is not, Sections 1--3 construct its first seam and `Q_1=V`.

Unwinding the chain gives the explicit 2-adic normal form

```text
U=2^r Q_r union {2^i h_i:0<=i<r}.                        (16a)
```

Hence a depth-`r` tower (`0<=r<=7`) has exactly one quotient speed of each
2-adic valuation `0,1,...,r-1`; all other `10-r` speeds are divisible by
`2^r`.  This is an equality of disjoint sets because every `h_i` is odd.
It follows immediately by substituting (16) recursively, and makes the
dyadic obstruction directly testable without reconstructing the intermediate
quotients.  The terminal core has at least three members: a primitive
divisor-complete singleton is impossible, while a hereditarily primitive
two-set would require both singleton deletions to equal `{1}`, contradicting
distinctness.  Since `|Q_r|=10-r`, this proves the stated depth bound.

### Inductive proof of the gcd descent

Suppose

```text
Q_(i-1)=2Q_i union {h_(i-1)}                              (17)
```

has the three stated invariants, and let `P=Q_i\{q}` have gcd `d>1`.
Primitivity of `Q_i` gives `gcd(d,q)=1`.  Also `d` cannot divide the preceding
guard `h_(i-1)`: otherwise

```text
Q_(i-1)\{2q}=2P union {h_(i-1)}
```

would be imprimitive, contradicting the uniqueness of the deletion of
`h_(i-1)` one level above.

At a maximizer of `P`, all `d` translation lifts are strictly `P`-safe by
the lower-dimensional LRC bound.  Each lift is covered either by the
`delta`-danger tooth of `q` or, when it lies in `G_(Q_i)`, by the
`2delta`-eligibility tooth of `h_(i-1)`.  The same capacities (7)--(8) force
`d=2`.  Thus `q` is the unique odd element,

```text
Q_i=2Q_(i+1) union {q},       gcd(Q_(i+1))=1.              (18)
```

The deletion of `q` is the sole imprimitive deletion by the same argument.

For `sigma in G_(Q_(i+1))`, the odd `q` can block at most one of its two
lifts.  It must block one: if both lifts lay in `G_(Q_i)`, the preceding odd
guard `h_(i-1)` would have to be `2/13`-eligible at two antipodal phases,
which is impossible because `4/13<1/2`.  Hence exactly one lift is blocked,
and `q` becomes the next persistent guard.  This regenerates all induction
invariants.  Cardinality drops at every seam, so the process terminates.

### Inductive proof of divisor transfer

It remains to justify that divisor-completeness survives every later seam.
Suppose `Q_(i+1)` misses `m in {2,...,12}`.  Every unit `a/m` lies in its
loose set, so the newest odd guard `h_i` is all-unit `2/13`-eligible.  The
same shell used above excludes even `m` and, for odd `m`, gives

```text
h_i=m r,       r odd.                                    (19)
```

For each unit `a mod m`, the lift blocked by `h_i` has even numerator over
`2m`; the unique safe lift has odd numerator.  As `a` varies, these safe
lifts are exactly all unit fractions modulo `2m`, and they lie in `G_(Q_i)`.
The preceding odd guard `h_(i-1)` would therefore be `2/13`-eligible at
every unit fraction of denominator

```text
2m in {6,10,14,18,22}.                                   (20)
```

This is impossible.  For completeness, the danger thresholds in these five
rows are `0,1,2,2,3`; according to the possible odd gcd orbit, the maximum
least residues are respectively

```text
{1,3}, {3,5}, {5,7}, {7,3,9}, {9,11},
```

all strictly larger than the corresponding threshold.  This elementary
five-row shell is beyond THM-772's stated denominator range and is included
here explicitly.  Thus every quotient in (16) remains divisor-complete.

## 5. Ownership, tournaments, and what the quotient preserves

The theorem-facing object is a saturated sheet-incidence hypergraph.  At the
first seam the four sheet vertices carry the literal disjoint `2+1+1`
partition (11).  At later seams, assign the newly blocked half to the newest
odd guard and recursively assign the safe half below it.  This gives a
canonical disjoint ownership partition on a binary sheet tree.

It is important not to overstate that conclusion: old guards may also be
dangerous on parts of the newly blocked half.  Beyond the first `Z/4` seam,
the assigned ownership is disjoint, but the raw danger incidences need not be.

A runner tournament may orient an edge by which guard owns the lower-labelled
sheet, with circle reflection as the switch/gauge and binary sheet order as
the tie Hamiltonian path.  Those tournaments are transitive at a fixed seam;
they discard the capacity saturation, the safe-half recursion, and the unit
columns used in the divisor transfer.  The binary sheet tree plus its
speed-incidence sidecar is the smallest quotient here known to preserve the
LRC covering predicate.

This also records the challenged assumption: the natural recursive vertices
are neither runners nor arcs, but proof obligations indexed by sheets.  The
absolute speeds and atom positions can be forgotten; parity ownership,
eligibility radii, and the safe-child map cannot.

## 6. Frontier effect

The two-sheet branch is no longer an arbitrary persistent parity cover.  Any
failure of hereditary primitivity is forced into a finite dyadic tower whose
quotients stay primitive and divisor-complete, ending at a hereditarily
primitive base.  This does not yet exclude a tower or its terminal base
uniformly.  It converts the remaining problem into one of proving that no
primitive divisor-complete terminal core can support the folded cover of
THM-774, or that the canonical binary ownership tree cannot persist at
unbounded height.
