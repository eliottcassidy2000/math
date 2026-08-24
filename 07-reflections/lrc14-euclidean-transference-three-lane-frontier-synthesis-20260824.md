---
status: >
  SESSION SYNTHESIS. THM-4009 is CITED + PROVED ALGEBRA + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. The parallel-class and divisor-anchor results
  are proved/finite-exact; the Fourier and x7 lanes are finite-exact bounded
  laboratories. The diagonal-polar refinement and the first-kind
  character-sensitive Foster theorem are PROVED ALGEBRA + VERIFIED-EXACT.
  Their extensions to owner separation and arbitrary lattices remain research
  targets. LRC(14) remains OPEN.
source: root LRC(14) frontier session, 2026-08-24
tags:
  - lonely-runner
  - lrc14
  - geometry-of-numbers
  - graver
  - transference
  - rank-eleven
  - fourier
  - lift-tower
---

# LRC(14): Euclidean transference and the three surviving sidecars

## Verdict

The session produced one new canonical reduction, two proved refinements, and
three sharply typed continuations. It did not prove LRC(14).

[THM-4009](../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md)
replaces the old `l1<=356` flatness row by a primitive Graver relation

```text
a dot n=0,       ||a||_2<14,       sum a_i^2<=195,
||a||_1<=50,     max |a_i|<=13.                              (1)
```

The exact speed-dependent threshold is

```text
||a||_2 < 14 sqrt(m^2+M^2)/(m+M),
m=min n_i,       M=max n_i.                                 (2)
```

The reduction leaves `55,459` primitive absolute coefficient histograms,
not the effectively unstructured `l1<=356` ball. Its support-two branch has
47 coprime ratios, 3,666 unoriented ratio/support packets, and 7,332 oriented
labelled assignments modulo global sign. An independent implementation
caught and repaired the orientation distinction and reproduced every other
constant, census, strictness step, parity statement, and rank-eleven cap.

The obstruction remaining after all four lanes is not another missing scalar.
It appears in three equivalent-looking but not yet identified forms:

```text
lattice geometry:   a prescribed odd character of the dual lattice,
rank-eleven code:   a crossing row with an endpoint owner,
Fourier/lift code:  a cross-order or cross-sheet phase word.                (3)
```

No implication among the three formulations is asserted without its map and
sidecars.

[THM-4014](../01-canon/theorems/THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression.md)
now improves the necessary relation further: there is a (possibly different)
Graver row with

```text
||a||_1<=49,      |a_M|<=7,
a_M!=0  =>  sum a_i^2<=193.                              (3a)
```

This is localization rather than closure. The sharp character theorem on
Voronoi-first-kind lattices identifies exactly what a successful parity
upgrade would look like, while the speed-weighted circuit explains why that
theorem does not yet apply to the LRC projected lattice.

## Inheritance pass

- **Closest proved mechanism.** THM-3743 identifies the quotient dual with
  the integer speed-relation lattice and converts zonotope width to `l1`.
  THM-2052 gives the rank-eleven/rank-twelve fork. THM-3818 identifies the
  two-component equality branch. THM-4002 shows that an unlabelled relation
  fibre can have opposite Fourier covariance.
- **Canonical hostile.** AP13 has the norm-square-six Graver relation
  `(1,-2,1,0,...)` and is nevertheless lonely at `t=1/14`. Short resonance is
  necessary, not sufficient.
- **Corrected near miss.** The orientation census is `3,666` only before
  assigning `p<q` to the two coordinate labels; the labelled oriented count
  is `7,332`. The stale THM-3818 two-component companion also used the old
  table-free common-scale filter: MISTAKE-486 restores two internal edges and
  changes `25 -> 27` without changing any theorem implication.
- **Least-used sidecar.** The zonotope center is a nonzero order-two point of
  the projected torus. Biduality therefore forces an odd-sum relation, but
  ordinary transference does not prove that the Euclidean-shortest relation
  is odd.

## Concept board

| live object | exact retained predicate | first information loss | cheapest hostile test |
|---|---|---|---|
| projected lonely zonotope | closed lattice-point avoidance | owner and phase | AP13 |
| square-norm-195 Graver atlas | one exact speed relation | odd character and component cut | AP13 / internal norm-10 rows |
| rank-eleven two-component code | all bounded support-three rows | a crossing edge becomes a coloop | canonical `2+11` hostile |
| prefix-order Fourier bouquet | sign, scale, order, cross-phase | endpoint atom | AP and Goddyn--Wong |
| first `x7` coefficient fibre | every new clock word | old clocks, gcd ownership, deck | `q=86` SAT / `q=382` gcd ghost |
| order-two torus character | existence of an odd relation | Euclidean length / Smith index | all-odd primitive rows are already lonely |

## 1. Why the Euclidean constant is exactly fourteen

Write `V=n^perp`, `Lambda=pi(Z^13)`, and

```text
Z(n)=c+(3/7)K,       c=(1/2)pi(1),
K=pi([-1,1]^13).                                           (4)
```

Every Euclidean unit vector of `V` belongs to the cube and projects to
itself, so `B_V(0,1) subset K`. More precisely,

```text
inrad(K)=min_(u in V,||u||2=1)||u||1
        =(m+M)/sqrt(m^2+M^2).                              (5)
```

The minimizer is supported on the minimum and maximum speeds. A counterexample
makes the closed ball about `c` disjoint from `Lambda`, hence its distance is
strictly greater than `(3/7)inrad(K)`. Banaszczyk's exact Euclidean inequality

```text
2 mu(L) lambda_1(L*) <= dim L=12                           (6)
```

then gives (1)--(2). Closed-boundary disjointness, not merely an empty
interior, supplies the strict `<14` needed for the integer cap `195`.

Two older scale facts locate the live rows around (2). Write `M_2` for the
second-largest speed. THM-1043 already makes
every row with `M/m<=13` lonely at `t=1/(14m)`. THM-1008 makes every row with
`M/M_2>=13` lonely by deleting the largest speed and perturbing an LRC(13)
witness. Thus a counterexample must satisfy

```text
M/m>13,                    M/M_2<13.                        (7)
```

It is globally spread but cannot be a single dominant-top extension. The
`R=13` line of THM-4009's bounded-spread table is therefore calibration, not
a surviving case.

## 2. Exact rank-eleven join

Let `W` be the span of THM-2052's support-at-most-three height-`91^6` rows.
In rank eleven, the short row `a` has only two possibilities:

```text
a outside W  => rank twelve and an explicit 134-digit speed cap,
a inside W   => W contains a primitive row of square norm <=195.            (8)
```

Support at most three is automatically inside `W`, because height `13` is far
below `91^6`. Hence only support at least four can be the missing rank
increment.

The two-component quotient does not turn (8) into a crossing theorem. Its
column matroid is exactly two parallel classes. Every triple circuit is
internal; adding a graphic edge across the cut raises rank and is a bridge,
hence belongs to no cycle. This proves why basis exchange and ordinary
augmentation stall: the desired row is transverse to the entire circuit
algebra.

The 16 scale-one residual tail pairs interact with (1) as follows. Fourteen
have tail relation square norm at most `195`; only `(8,21)` and `(9,11)` are
larger, with squares `505` and `202`. The scale-two `(1,9)` tail has square
`82`. None closes, because the eleven-body component can absorb a still
shorter internal row; explicit controls have internal minimum square `10`.

There is nevertheless one new crossing lemma. For same-component speeds
`A=sd`, `M=sU` and opposite-component `N=tp`, if

```text
A divides M,       A divides N,       0<N<M,       M/A<=Q^2,                (9)
```

Dirichlet pigeonholing gives `aA+bM-cN=0` with coefficient height at most
`Q`. In the rank-eleven equality branch this is forbidden. In original
coordinates the load-bearing conditions are `d|U` and `sd|tp`; omitting
`A|M` makes the claimed coefficient nonintegral. The scale-two parity escape
is visible exactly as the failure `2|t` for its odd `t`.

The incoming THM-4003/4004 component-erosion and divisor-profile theorems were
intersected with the new metric, not merely appended. In the scale-two row
`2u direct-sum t(1,9)`, the global maximum and runner-up are `9t` and `2U`,
so THM-4014 automatically improves to `||a||_1<=48`. Its selected top-pair
support-two subbranch is empty: `p+q<=13` would be a height-at-most-twelve
crossing row, forbidden by the two-component code. This deletes a support
branch but no `(t,U)` cell, because the tail-internal relation has `l1=10`
and weighted energy below `94`. In the `t<U` divisor lane, a body containing
`1,2,3` can similarly absorb the metric through `1+2-3=0` while every prime
of `t` hits no body coordinate. The exact missing statistic is therefore

```text
min(lambda_H(body),lambda_H(tail));                    (9a)
```

only a lower bound of `196` on both internal weighted minima would force the
transference-selected row to cross.

## 3. Fourier verdict: order matters, owner still matters more

On the strict row `V_38={1,...,11,13,38}`, the 30 Euclidean-shortest rows are
the additive triples `x+y=z`, all of square norm three. Every one of the 60
single-prefix-order bouquets is rigorously positive definite, and centered
Fejer orders one through seventeen stay positive. Retaining both prefix
orders creates exactly two negative matrices. The clean certificate belongs
to `2+6=8` and has

```text
Q=-5.911372570568...<0,       integer coefficient square norm=394.          (10)
```

Thus the gain comes from cross-order phase, not relation shortness. AP and
Goddyn--Wong remain positive on all 30 two-order bouquets. Their valid
boundary witness is an atomic endpoint event with opposite-slope owners, so
an absolutely continuous Toeplitz quadratic cannot see it. THM-4002's two
rows with identical norm-three relation fibre and opposite covariance remain
a second hostile.

The next lawful Fourier object is therefore an owner-labelled, all-order
prefix SDP, not a larger one-order scalar bouquet.

## 4. The first `x7` fibre is a grouped covering code

For a fixed parent `v mod q`, `gcd(q,7)=1`, write children as

```text
x_i=v_i+q a_i,                     a_i in F_7.             (11)
```

Every new clock `h/(7q)`, `7 not|h`, forbids exactly one digit in each
coordinate, except at an exact endpoint where it forbids none. An improper
child is therefore a grouped positive-CNF word hitting every clock clause.
The separate divisor-seven condition supplies one distinguished digit per
coordinate and forbids any twelve of the thirteen from being selected
simultaneously.

The finite controls separate the roles sharply:

- the AP fibre at `q=86` is SAT and has an explicit improper ansatz child;
- the AP fibre at `q=382` is UNSAT in two engines on 1,267 clauses;
- deleting the divisor sidecar produces an all-divisible false survivor.

The scalar best-column capacity and shortest-clause census are both blind:
the SAT and UNSAT controls have the same boundary-clause profile. As in the
Fourier lane, full cross-phase incidence plus an ownership sidecar is the
decisive object.

The source boundary is important. The published Sungkawichai--
Trakulthongchai paper explicitly identifies computation of `I(13,p,1)` as
the extension bottleneck and gives the general `c^k` lift cost. It does not
publish a fourteen-runner lift diagram or state that a particular `x7` step
is unavoidable. The `x7` fibre is a repository inference from the general
architecture and `14=2*7`, not a quoted theorem.

## 5. A new exact polar/ellipsoid lane

[THM-4014](../01-canon/theorems/THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression.md)
canonizes this lane.

The projected cube has the exact polar

```text
K polar = {u in V : ||u||_1<=1}.                         (12)
```

Because all speeds are positive, its vertices are precisely

```text
u_ij=+/- (n_j e_i-n_i e_j)/(n_i+n_j),       i<j.        (13)
```

This gives a component-sensitive replacement for the scalar Euclidean ball.
Let `H=diag(h_1,...,h_13)` be positive and suppose every pair satisfies

```text
(h_i-1)n_j^2+(h_j-1)n_i^2 <= 2 n_i n_j.                 (14)
```

Then all vertices (13) lie in `u^T H u<=1`; by convexity, so does the entire
polar. Polar reversal gives an ellipsoid inside `K`. The operator to invert is
the compression `A=P_V H|_V`, not the ambient diagonal matrix; this distinction
is load-bearing. The exact identity

```text
<A^(-1)x,x>=x^T H^(-1)x
              -(n^T H^(-1)x)^2/(n^T H^(-1)n)           (14a)
```

repairs the tempting but false ambient-inverse identification. Applying
Euclidean transference after `A^(-1/2)` proves that a counterexample forces a
Graver relation with

```text
sum_i h_i a_i^2 <196.                                   (15)
```

Two useful feasible banks are

```text
h_i=1+2 n_i/(M+M_2)                         for all i,   (16)

h_i=1+2 n_i/max_(j!=i)n_j,  h_j=1 (j!=i)    for one i.  (17)
```

Here `M_2` is the second-largest speed. For the maximum coordinate, (17)
gives

```text
a_M=0, or  sum a_i^2 + 2(M/M_2)a_M^2 <196.              (18)
```

The kernel equation and exact integer optimization sharpen (18) to (3a).
Moreover `||a||_1=49` is possible only in the rigid shape

```text
|a_M|=1,       |a_j|=4 for all j!=M,
M/M_2<3/2,     4 divides M.                              (18a)
```

Thus `M/M_2>=3/2` or `4` not dividing `M` improves the selected relation to
`||a||_1<=48`. The simultaneous bank (16) gives one row satisfying

```text
(M+M_2) sum a_i^2+2 sum n_i a_i^2<=196(M+M_2)-1,       (18b)
```

and reduces its top-pair support-two branch from 47 to 28 ratios. These are
proved necessary relation banks, not crossing or owner theorems: the same
very short internal component row may absorb different coordinate boosts.
The next experiment is an exact separation program whose objective penalizes
all internal relations while (13) supplies the 78 containment constraints.

## 6. The character-sensitive transference target

Since `2c in Lambda` but `c notin Lambda`, define

```text
lambda_odd=min{||b||_2 : b in Lambda*, <b,c>=1/2 mod Z}. (19)
```

THM-4009 proves that (19) is finite, but not that it is below fourteen. The
clean desired inequality is a relative/coset analogue of Banaszczyk:

```text
dist(c,Lambda) lambda_odd <= d/2.                         (20)
```

In coordinates, if `Lambda=B Z^d` and `c=Bk/2`, it asks for a bound between
the shortest vector in the primal parity coset `k+2Z^d` and the shortest
dual vector `m` with `m dot k` odd. Ordinary transference can select the
trivial character, so it does not prove (20). Applying transference merely to
the index-two over- or sublattice also loses the quotient sector. A proof
would require relative minima of a nested lattice pair or a signed Gaussian
argument; a counterexample would precisely explain why parity and shortness
cannot be fused generically.

This is the highest-value theorem search because a positive answer with the
same constant would upgrade (1) to an odd-sum short Graver row immediately.
Until a primary theorem, proof, or exact counterexample is found, (20) is
**OPEN**, not an imported consequence of Banaszczyk.

There is now a sharp positive model. If `Lambda` has an obtuse superbase, its
Selling parameters are conductances of a connected graph. A parity class is
a cut of capacity `C`; every crossing edge gives an odd dual vector of squared
norm equal to its effective resistance. Weighted Foster,

```text
sum_e conductance(e) R_eff(e)=d,                         (20a)
```

then proves

```text
dist(c,Lambda) lambda_odd <= sqrt(d)/2.                  (20b)
```

The constant is sharp for `Z^d` with the all-half centre, and the theorem
covers every lattice in dimensions at most three. The LRC frame is close but
not graphic in the required integral sense:

```text
g_i=pi(e_i),     g_i dot g_j<0,
sum_i n_i g_i=0.                                        (20c)
```

The speed-weighted circuit (20c) is not an unweighted unimodular obtuse
superbase. Its Smith-index payload is exactly where the crossing-edge
functional can fail to lift to an integral odd character. A general proof of
(20), or a weighted-circuit Foster theorem retaining that index, remains open.

[THM-4015](../01-canon/theorems/THM-4015-first-kind-character-sensitive-foster-transference.md)
canonizes the first-kind theorem and the exact stopping boundary. The natural
vectors `n_i pi(e_i)` do form a strict obtuse superbase, but only for an
index-`product n_i` sublattice. Matrix-tree detects exactly the squared index
factor. This is not a nonexistence proof for a different full superbase: the
projection for `n=(1,2,3)` has an explicit mixed-coordinate full superbase
even though its natural one has index six.

The stronger conjecture that the same `sqrt(d)/2` constant holds for every
lattice survived a second exact hostile pass: `D_4` and `D_5` attain `d` in
the squared product, `A_4` gives `16/5`, and 5,728 certified random
characteristic instances in ranks four and five never exceed `d`. This is
finite evidence, not a theorem. The exact proof wall is

```text
p/2 in Vor(L) and p primitive
  does not yet control
min{||y||:y in L*, <p,y> odd}.                          (20d)
```

## 7. Incoming mathematics and the observer lesson

The incoming main-branch promotion of THM-4007 was inspected as mathematical
signal. Its lawful content is in a different problem: a next source-normal
row forces a support coordinate missed by earlier rows. There is no map from
that Jacobian jet to LRC(14). The reusable move is narrower: after a low-order
observer saturates, add the next observer only if it retains the coordinate
the earlier quotient erased. In the present session that principle manifests
concretely as prefix order, endpoint owner, divisor digit, or odd character.

The later incoming THM-4003/4004 promotions are directly in scope and were
integrated in Section 2: they turn “internal absorption” into the exact pair
of componentwise weighted minima (9a), while deleting the top-pair
support-two subbranch. Incoming THM-4010 is unrelated; reserved THM-4011/4012
are not dependencies here.

## Connection contracts

### Euclidean transference to the rank code

```text
source:      closed lattice-free projected zonotope
target:      primitive Graver row in the speed kernel
map:         exact inball -> covering radius -> dual minimum
preserved:   Euclidean coefficient mass and exact speed equation
destroyed:   parity, support, component cut, owner, phase, arrival
sidecars:    order-two character and two-component endpoint word
test:        intersect 55,459 histograms with rank-eleven star spaces
```

### Relation row to Fourier bouquet

```text
source:      signed ordered scalar relation
target:      Hermitian danger quadratic
map:         ordered partial sums -> frequency bouquet
preserved:   sign, scale, ordering, cross-order phase
destroyed:   endpoint atom and physical component
sidecar:     signed owner-address word
test:        47 ratio types x 17 residual types, all orders
```

### Lift fibre to covering code

```text
source:      one fixed parent and its first literal x7 coefficient fibre
target:      13-group positive CNF
map:         child digit word; new clock forbidden word
preserved:   every genuinely new ansatz clock, exactly
destroyed:   old clocks, global deck, divisor ownership
sidecar:     parent status, thirteen special gcd digits, endpoint stars
test:        authentic I(13,191,1) parents rather than the AP control alone
```

## Ranked continuation

1. Decide the relative order-two transference inequality (20) beyond
   first-kind lattices: prove it, find a counterexample, or extend Foster to
   the speed-weighted circuit while retaining its Smith index.
2. Solve the 78-constraint exact ellipsoid separation problem against internal
   component relations. A success must force a crossing row, not merely lower
   a scalar norm.
3. Traverse the 55,459 coefficient histograms inside rank-eleven star spaces,
   retaining component cut and owner word. Compute each component's internal
   weighted minimum before claiming the short row crosses; THM-4003's tail
   absorber and the `1+2-3` body hostile are mandatory controls.
4. Classify the divisor anchors `d|U` and `sd|tp` in the 16 scale-one residual
   shapes. Intersect with owner-labelled endpoint words; scalar divisibility
   alone gives rank, not a lonely time.
5. Run the all-order owner-labelled Fourier SDP on the 47 support-two ratios
   crossed with the 17 live residual types. Keep AP, Goddyn--Wong, and the two
   THM-4002 phase twins as mandatory hostiles.
6. Export authentic `I(13,191,1)` parents and compile their grouped `x7` CNFs
   incrementally. Seek small independently checkable UNSAT cores rather than
   treating two solver verdicts as kernel certificates.

The main gain is a change in the finite object. Before this session the
short-relation route was a huge coefficient ball. It is now a 55,459-shape
Graver atlas, a localized `l1<=49` relation bank, a separate order-two
character with an exact graphic model, and three concrete ways of retaining
the missing transverse information. The proof is not finished, but the
bottleneck is no longer unnamed.
