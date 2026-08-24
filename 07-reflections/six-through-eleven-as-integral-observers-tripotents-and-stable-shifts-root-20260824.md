# Six through eleven as integral observers, tripotents, and stable shifts

> **RESEARCH SYNTHESIS, 2026-08-24.**  Proof claims in this reflection are
> controlled by
> [THM-4000](../01-canon/theorems/THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas.md),
> its [exact companion](../04-computation/six_to_eleven_number_machinery_thm4000.py),
> and the cited canonical files.  Historical number-pattern reflections were
> used as query generators, not as truth sources.  No consequence for LRC(14),
> class-group rank, Bott periodicity, or the complex-structure problem on `S6`
> is asserted without a typed map.

## Executive answer

The productive object is not the list `6,7,8,9,10,11`.  It is the evaluation
of an integer polynomial at consecutive nodes.

```text
observer node a:       -1    0    1    2    3    4
distance from base 10: 11   10    9    8    7    6
```

For every `P in Z[X]`, the six exact samples compile `P(10)` modulo

```text
6*7*8*9*10*11=332640.
```

The central quotient is `Z[X]/(X^3-X)`, observed at `-1,0,1`.  It treats

```text
6,7,8  as  7-1,7,7+1,
9,10,11 as 10-1,10,10+1,
```

with one integral parity invoice at the two endpoints.  From that one object
come exact CRT formulas, the classification of tripotents `x^3=x`, a Farey
cusp ray, and Pell towers.

There are two other valid but different appearances of eight:

```text
arithmetic clock:       Z[X]/(X^4+1), where the base has order 8;
lattice stabilization: L |-> L direct_sum E8, preserving the discriminant form.
```

They must not be identified with each other or with Bott periodicity.  In
particular, the honest period-eight statement sends ranks `9,10,11` back to
`1,2,3`, not all three back to one.

## 1. Inheritance pass and session portfolio

| lane | nearest proved mechanism | canonical hostile | corrected near miss | least-used sidecar |
|---|---|---|---|---|
| **Anchor: integral observer** | the split fiber `X^3-X` in [THM-3555](../01-canon/theorems/THM-3555-catalan-thickening-universal-cubic-root-cover.md) | residues forget exact lifts | “the lcm is the whole modulus” | factorial Smith coordinates |
| **Niche: CRT/Farey/Pell** | the parabolic spine of [THM-3334](../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md) | modulo `10` there are six, not three, tripotents | “`-1,0,1` exhaust every modulus” | idempotent support plus corner involution |
| **Wildcard: eight-shifts** | cyclotomic clocks and `E8` orthogonal sum | `Z^8` and `E8` have different norm shells | “rank eight means the same theta series” | discriminant form and theta normalization |
| **Repo pull: scalar observers** | OCF fugacity in [THM-3380](../01-canon/theorems/THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries.md) | the kernel `4e_1-e_3` first appears at covered size nine | “one scalar retains the packing vector” | `b_333,b_335` |
| **Repo pull: min-plus** | Moon--Busch floor imported by [THM-1370](../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) | quadratic and Fibonacci fits both break | “a short sequence determines its law” | the optimizer `(a,b)` / block word |

The old reflections
[Numbers as Patterns of Three Types](numbers-as-type-patterns.md) and
[Representation Theory of the Natural Numbers](representation-theory-of-numbers.md)
contain many untyped identifications.  Their useful residue was a list of
search terms.  Every connection retained below has an object, a map, a
preserved predicate, and a hostile.

## 2. The exact octonionic `6 -> 7 -> 8` ladder

The previous
[octonion/Hopf audit](octonions-at-the-center-s6-bott-lyapunov-and-class-rank-frontiers-root-20260824.md)
already isolates the genuine geometric operation:

```text
O = R direct_sum Im(O),       dim O=8,       dim Im(O)=7,
S6 = {x in Im(O): |x|=1}.
```

Thus the maps are a norm-one hypersurface and unitalization:

```text
S(Im(O))  ->  Im(O)  ->  R direct_sum Im(O)=O,
dimension 6           7                         8.
```

The first inclusion preserves the norm and the `G2` action but forgets radial
scale.  The second preserves the imaginary representation and adds a scalar
unit.  The octonion cross product, or equivalently its positive three-form, is
the load-bearing sidecar: the same dimension ladder exists for every
seven-dimensional Euclidean vector space and is not octonionic by itself.

This exact ladder explains why `G2/SU(3)=S6`, `Spin(7)/G2=S7`, and the
octonionic Hopf geometry cluster around dimensions six through eight.  It does
not make the standard almost-complex structure on `S6` integrable.  The
unrestricted complex-`S6` question remains **OPEN**.

## 3. The integral split cubic makes both triples one object

Evaluation gives an injective finite-index map

```text
Z[X]/(X^3-X) -> Z^3,
P |-> (P(-1),P(0),P(1)).
```

Its image is exactly the triples `(u_-,u_0,u_+)` satisfying

```text
u_- == u_+ (mod 2),
```

and its Smith form is `(1,1,2)`.  Therefore

```text
P(b) == b(b+1)/2 P(1) +(1-b^2)P(0)+b(b-1)/2 P(-1)
        (mod b^3-b).                                    (A)
```

At base seven, `(A)` is

```text
P(7) == 28P(1)-48P(0)+21P(-1)  (mod 336),
```

where the separate moduli `6,7,8` reach only `168`; exact endpoint parity
supplies the second factor two.  At base ten,

```text
P(10) == 55P(1)-99P(0)+45P(-1) (mod 990),
```

and `9,10,11` are pairwise coprime, so ordinary CRT already reaches `990`.

This gives a precise reading of “nine, ten, and eleven are like one.”  They are
not three copies of one.  Relative to decimal evaluation they are the three
fibers at the observer coordinates

```text
modulus:             9    10    11
observer 10-modulus: 1     0    -1.
```

The quotient forgets every multiple of `X^3-X`; the exact sample lifts and
their parity are the sidecar needed to reconstruct the strongest modulus.

## 4. Six observers reveal a factorial carry lattice

For consecutive nodes `a,...,a+m`, Newton interpolation gives

```text
P(B) == sum_(j=0)^m binom(B-a,j) Delta^jP(a)
        (mod product_(j=0)^m(B-a-j)).                   (B)
```

The modulus is optimal over `Z[X]`, because the node polynomial itself attains
it.  The evaluation lattice has Smith diagonal

```text
0!,1!,2!,...,m!.
```

For nodes `-1,...,4` and `B=10`, `(B)` becomes

```text
P(10) == -252P(-1)+1386P(0)-3080P(1)
          +3465P(2)-1980P(3)+462P(4)  (mod 332640).
```

The individual residues modulo `6,...,11` retain only their lcm `27720`.
The exact samples retain a further factor `12`.  That factor is not mystical
and is not the cokernel order: it is the overlap debt between product and lcm,
paid by the factorial compatibility data

```text
1,1,2,6,24,120.
```

This is reusable number-theoretic machinery: whenever several modular tests
come from evaluations of one integral object, compute the evaluation lattice
before applying CRT.  Exact lifts can contain prime-power information which
the separated residues have destroyed.

## 5. Tripotents distinguish the six moduli sharply

Let `Trip(n)={x mod n:x^3=x}`.  CRT and Hensel lifting give

```text
n   Trip(n)
6   {0,1,2,3,4,5}
7   {0,1,6}
8   {0,1,3,5,7}
9   {0,1,8}
10  {0,1,4,5,6,9}
11  {0,1,10}.
```

For `n>=3`, the centered roots `0,+1,-1` exhaust precisely the odd prime
powers and `n=4`.  Thus `7,9,11` really do have the same pure trichotomy, but
`8` and `10` are decisive hostiles.  Modulo ten, the decimal centered triple
is only the diagonal part of a six-state CRT algebra.

The correct ring-theoretic coordinate is also revealing.  If `x^3=x`, then

```text
e=x^2,       e^2=e,       ex=x.
```

So `x` is an involution in the corner ring `eR`.  Only when two is invertible
does it split as a difference of orthogonal idempotents
`(x^2+x)/2-(x^2-x)/2`.  The hostile `x=3 mod 8` explains why the loose phrase
“signed idempotent” loses exactly the interesting 2-adic states.

Six has a separate universal role:

```text
gcd_{z in Z}(z^3-z)=6,
x^3==x (mod n) for every x  iff  n divides 6.
```

The factor two comes from the endpoint resultant collision; the factor three
comes from the three roots exhausting `F_3`.  Even inside one integer, the two
prime factors have different mechanisms.

## 6. Translation gives a Farey ray; iteration gives Pell towers

Define

```text
A_b=[b b-1; b+1 b],       C=[0 1; -1 2].
```

Then `det(A_b)=1` and `C A_b=A_(b+1)`.  Therefore

```text
6/7 < 7/8,        9/10 < 10/11
```

are two edges on the same parabolic ray toward the cusp one.  Gaussian
squaring sends `(b,b-1)` to

```text
(2b-1,2b(b-1),2b(b-1)+1),
```

giving `(13,84,85)` at `b=7` and `(19,180,181)` at `b=10`.  This is the exact
pre-Gaussian carrier of the Berggren spine in THM-3334.

Powering at a fixed base is a different operation.  Writing `d=b^2-1`,

```text
A_b^n=[T_n(b), (b-1)U_(n-1)(b);
       (b+1)U_(n-1)(b), T_n(b)],
T_n(b)^2-dU_(n-1)(b)^2=1.
```

At odd depth, the endpoints factor as squares:

```text
T_(2m+1)(b)-1=(b-1)(U_m+U_(m-1))^2,
T_(2m+1)(b)+1=(b+1)(U_m-U_(m-1))^2.
```

Hence the decimal packet seeds

```text
X=10,3970,1580050,...,
X-1=9A^2,       X+1=11B^2,
```

while the base-seven packet seeds `7,1351,262087,...` with square-scaled
endpoints six and eight.  Translation along the cusp and Pell depth at a fixed
base commute with neither interpretation automatically; retaining both
operations prevents a pretty but false recurrence.

## 7. Three meanings of “period eight,” only two developed here

### 7.1 Stable rank

Orthogonal sum with the even unimodular lattice `E8` gives

```text
discform(L direct_sum E8)=discform(L),
Theta_(L direct_sum E8)=Theta_L E4.
```

For the root lattices this yields

```text
A1 -> A1 direct_sum E8: rank 1 -> 9,  |disc|=2, roots 2 -> 242;
A2 -> A2 direct_sum E8: rank 2 -> 10, |disc|=3, roots 6 -> 246;
A3 -> A3 direct_sum E8: rank 3 -> 11, |disc|=4, roots 12 -> 252.
```

This is the clean stable sense in which `9,10,11` resemble smaller numbers:

| rank | centered decimal coordinate | ancestor after an `E8` shift | tripotent behavior |
|---:|---:|---:|---|
| `9` | `+1` | `1` | pure `0,+1,-1` |
| `10` | `0` | `2` | six mixed CRT states |
| `11` | `-1` | `3` | pure `0,+1,-1` |

Only nine is congruent to one modulo eight.  Ten and eleven return to two and
three.  Moreover stable discriminant form does not preserve rank, root count,
or theta series.

[MISTAKE-484](../01-canon/MISTAKES.md) records the hostile which forced this
typing: in common norm shells, `Z^8` has `112` and `1136` vectors of norms two
and four, while `E8` has `240` and `2160`.  Equal rank is not lattice equality.

### 7.2 Arithmetic order eight

The separate cyclotomic quotient

```text
Phi_8(b)=b^4+1
```

makes `b` have exact multiplicative order eight.  Every odd prime divisor is
one modulo eight, and polynomial reduction alternates four-digit base-`b`
blocks.  In particular,

```text
7^4+1=2*1201,       ord_1201(7)=8,
10^4+1=73*137,      ord_73(10)=ord_137(10)=8.
```

This is an arithmetic clock, not Bott periodicity.  Bott uses Clifford Morita
stabilization; the cyclotomic observer uses multiplicative order in a changing
finite quotient; `E8` stabilization uses orthogonal direct sum.  The common
integer eight is output, not a functor between the three categories.

## 8. The repo's real `6 -> 7 -> 8 -> 9` boundary ladder

Several seemingly unrelated tournament results align because the available
support configurations grow with the vertex budget.  This is a genuine
invariant-resolution ladder, not an identification of the integers.

| first order | exact event | mechanism |
|---:|---|---|
| `6` | `H` ceases to be spectral | two disjoint odd cycles create the first support-correlation coordinate; [THM-1780, spectral slug](../01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md) |
| `7` | total odd-cycle count ceases to be spectral; signed Redei data becomes independent | a Hamiltonian seven-cycle mixes with overlapping compound walks; [THM-500](../01-canon/theorems/THM-500-the-second-spectral-boundary-odd-cycle-count-non-spectral-from-n7.md), [THM-1966](../01-canon/theorems/THM-1966-signed-redei-count-independent-invariant-n7.md) |
| `8` | the odd-cycle conflict graph can cease to be perfect, while real-rootedness survives | an induced five-hole fits, but a claw still needs nine vertices; [THM-019](../01-canon/theorems/THM-019-omega-perfectness.md), [THM-020](../01-canon/theorems/THM-020-real-roots.md) |
| `9` | real-rootedness can fail and fugacity evaluation first has a length kernel | three disjoint triangles fit; [THM-025](../01-canon/theorems/THM-025-real-roots-disproved.md), THM-3380 |

The OCF packing types make the continuation through eleven exact:

```text
covered size   odd-part types       cycle lengths
6              3+3                  {2}
7              7                    {1}
8              3+5                  {2}
9              9,3+3+3              {1,3}
10             3+7,5+5              {2}
11             11,3+3+5             {1,3}.
```

Evaluation at fugacity two maps `e_ell` to `2^ell`.  At covered sizes nine and
eleven its primitive kernel is

```text
4e_1-e_3.
```

Four single cycles can trade for one three-cycle packing without changing the
scalar.  Covered size ten has two packing types but both have length two, so
there is no new length kernel there.  An order-ten tournament can still inherit
the covered-size-nine loss with one unused vertex.

The immediate new consequence is **PROVED** in THM-4000:

```text
(D_ham,b_333,b_335) determines the full bivariate OCF polynomial at order 11.
```

Necessity of `b_333` is witnessed.  Necessity of `b_335` among realizable
order-eleven tournaments remains **OPEN**.

## 9. Quadratic characters connect seven to eleven but reject nine and ten

For a finite field of odd order `q`, the Paley rule

```text
x -> y  iff  y-x is a nonzero square
```

is a tournament exactly when `-1` is a nonsquare, equivalently `q==3 mod 4`.
This produces the exact Paley tournaments at seven and eleven.  At nine,
`-1` is a square, so the same character produces an undirected Paley graph,
not a tournament.  There is no field of order ten, and every even-order
tournament is nonregular anyway.

At seven, [THM-211](../01-canon/theorems/THM-211-paley-design.md) identifies the
directed triangles as two cyclic Fano planes.  At eleven,
[THM-132](../01-canon/theorems/THM-132-z11-ocf-alpha1-not-monotone.md) gives the
exact OCF census for all circulant tournaments.  The preserved coordinate is
the additive translation action plus the quadratic character; arbitrary
tournament structure is lost.  The `q==3 mod 4` gate is the cheapest hostile
test and blocks any proposed `7,9,11` Paley periodicity.

## 10. Moon--Busch turns the same window into min-plus arithmetic

The strong-tournament Hamiltonian-path floor is the cited Moon--Busch formula

```text
f(n)=min_{2a+3b=n-1} 3^a 5^b.
```

Because replacing `(a,b)` by `(a-3,b+2)` multiplies the objective by
`25/27<1`, the optimizer uses the largest allowed `b` of the required parity.
Thus

```text
n        6   7   8   9   10  11
f(n)    15  25  45  75  125 225
(a,b)  1,1 0,2 2,1 1,2 0,3 2,2.
```

This is a true number machine extracted from a seemingly unrelated theorem:
the additive numerical semigroup `2a+3b=n-1` is optimized in the tropical
weight `a log 3+b log 5`.  The optimizer records parameters for an attaining
Moon construction, but the value formula does not classify all extremal
tournaments or block placements.  The repo's earlier quadratic fit and then
its Fibonacci-style replacement both failed at the first new cases.
MISTAKE-055 is therefore part of the mechanism, not historical decoration.

## 11. Concept board and connection contracts

The session ended with seven live objects:

1. the consecutive evaluation lattice;
2. the split cubic and its resultant-two invoice;
3. tripotent CRT support/involution coordinates;
4. the parabolic Farey matrix and its Pell powers;
5. cyclotomic order-eight quotients;
6. `E8` stable discriminant data; and
7. scalarization kernels in OCF and min-plus optimization.

Their contracts are deliberately different.

| source | target / map | preserved | destroyed | sidecar / decisive test |
|---|---|---|---|---|
| exact samples of `P` | residue compiler via Newton interpolation | `P(B)` modulo the node product | coefficients above degree `m` | factorial Smith coordinates; test `P=F_(a,m)` |
| `Z[X]/(X^3-X)` | centered values at `-1,0,1` | quadratic remainder | multiples of the split cubic | endpoint parity; test base seven |
| tripotent | support `e=x^2` and corner involution | the equation `x^3=x` | sign splitting at two | test `x=3 mod 8` |
| centered base | `A_b in SL_2(Z)` and its projective interval | determinant and Farey adjacency | digit/residue meaning | retain base versus depth; compare translation with powering |
| polynomial modulo `X^4+1` | order-eight finite quotient | alternating block fold | Bott/Clifford class | factor `b^4+1` and compute orders |
| lattice `L` | `L direct_sum E8` | discriminant form | rank, roots, theta series | normalized theta shells |
| OCF packing vector | evaluation at `x=2` | Hamiltonian count | length composition | kernel basis `4e_1-e_3` and `b_333,b_335` |
| order `n` plus Moon construction | min-plus lattice point `(a,b)` | extremal `H` value | classification of all extremizers | attaining construction parameters; hostile next term |

The board changed the reading of the original question.  “Number” here is not
an intrinsic personality attached to each integer.  It is the intersection of
typed operations in which the integer participates: translation distance,
residue modulus, stable rank, support budget, finite-field order, or
optimization weight.

## 12. Frontiers opened by the machinery

### 12.1 Confluent observers

[THM-4010](../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)
now closes the first three layers: `k` Hasse jets have exact kernel `(F^k)`,
give the optimal modulus `|F(B)|^k`, and define a lattice of index
`product(j!)^(k^2)`. The missing coordinate is finer than the determinant:
the repeated-factorial Smith guess first breaks at `(m,k)=(3,2)`, changing
`36,36` into `12,108`. What survives is the mod-`p` rank
`k min(m+1,p)` and the total `p`-valuation; the higher-`p^e` distribution is
still open.

This also sharpens the comparison with the Rule 30 observer. Polynomial jets
live in one fixed integral quotient, so the monic kernel law makes a new
evaluation legal without a routing choice. A moving phase action does not:
THM-4006's `943/951` hostile needs its physical basepoint even after a long
static arithmetic shadow is retained. The shared move is “compute the kernel
of observation”; the destroyed data differ—higher congruence layers here,
chronological phase ownership there.

### 12.2 Tripotent functoriality beyond square-free intuition

Classify the fibers of `x^m=x` over prime powers using support idempotents and
corner torsion units.  The tripotent case shows why raw root counts are too
coarse at two.  Cheap test: compare `m=3,5` over `2^e`, retain the square map,
and ask which centered roots survive every lift.

### 12.3 Primitive primes along the Pell ray

The odd-depth identities make `T_(2m+1)(b)-1` and `T_(2m+1)(b)+1` controlled
square multiples of `b-1,b+1`.  Factor their primitive parts and compare their
multiplicative orders with the separate cyclotomic clock.  The required
firewall is gcd/resultant data: a shared prime in two displayed integers is not
yet a shared dynamical operation.

### 12.4 Quadratic-character versus OCF resolution

For prime powers `q==3 mod 4`, compute the full packing-length vector of the
Paley tournament before evaluating at two.  The seven- and eleven-vertex data
suggest a precise question: which Gauss-sum moments determine each OCF layer,
and where does a correlation sidecar first escape the character spectrum?
Nine supplies the hostile because the orientation itself disappears.

### 12.5 The order-eleven kernel witness

Search for two eleven-vertex tournaments with equal `D_ham` and `b_333` but
different `b_335`.  This is finite and decisive.  A positive pair proves the
third coordinate in THM-4000 is necessary; exhaustive failure would suggest a
realizability relation not visible in the free packing lattice.

### 12.6 Stable lattice arithmetic

Within a fixed discriminant form, compare genera, theta series, root systems,
and automorphism groups under repeated `E8` stabilization.  This develops the
valid `r -> r+8` number theory while keeping the invariant that actually
survives.  It does not by itself produce imaginary-quadratic class-group
generators; that challenge still requires explicit forms and relation tests.

## Final judgement

The integers six through eleven do not form one secret algebra.  They sit at
the crossing of several exact operations:

```text
6,7,8:    sphere / imaginary space / octonion;
          one centered base packet;
9,10,11: one centered decimal packet;
          ranks 1,2,3 after an E8 shift;
6..11:   one optimal consecutive-sample congruence;
7,9,11:  pure odd-primary tripotent trichotomy;
7,11:    quadratic-character tournaments;
8:       cyclotomic order and stable-rank shift;
9,11:    first OCF length-scalarization kernels.
```

The frontier lesson is more useful than a slogan: when integers appear to
repeat, ask for the operation producing the repetition, compute its kernel,
and retain the smallest sidecar that restores what the operation forgets.
